"""Local inference on records from a checkpoint-compatible tokenizer."""

import argparse
import json
import tempfile
from pathlib import Path

from .cli import sha


def load_model(kind, folder, device):
    import torch

    if kind in {"v1", "v2"}:
        from transformers import BertModel

        model, info = BertModel.from_pretrained(
            folder,
            local_files_only=True,
            add_pooling_layer=False,
            output_loading_info=True,
        )
        if info["missing_keys"] or info["mismatched_keys"]:
            raise ValueError("Incomplete Geneformer encoder")
        maximum = 2048 if kind == "v1" else 4096
        if model.config.max_position_embeddings != maximum:
            raise ValueError(
                "Checkpoint context does not match the selected Geneformer version"
            )

        def forward(tokens, values):
            hidden = model(tokens, output_hidden_states=True).hidden_states[-2][0]
            if kind == "v2":
                hidden = hidden[1:-1]
            return hidden.mean(0)
    else:
        from scgpt.model import TransformerModel

        cfg = json.loads((folder / "args.json").read_text())
        vocab = json.loads((folder / "vocab.json").read_text())
        model = TransformerModel(
            ntoken=len(vocab),
            d_model=cfg["embsize"],
            nhead=cfg["nheads"],
            d_hid=cfg["d_hid"],
            nlayers=cfg["nlayers"],
            nlayers_cls=cfg["n_layers_cls"],
            n_cls=1,
            vocab=vocab,
            dropout=cfg["dropout"],
            pad_token=cfg["pad_token"],
            pad_value=cfg["pad_value"],
            do_mvc=True,
            use_fast_transformer=False,
            pre_norm=False,
        )
        weights = torch.load(
            folder / "best_model.pt", map_location="cpu", weights_only=True
        )
        weights = {k.replace("Wqkv.", "in_proj_"): v for k, v in weights.items()}
        keys = [
            k
            for k in model.state_dict()
            if k.startswith(("encoder.", "value_encoder.", "transformer_encoder."))
        ]
        model.load_state_dict({k: weights[k] for k in keys}, strict=False)

        def forward(tokens, values):
            return model._encode(tokens, values, tokens.eq(vocab[cfg["pad_token"]]))[
                0, 0
            ]

    model.to(device).eval()
    return forward


def score_record(r, kind, encode):
    ids = r["tokens"]
    genes = r["genes"]
    values = r.get("values", [])
    if (
        not ids
        or len(ids) != len(genes)
        or len(set(ids)) != len(ids)
        or len(set(genes)) != len(genes)
    ):
        raise ValueError(
            "Expected unique, unpadded gene tokens with parallel gene identities"
        )
    if kind == "v2" and (
        genes[0] != "<cls>" or genes[-1] != "<eos>" or len(ids) > 4096
    ):
        raise ValueError("V2 requires native CLS and EOS and at most 4096 tokens")
    if kind == "v1" and (any(g.startswith("<") for g in genes) or len(ids) > 2048):
        raise ValueError("V1 requires ranked genes only, maximum 2048")
    if kind == "scgpt" and (len(values) != len(ids) or genes[0] != "<cls>"):
        raise ValueError("scGPT requires paired expression values and initial CLS")
    result = {k: r[k] for k in ("cell", "donor", "role", "state")}
    result["original"] = encode(ids, values)
    if r["role"] == "reference":
        return result
    target = r["target_id"]
    controls = r["control_ids"]
    result.update(target_id=target, control_ids=controls)
    if not controls or target in controls or len(set(controls)) != len(controls):
        raise ValueError("Choose distinct matched controls")
    members = r.get("target_genes", [target])
    nulls = r.get("control_sets", [[c] for c in controls])
    if len(nulls) != len(controls):
        raise ValueError("Control sets must align with control IDs")
    for group in [members, *nulls]:
        if (
            not isinstance(group, list)
            or not group
            or len(set(group)) != len(group)
            or any(not isinstance(g, str) or not g or g.startswith("<") for g in group)
        ):
            raise ValueError("Expected unique gene identities for each edit")
    if any(len(group) != len(members) or set(group) & set(members) for group in nulls):
        raise ValueError(
            "Match edit cardinality and exclude target genes from controls"
        )
    result.update(target_genes=members, control_sets=nulls)
    present = set(genes)
    if all(set(group) <= present for group in [members, *nulls]) and sum(
        not g.startswith("<") for g in genes
    ) > len(members):

        def remove(group):
            keep = [i for i, g in enumerate(genes) if g not in group]
            return encode(
                [ids[i] for i in keep], [values[i] for i in keep] if values else []
            )

        result.update(
            target=remove(members), controls=[remove(group) for group in nulls]
        )
    return result


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("model", choices=["v1", "v2", "scgpt"])
    p.add_argument("checkpoint", type=Path)
    p.add_argument("tokens", type=Path)
    p.add_argument("output", type=Path, help="New run directory")
    p.add_argument("--device", default="cpu")
    a = p.parse_args()
    from importlib.metadata import version

    import torch

    if a.output.exists():
        p.error("Choose a new output directory")
    a.output.parent.mkdir(parents=True, exist_ok=True)
    infer = load_model(a.model, a.checkpoint, a.device)

    @torch.inference_mode()
    def encode(ids, values):
        x = torch.tensor([ids], dtype=torch.long, device=a.device)
        v = torch.tensor([values], dtype=torch.float32, device=a.device)
        return infer(x, v).float().cpu().tolist()

    with tempfile.TemporaryDirectory(dir=a.output.parent) as tmp:
        stage = Path(tmp) / "run"
        stage.mkdir()
        path = stage / "embeddings.jsonl"
        replay_errors = []
        with a.tokens.open() as src, path.open("w") as dest:
            for line in src:
                if line.strip():
                    record = json.loads(line)
                    result = score_record(record, a.model, encode)
                    if len(replay_errors) < 3:
                        replay = encode(record["tokens"], record.get("values", []))
                        replay_errors.append(
                            max(abs(x - y) for x, y in zip(replay, result["original"]))
                        )
                    dest.write(json.dumps(result, allow_nan=False) + "\n")
        if not replay_errors:
            raise ValueError("Empty native token file")
        files = [
            p
            for p in a.checkpoint.iterdir()
            if p.is_file() and p.suffix in {".json", ".pt", ".bin", ".safetensors"}
        ]
        receipt = {
            "model": a.model,
            "device": a.device,
            "noop_replay_cells": len(replay_errors),
            "noop_max_absolute_coordinate_error": max(replay_errors),
            "cuda_runtime": torch.version.cuda,
            "device_name": torch.cuda.get_device_name(a.device)
            if torch.device(a.device).type == "cuda"
            else a.device,
            "torch": str(torch.__version__),
            "adapter_version": version(
                "scgpt" if a.model == "scgpt" else "transformers"
            ),
            "checkpoint": {p.name: sha(p) for p in files},
            "input_sha256": sha(a.tokens),
            "output_sha256": sha(path),
            "adapter_sha256": sha(Path(__file__)),
        }
        (stage / "model.json").write_text(json.dumps(receipt, indent=2) + "\n")
        stage.rename(a.output)
