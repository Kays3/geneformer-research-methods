import tempfile
import unittest
from pathlib import Path

from score import cell_scores, cis_scores

HERE = Path(__file__).resolve().parent


class Scores(unittest.TestCase):
    def test_unequal_cells_equal_donors(self):
        result = cell_scores(HERE / "cohort.tsv", HERE / "cell.tsv")
        self.assertEqual((result["cells"], result["donors"], result["joint"]), (3, 2, 1))
        self.assertAlmostEqual(result["equal_donor_nf"], 0.125)
        self.assertAlmostEqual(result["equal_donor_hcm"], 0.025)
        self.assertEqual(cell_scores(HERE / "cohort.tsv", HERE / "cell.tsv", 0.2)["joint"], 0)

    def test_pairing_failures(self):
        lines = (HERE / "cell.tsv").read_text().splitlines()
        for bad in [lines[:-1], lines + lines[-1:], [*lines[:-1], lines[-1].replace("B", "A")],
                    [*lines[:-1], lines[-1].replace("0.2", "nan")]]:
            with tempfile.TemporaryDirectory() as tmp:
                path = Path(tmp) / "bad.tsv"
                path.write_text("\n".join(bad) + "\n")
                with self.assertRaises(ValueError):
                    cell_scores(HERE / "cohort.tsv", path)

    def test_loci_and_missingness(self):
        result = cis_scores(HERE / "alleles.tsv")
        self.assertEqual((result["registered_loci"], result["directional_loci"], result["same_sign"]), (3, 2, 1))
        text = (HERE / "alleles.tsv").read_text()
        for bad in [text.replace("L2", "L1"), text.replace("NA\tNA", "0.2\tNA")]:
            with tempfile.TemporaryDirectory() as tmp:
                path = Path(tmp) / "bad.tsv"
                path.write_text(bad)
                with self.assertRaises(ValueError):
                    cis_scores(path)


if __name__ == "__main__":
    unittest.main()
