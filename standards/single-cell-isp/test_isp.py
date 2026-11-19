import json
import shutil
import tempfile
import unittest
from pathlib import Path

from isp import HERE, fingerprint, inspect, score


class Workflow(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.run = Path(self.tmp.name) / 'study'
        shutil.copytree(HERE / 'example', self.run)
        (self.run / 'lock.json').write_text(json.dumps(fingerprint(self.run)))
        self.shifts = self.run / 'shifts.tsv'

    def test_donor_intersection(self):
        out = score(self.run, self.shifts)
        self.assertEqual(out['passing_in_all_models'], ['X'])
        self.assertEqual(out['models']['v1']['passing'], 2)
        self.assertEqual(out['models']['scgpt']['passing'], 1)

    def test_incomplete_or_wrong_controls(self):
        original = self.shifts.read_text()
        for wrong in ['\n'.join(original.splitlines()[:-1]), original.replace('G2,G3','G1,G2')]:
            self.shifts.write_text(wrong)
            with self.assertRaises(ValueError):
                score(self.run, self.shifts)

    def test_split_leakage(self):
        with (self.run / 'donors.tsv').open('a') as f:
            f.write('X\ttrain\n')
        with self.assertRaises(ValueError):
            inspect(self.run)

    def test_modified_design(self):
        with (self.run / 'study.json').open('a') as f:
            f.write('\n')
        with self.assertRaises(ValueError):
            score(self.run, self.shifts)


class CisWorkflow(unittest.TestCase):
    def test_registered_locus_identity_and_missingness(self):
        from cis import evaluate, hashes
        with tempfile.TemporaryDirectory() as tmp:
            folder = Path(tmp) / 'study'
            shutil.copytree(HERE / 'example', folder)
            (folder / 'cis-lock.json').write_text(json.dumps(hashes(folder)))
            path = folder / 'alleles.tsv'
            result = evaluate(folder, path)
            self.assertEqual((result['registered_loci'], result['measured_loci'], result['same_sign']), (2, 1, 1))
            path.write_text(path.read_text().replace('mock2', 'wrong'))
            with self.assertRaises(ValueError):
                evaluate(folder, path)
