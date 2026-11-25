import unittest

from balance import MODELS, balance


def example(failures):
    return [dict(model=m, donor=str(i), nf=1, hcm=0 if (m, i) in failures else -1)
            for m in MODELS for i in range(9)]


class BalanceTest(unittest.TestCase):
    def test_same_donor_missed(self):
        self.assertTrue(balance(example({(m, 1) for m in MODELS}))["target_met"])

    def test_equal_counts_can_hide_disagreement(self):
        result = balance(example({("v1", 1), ("v2", 2)}))
        self.assertEqual(result["weakest_model_count"], 8)
        self.assertEqual(result["joint_count"], 7)
        self.assertFalse(result["target_met"])

    def test_missing_row(self):
        with self.assertRaises(ValueError):
            balance(example(set())[:-1])
