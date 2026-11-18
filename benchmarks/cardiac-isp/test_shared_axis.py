import unittest
from shared_axis import shared_axis


class CommonDirection(unittest.TestCase):
    def test_models_must_agree_on_which_axis(self):
        rows = [dict(model=m, donor="A", nf=x, hcm=y) for m, x, y in
                [("v1", 1, -1), ("v2", 1, 1), ("scgpt", -1, -1)]]
        self.assertEqual(shared_axis(rows, 1)["either"], [])
        rows[-1]["nf"] = 1
        self.assertEqual(shared_axis(rows, 1)["either"], ["A"])
        self.assertEqual(shared_axis(rows, 1)["both"], [])
        with self.assertRaises(ValueError):
            shared_axis(rows[:-1], 1)
