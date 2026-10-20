"""Negative controls for independent named-pattern replay."""

import copy
import json
from pathlib import Path
import unittest
from verify_named_patterns import verify_named
from verify_technique_examples import verify


class NamedPatternTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.corpus = json.loads(
            (
                Path(__file__).resolve().parent.parent
                / "tests/fixtures/technique_examples.json"
            ).read_text()
        )

    def replay(self, e, p):
        cert = e["certificate"]["Certified"]
        verify_named(
            p,
            e["state"],
            e["hint"],
            self.corpus["puzzles"][e["puzzle"]]["givens"],
            cert["deduction"],
        )

    def test_all_48_patterns_and_deductions(self):
        verify(self.corpus)

    def test_a_sound_endpoint_does_not_justify_a_missing_named_proof(self):
        c = copy.deepcopy(self.corpus)
        del c["examples"][0]["certificate"]["Certified"]["named"]
        with self.assertRaises(ValueError):
            verify(c)

    def test_binding_and_structural_mutations(self):
        for e in self.corpus["examples"]:
            original = e["certificate"]["Certified"]["named"]
            with self.subTest(method=e["hint"]["technique"]):
                self.replay(e, original)
                for field, value in [
                    ("version", 2),
                    ("state_hash", ""),
                    ("technique", "Unknown"),
                    ("conclusion", "Unique"),
                    ("unchecked", True),
                ]:
                    p = copy.deepcopy(original)
                    p[field] = value
                    with self.assertRaises(ValueError):
                        self.replay(e, p)
                p = copy.deepcopy(original)
                w = p["witness"]
                family = w if isinstance(w, str) else next(iter(w))
                a = {} if isinstance(w, str) else w[family]
                if family == "Subset":
                    a["cells"] = []
                elif family == "Fish":
                    a["parts"][0]["covers"] = []
                elif family == "Wing":
                    a["cells"][0] = a["cells"][1]
                elif family == "Als":
                    a["sets"][0]["digits"] = []
                elif family == "Chain":
                    a["nodes"][0] = 729
                elif family == "Trade":
                    a["cells"][0] = 81
                elif family == "Bug":
                    p["uniqueness"]["proof"]["nodes"] = []
                elif family == "Aligned":
                    a["exclusions"] = []
                elif family == "Forcing":
                    a["branches"][0]["steps"] = []
                elif family == "Native":
                    native = a["proof"][next(iter(a["proof"]))]
                    native["state_hash"] = ""
                elif family == "Search":
                    a["proof"]["nodes"] = []
                else:
                    self.fail("Unhandled family")
                with self.assertRaises((ValueError, IndexError)):
                    self.replay(e, p)

    def test_relabeling_a_valid_endpoint_cannot_bypass_the_pattern_rule(self):
        e = copy.deepcopy(self.corpus["examples"][0])
        p = copy.deepcopy(e["certificate"]["Certified"]["named"])
        e["hint"]["technique"] = "XYWing"
        p["technique"] = "XYWing"
        with self.assertRaises(ValueError):
            self.replay(e, p)


if __name__ == "__main__":
    unittest.main()
