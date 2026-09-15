"""Regression tests for the production coverage gate; no Rust compiler needed."""

import importlib.util
from pathlib import Path
import tempfile
import unittest

spec = importlib.util.spec_from_file_location("sudoku_coverage", Path(__file__).with_name("coverage.py"))
coverage = importlib.util.module_from_spec(spec)
spec.loader.exec_module(coverage)


class CoverageGateTests(unittest.TestCase):
    def setUp(self):
        self.directory = tempfile.TemporaryDirectory()
        self.addCleanup(self.directory.cleanup)
        self.root = Path(self.directory.name)
        (self.root / "Cargo.toml").write_text('[package]\nname = "coverage-fixture"\n')
        self.paths = sorted(set(coverage.FOUNDATIONS + coverage.SOLVER + coverage.GENERATION + coverage.ADVANCED))
        for path in self.paths:
            self.source(path, "fn example() {\n    work();\n}\n")
        self.lcov = self.root / "coverage.info"

    def source(self, relative, contents):
        path = self.root / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(contents)
        return path

    def record(self, relative, hits=(1, 1), functions=(1,)):
        return "\n".join([
            f"SF:{self.root / relative}",
            *(f"FN:{line},function_{line}" for line in functions),
            *(f"DA:{index},{hit}" for index, hit in enumerate(hits, 1)),
            "end_of_record",
        ]) + "\n"

    def report(self, extra="", omit=(), overrides=None):
        overrides = overrides or {}
        self.lcov.write_text("".join(
            self.record(path, overrides.get(path, (1, 1)))
            for path in self.paths if path not in omit
        ) + extra)
        return coverage.evaluate(self.root, self.lcov)

    def test_duplicate_records_union_hits_without_inflating_lines(self):
        report = self.report(self.record(self.paths[0], (0, 1)), overrides={self.paths[0]: (1, 0)})
        self.assertTrue(report["passed"])
        self.assertEqual(report["checks"]["production"]["measured"], len(self.paths) * 2)
        self.assertEqual(report["checks"]["production"]["covered"], len(self.paths) * 2)

    def test_zero_hits_stay_in_denominator_and_engine_minimum_is_enforced(self):
        report = self.report(overrides={coverage.ADVANCED[0]: (1, 0)})
        self.assertFalse(report["passed"])
        self.assertFalse(report["checks"]["aic_engine"]["passed"])
        self.assertEqual(report["checks"]["production"]["measured"], len(self.paths) * 2)

    def test_group_thresholds_use_all_members(self):
        for group, paths in (("foundations", coverage.FOUNDATIONS), ("solver", coverage.SOLVER),
                             ("generation", coverage.GENERATION), ("advanced", coverage.ADVANCED)):
            with self.subTest(group=group):
                report = self.report(overrides={path: (1, 0) for path in paths})
                self.assertFalse(report["checks"][group]["passed"])
                self.assertEqual(report["checks"][group]["measured"], 2 * len(paths))
                self.assertEqual(report["checks"][group]["percent"], 50)

    def test_missing_production_file_and_missing_function_fail_closed(self):
        report = self.report(omit=[self.paths[0]])
        self.assertFalse(report["passed"])
        self.assertEqual(report["missing_production_functions"][0]["function_lines"], [1])
        self.source(self.paths[0], "fn first() {}\nfn omitted() {}\n")
        report = self.report()
        self.assertFalse(report["passed"])
        self.assertEqual(report["missing_production_functions"][0]["function_lines"], [2])

    def test_new_unmapped_production_module_cannot_be_silently_excluded(self):
        self.source("src/new_engine.rs", "fn new_behavior() {}\n")
        report = self.report()
        self.assertFalse(report["passed"])
        self.assertIn("src/new_engine.rs", [item["file"] for item in report["missing_production_functions"]])

    def test_only_test_sources_and_external_dependencies_are_excluded(self):
        excluded = ["src/grid_tests.rs", "src/solver/tests.rs", "src/tests/fixtures.rs",
                    "tests/public_api.rs", "vendor/dependency.rs"]
        for path in excluded:
            self.source(path, "fn uncovered_test() {}")
        report = self.report("".join(self.record(path, (0, 0)) for path in excluded))
        self.assertTrue(report["passed"])
        self.assertEqual(report["checks"]["production"]["measured"], len(self.paths) * 2)

    def test_declarations_comments_strings_and_trait_signatures_need_no_mapping(self):
        self.source("src/declarations.rs", '''// fn example() {}
/* fn example() {} */
const EXAMPLE: &str = r#"fn example() {}"#;
pub struct Model { pub value: i32 }
trait Read { fn read(&self) -> i32; }
''')
        report = self.report()
        self.assertTrue(report["passed"])
        self.assertIn("src/declarations.rs", report["sources_without_executable_mappings"])

    def test_inline_tests_are_rejected_instead_of_counted_as_production(self):
        self.source(self.paths[0], "fn example() {}\n#[cfg(test)]\nmod tests { fn test() {} }\n")
        with self.assertRaisesRegex(ValueError, "Move inline tests"):
            self.report()

    def test_empty_and_invalid_reports_fail(self):
        self.lcov.write_text("")
        self.assertFalse(coverage.evaluate(self.root, self.lcov)["passed"])
        self.lcov.write_text(self.record(self.paths[0], (1, -1)))
        with self.assertRaisesRegex(ValueError, "Invalid DA"):
            coverage.evaluate(self.root, self.lcov)


if __name__ == "__main__":
    unittest.main()
