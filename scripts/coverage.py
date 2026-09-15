#!/usr/bin/env python3
"""Enforce coverage of authored Rust using unique LCOV source-line entries.

LLVM's aggregate line totals are deliberately unused: a source line can occur in
several binaries or generic instantiations. It counts once here, with any hit
marking it covered. The report must also contain every production function.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import re
import sys
import tomllib


FOUNDATIONS = (
    "src/bitset.rs", "src/cell.rs", "src/constraint.rs", "src/grid.rs",
    "src/hash.rs", "src/position.rs", "src/puzzle_id.rs",
    "src/solver/fabric.rs", "src/solver/types.rs", "src/solver/explain.rs",
)
SOLVER = ("src/solver/mod.rs", "src/solver/basic.rs", "src/solver/backtrack.rs")
GENERATION = ("src/generator.rs", "src/diversity.rs")
ADVANCED = (
    "src/solver/aic_engine.rs", "src/solver/als_engine.rs",
    "src/solver/arithmetic.rs", "src/solver/fish_engine.rs", "src/solver/uniqueness.rs",
)
THRESHOLDS = {
    "production": (("",), 75.0),
    "foundations": (FOUNDATIONS, 90.0),
    "solver": (SOLVER, 90.0),
    "generation": (GENERATION, 80.0),
    "advanced": (ADVANCED, 75.0),
    **{Path(path).stem: ((path,), 60.0) for path in ADVANCED},
}


def is_test_source(path: Path) -> bool:
    return "tests" in path.parts or path.name == "tests.rs" or path.name.endswith("_tests.rs")


def production_sources(root: Path) -> set[Path]:
    manifest = tomllib.loads((root / "Cargo.toml").read_text())
    workspace = manifest.get("workspace", {})
    packages = {root} if "package" in manifest else set()
    for pattern in workspace.get("members", []):
        packages.update(path for path in root.glob(pattern) if (path / "Cargo.toml").is_file())
    for pattern in workspace.get("exclude", []):
        packages.difference_update(root.glob(pattern))
    return {
        path.resolve()
        for package in packages
        for path in (package / "src").rglob("*.rs")
        if not is_test_source(path.relative_to(root))
    }


def function_lines(path: Path) -> set[int]:
    """Find authored function definitions, leaving trait declarations out.

    Mask comments and string literals while preserving line numbers. The source
    inventory is independent of LCOV so an omitted file or function fails closed.
    """
    source = path.read_text()
    token = re.compile(r'//[^\n]*|/\*.*?\*/|r(?P<hashes>\#*)".*?"(?P=hashes)|"(?:\\.|[^"\\])*"', re.S)
    code = token.sub(lambda match: re.sub(r"[^\n]", " ", match.group()), source)
    if re.search(r"#\s*\[\s*cfg\s*\(\s*test\s*\)\s*\]\s*mod\s+\w+\s*\{", code):
        raise ValueError(f"Move inline tests into a separate test source: {path}")
    return {
        code.count("\n", 0, match.start()) + 1
        for match in re.finditer(r"\bfn\s+(?:r#)?[A-Za-z_]\w*\b[^;{]*\{", code)
    }


def read_lcov(path: Path, root: Path, sources: set[Path]) -> tuple[dict, dict]:
    lines: dict[Path, dict[int, bool]] = {}
    functions: dict[Path, set[int]] = {}
    current = None
    for number, raw in enumerate(path.read_text().splitlines(), 1):
        if raw.startswith("SF:"):
            candidate = Path(raw[3:])
            current = (candidate if candidate.is_absolute() else root / candidate).resolve()
        elif raw == "end_of_record":
            current = None
        elif current in sources and raw.startswith("DA:"):
            parts = raw[3:].split(",")
            try:
                line, hits = int(parts[0]), int(parts[1])
                if line < 1 or hits < 0:
                    raise ValueError()
            except (ValueError, IndexError) as error:
                raise ValueError(f"Invalid DA entry at {path}:{number}: {raw}") from error
            previous = lines.setdefault(current, {}).get(line, False)
            lines[current][line] = previous or hits > 0
        elif current in sources and raw.startswith("FN:"):
            try:
                start = int(raw[3:].split(",", 1)[0])
                if start < 1:
                    raise ValueError()
            except ValueError as error:
                raise ValueError(f"Invalid FN entry at {path}:{number}: {raw}") from error
            functions.setdefault(current, set()).add(start)
    return lines, functions


def totals(files: dict) -> dict:
    covered = sum(item["covered"] for item in files.values())
    measured = sum(item["measured"] for item in files.values())
    return {
        "covered": covered,
        "measured": measured,
        "percent": 100.0 * covered / measured if measured else None,
    }


def evaluate(root: Path, lcov: Path) -> dict:
    root = root.resolve()
    sources = production_sources(root)
    if not sources:
        raise ValueError("No workspace production Rust source files found")
    lines, functions = read_lcov(lcov, root, sources)
    files = {}
    missing = []
    declarative = []
    for source in sorted(sources):
        relative = source.relative_to(root).as_posix()
        expected = function_lines(source)
        missing_functions = sorted(expected - functions.get(source, set()))
        if missing_functions:
            missing.append({"file": relative, "function_lines": missing_functions})
        entries = lines.get(source, {})
        if expected and not entries and not missing_functions:
            missing.append({"file": relative, "reason": "No executable line mappings"})
        if not entries:
            declarative.append(relative)
            continue
        files[relative] = {
            "covered": sum(entries.values()),
            "measured": len(entries),
            "percent": 100.0 * sum(entries.values()) / len(entries),
            "uncovered_lines": sorted(line for line, hit in entries.items() if not hit),
        }
    checks = {}
    for name, (paths, minimum) in THRESHOLDS.items():
        aggregate = totals({path: data for path, data in files.items() if "" in paths or path in paths})
        checks[name] = {
            **aggregate,
            "minimum_percent": minimum,
            "passed": aggregate["percent"] is not None and aggregate["percent"] >= minimum,
        }
    return {
        "schema_version": 1,
        "metric": "Unique production LCOV DA source lines; covered if any execution hits the line",
        "passed": not missing and all(check["passed"] for check in checks.values()),
        "checks": checks,
        "files": files,
        "missing_production_functions": missing,
        "sources_without_executable_mappings": declarative,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--lcov", type=Path, default=Path("target/coverage/lcov.info"))
    parser.add_argument("--output", type=Path, default=Path("target/coverage/summary-production.json"))
    args = parser.parse_args()
    try:
        report = evaluate(args.root, args.lcov)
    except (OSError, ValueError) as error:
        report = {"schema_version": 1, "passed": False, "error": str(error)}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + "\n")
    if "error" in report:
        print(report["error"], file=sys.stderr)
    for name, check in report.get("checks", {}).items():
        percent = f'{check["percent"]:.2f}%' if check["percent"] is not None else "unmeasured"
        print(f'{name}: {check["covered"]}/{check["measured"]} = {percent}; minimum {check["minimum_percent"]:g}% [{"PASS" if check["passed"] else "FAIL"}]')
    for missing in report.get("missing_production_functions", []):
        print(f"Missing production coverage mappings: {missing}", file=sys.stderr)
    return 0 if report["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
