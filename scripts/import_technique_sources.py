#!/usr/bin/env python3
"""Import HoDoKu factual puzzle states or compare them with the frozen examples."""
import argparse
import hashlib
import json
from pathlib import Path


def records(raw):
    result = []
    for line_no, line in enumerate(raw.decode("utf-8-sig").splitlines(), 1):
        if not line.startswith(":"):
            continue
        fields = line.split(":")
        if "-x" in fields[1]:
            continue
        values, given, placed = [], [], False
        for c in fields[3]:
            if c == "+":
                placed = True
                continue
            digit = int(c) if c in "123456789" else 0
            values.append(digit)
            given.append(str(digit) if digit and not placed else ".")
            placed = False
        if len(values) != 81:
            raise ValueError(f"Invalid board size at source line {line_no}")
        domains = []
        for cell, value in enumerate(values):
            mask = 1 << (value - 1) if value else 511
            if not value:
                for other, digit in enumerate(values):
                    if digit and (
                        cell // 9 == other // 9
                        or cell % 9 == other % 9
                        or (cell // 27, cell % 9 // 3) == (other // 27, other % 9 // 3)
                    ):
                        mask &= ~(1 << (digit - 1))
            domains.append(mask)
        for token in fields[4].split():
            if len(token) != 3 or any(d not in "123456789" for d in token):
                raise ValueError(f"Invalid candidate token at source line {line_no}")
            d, r, c = map(int, token)
            domains[(r - 1) * 9 + c - 1] &= ~(1 << (d - 1))
        result.append(
            {
                "id": f"hodoku-reglib-1.3-line-{line_no}",
                "source_technique": fields[1],
                "puzzle": "".join(given),
                "source_url": "https://hodoku.sourceforge.net/libraries/reglib-1.3.txt",
                "source_index_1based": line_no,
                "source_sha256": hashlib.sha256(raw).hexdigest(),
                "attribution": "Bernhard Hobiger, HoDoKu regression library; individual puzzle author unspecified",
                "state": {"version": 1, "values": values, "domains": domains},
            }
        )
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source")
    parser.add_argument("--output")
    parser.add_argument("--check")
    parser.add_argument(
        "--include-constructed",
        action="store_true",
        help="Include additional checked example sources from the frozen corpus",
    )
    args = parser.parse_args()
    data = records(Path(args.source).read_bytes())
    if args.include_constructed:
        from verify_technique_examples import from_givens

        corpus = json.loads(
            (
                Path(__file__).resolve().parent.parent
                / "tests/fixtures/technique_examples.json"
            ).read_text()
        )
        for example in corpus["examples"]:
            puzzle = corpus["puzzles"][example["puzzle"]]
            if "hodoku.sourceforge.net" in puzzle["source_url"]:
                continue
            data.append(
                {
                    "id": puzzle["id"],
                    "puzzle": puzzle["givens"],
                    "state": example["state"],
                    "source_url": puzzle["source_url"],
                    "attribution": puzzle["attribution"],
                }
            )
    if args.output:
        Path(args.output).write_text(json.dumps(data, indent=2) + "\n")
    if args.check:
        corpus = json.loads(Path(args.check).read_text())
        by_line = {r["source_index_1based"]: r for r in data}
        checked = 0
        for example in corpus["examples"]:
            puzzle = corpus["puzzles"][example["puzzle"]]
            if "hodoku.sourceforge.net" not in puzzle["source_url"]:
                continue
            source = by_line[puzzle["source_index_1based"]]
            if (
                puzzle["source_sha256"] != source["source_sha256"]
                or puzzle["givens"] != source["puzzle"]
                or example["state"] != source["state"]
            ):
                raise ValueError("Frozen example differs from its source")
            checked += 1
        if checked != 47:
            raise ValueError(f"Expected 47 published-source examples, got {checked}")
        print(
            f"Checked {checked} example source states against the original source bytes"
        )
    else:
        print(f"Imported {len(data)} source states")
