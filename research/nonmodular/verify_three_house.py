"""Independent permutation and Hall checks for the three-house research fixture.

Standard library only. This is a finite research check, not a production solver.
Candidate domains are premises; uniqueness is neither assumed nor required.
"""
import itertools
import json
from pathlib import Path

if not __debug__:
    raise SystemExit("Run without -O: this research checker uses assertions.")

fixture = json.loads(Path(__file__).with_name("three-house-847.json").read_text())
assert fixture["schema_version"] == 1
assert fixture["size"] == 9 and fixture["box_size"] == 3
masks = fixture["domains"]
assert len(masks) == 81 and all(0 < m <= 511 for m in masks)
target = fixture["target"]
assert target == {"cell": 3, "digit": 1, "inference": "elimination"}

rows = [[9 * r + c for c in range(9)] for r in range(9)]
columns = [[9 * r + c for r in range(9)] for c in range(9)]
boxes = [
    [9 * (3 * br + r) + 3 * bc + c for r in range(3) for c in range(3)]
    for br in range(3) for bc in range(3)
]
houses = rows + columns + boxes


def domain(cell):
    return {d for d in range(1, 10) if masks[cell] & (1 << (d - 1))}


def check_completion(board):
    assert len(board) == 81
    assert all(board[c] in domain(c) for c in range(81))
    assert all({board[c] for c in h} == set(range(1, 10)) for h in houses)


completion = list(map(int, fixture["completion"]))
check_completion(completion)
# An explicit second completion confirms that uniqueness rules are inapplicable.
alternative = completion.copy()
alternative[54:63], alternative[72:81] = completion[72:81], completion[54:63]
check_completion(alternative)
assert alternative != completion

assumed = masks.copy()
assert assumed[3] & 1
assumed[3] = 1  # r1c4=1
permutations = list(itertools.permutations(range(1, 10)))
valid = {}
for name, house in [("row", rows[0]), ("column", columns[0]), ("box", boxes[0])]:
    valid[name] = [
        p for p in permutations
        if all(assumed[c] & (1 << (d - 1)) for c, d in zip(house, p))
    ]

row_boundary = {p[:3] for p in valid["row"]}
column_boundary = {p[:3] for p in valid["column"]}
row_box = [p for p in valid["box"] if p[:3] in row_boundary]
column_box = [p for p in valid["box"] if (p[0], p[3], p[6]) in column_boundary]
row_column = {p[0] for p in valid["row"]} & {p[0] for p in valid["column"]}
joint = [p for p in row_box if (p[0], p[3], p[6]) in column_boundary]
counts = {name: len(values) for name, values in valid.items()}
assert counts == {"row": 28, "column": 300, "box": 180}
assert row_box and column_box and row_column
assert not joint

# A compact proof independent of the permutation enumeration.
# Under the target assumption, these row cells must occupy digits 2 and 3.
assert domain(1) - {1} == domain(2) - {1} == {2, 3}
# They share box 1 with r2c1, which must therefore be 4.
assert domain(9) - {2, 3} == {4}
# The same row and box then force r1c1=8.
assert domain(0) - {1, 2, 3, 4} == {8}
# Column 1 loses 4 and 8. Its r3c1 also loses 2 and 3 in box 1.
hall_cells = [18, 36, 45, 54, 72]
available = {c: domain(c) - {4, 8} for c in hall_cells}
available[18] -= {2, 3}
union = set().union(*available.values())
assert union == {3, 5, 7, 9}
assert len(union) < len(hall_cells)

print("Two distinct full completions respect all recorded domains.")
print("Under r1c4=1, individual-house assignment counts:", counts)
print("All three pairs remain feasible; compatible three-house assignments:", len(joint))
print("Independent Hall proof: five column cells have only digits", sorted(union))
print("Verified elimination: r1c4 != 1; no uniqueness premise.")
