"""Combined exact acceptance of the full-old-engine-stall benchmark.

Uses only the standard library. Search, numerical optimization, the Rust engine
and historical-binary execution are separate from these certificate checks.
"""
import json
from fractions import Fraction
from pathlib import Path

from audit_global_certificates import audit_mutations, verify as verify_graph
from linear_witness import verify as verify_linear
from nonrank_audit_template import known_pattern, verify as verify_template
from real_verify import verify as verify_puzzle

HERE = Path(__file__).resolve().parent


def main():
    path = HERE / "nonrank_full_old_stall.json"
    fixture = json.loads(path.read_text())
    templates = json.loads((HERE / "nonrank_template_primal.json").read_text())
    native = json.loads((HERE / "linear_real_nonrank_full_old_stall.json").read_text())
    replay = json.loads((HERE / "replay_nonrank_full_old_stall.json").read_text())
    puzzle_result = verify_puzzle(path, all_triples=True)
    template_result = verify_template(fixture, templates)
    native_result = verify_linear(fixture, native)
    graph_result = verify_graph(replay)
    assert replay["state"] == fixture["state"]
    assert replay["proof"]["cells"] == fixture["graph_cells"]
    assert {key: replay["proof"][key] for key in ("cell", "digit", "value")} == {
        "cell": 35, "digit": 9, "value": True}
    assert fixture["target"] == {"cell": 35, "digit": 9, "value": False}
    assert len(replay["proof"]["nodes"]) == 169
    assert sum(char in "123456789" for char in fixture["puzzle"]) == 28

    # The native witness is exactly the projection of the stronger template
    # witness, rather than an unrelated feasible point on different premises.
    projected = [[Fraction(0) for _ in range(9)] for _ in range(81)]
    for term in templates["terms"]:
        for cell in term["cells"]:
            projected[cell][term["digit"] - 1] += Fraction(term["numerator"], term["denominator"])
    dense = [[Fraction(0) for _ in range(9)] for _ in range(81)]
    for item in native["nonzero"]:
        dense[item["cell"]][item["digit"] - 1] = Fraction(item["numerator"], item["denominator"])
    assert projected == dense
    print(json.dumps({"puzzle": puzzle_result, "template": template_result,
                      "native_projection": native_result, "certificate": graph_result,
                      "certificate_mutations_rejected": audit_mutations(replay),
                      "classification": known_pattern(fixture)}, sort_keys=True))


if __name__ == "__main__":
    if not __debug__:
        raise SystemExit("Run without -O: exact acceptance uses assertions.")
    main()
