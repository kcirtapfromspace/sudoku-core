"""Exact acceptance of the unique four-digit nonforcing-gap control.

This is deliberately NOT a full-old-engine stall: its saved full dispatcher
finds Nishio r4c6 != 4. Recompiling the frozen probe reproduces that comparison;
this standard-library check validates the mathematical state and certificate.
"""
import json
import hashlib
from pathlib import Path

from audit_global_certificates import verify as verify_certificate
from real_verify import replay_trace, verify as verify_game

HERE = Path(__file__).resolve().parent


def main():
    if not __debug__:
        raise SystemExit("Run without -O: assertions are part of exact acceptance")
    path = HERE / "expanded_unique_nonforcing.json"
    case = json.loads(path.read_text())
    result = verify_game(path)
    assert result["solutions_capped_at_2"] == 1
    assert sum(char in "123456789" for char in case["puzzle"]) == 29
    assert case["state"]["values"][16] == 0
    assert case["state"]["domains"][16] == sum(1 << (d - 1) for d in (1, 4, 6, 7, 8))
    for field in ("old_fixedpoint", "old_full_comparison"):
        comparison = case[field]
        assert comparison["state"] == case["state"] and comparison["trace"] == []
        assert comparison["uniqueness_warranted"] is True
        assert comparison["arithmetic"]["hint"] is None
        assert comparison["arithmetic"]["error"] is None
        assert comparison["three_house"]["hint"] is None
        assert comparison["three_house"]["centers_checked"] == 81
    assert case["old_full_comparison"]["full_current_finding"] == ["NishioForcingChain", 32, [4], False]
    assert int(case["solution"][32]) != 4
    replay = json.loads((HERE / "replay_expanded_unique_nonforcing.json").read_text())
    assert replay["state"] == case["state"]
    proof = replay["proof"]
    assert (proof["cell"], proof["digit"], proof["value"]) == (16, 1, True)
    result["production_certificate"] = verify_certificate(replay)
    old = json.loads((HERE / "expanded_cover_old_forcing_trace_29.json").read_text())
    assert old["source_fixture_sha256"] == hashlib.sha256(path.read_bytes()).hexdigest()
    assert old["puzzle"] == case["puzzle"] and old["solution"] == case["solution"]
    assert old["start_state"] == case["state"]
    assert old["trace_from_givens"][:len(case["base_trace"])] == case["base_trace"]
    end = replay_trace(case["initial_domains"], old["trace_from_givens"],
                       list(map(int, case["solution"])))
    assert end == old["end_domains"] == old["end_state"]["domains"]
    assert old["result"] == "guardian-placed-by-older-methods" and end[16] == 1
    assert old["technique_counts_after_start"] == {"NishioForcingChain": 5}
    result["old_forcing_trace"] = "five Nishio steps plus replayed local closure place the guardian"
    result["comparison_tier"] = "nonforcing plus default arithmetic gap; full old Nishio finds a step"
    print(json.dumps(result, sort_keys=True))


if __name__ == "__main__":
    main()
