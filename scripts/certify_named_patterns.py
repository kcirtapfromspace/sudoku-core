#!/usr/bin/env python3
"""Offline witness discovery. Acceptance belongs to the independent replay checkers."""
import itertools as it
import json
from pathlib import Path
from verify_technique_examples import peers, from_givens

ROOT = Path(__file__).resolve().parent.parent
HOUSE = [
    [
        c
        for c in range(81)
        if (
            c // 9 == h
            if h < 9
            else c % 9 == h - 9 if h < 18 else c // 27 * 3 + c % 9 // 3 == h - 18
        )
    ]
    for h in range(27)
]


def digits(m):
    return [d for d in range(1, 10) if m >> (d - 1) & 1]


def ep(h):
    placement = "SetValue" in h["hint_type"]
    a = h["hint_type"]["SetValue" if placement else "EliminateCandidates"]
    return (
        a["pos"]["row"] * 9 + a["pos"]["col"],
        [a["value"]] if placement else a["values"],
        placement,
    )


class State:
    def __init__(self, s):
        self.s = s
        self.dom = s["domains"]
        self.val = s["values"]
        self.empty = [c for c in range(81) if not self.val[c]]

    def live(self, c, d):
        return c in self.empty and d in digits(self.dom[c])

    def loc(self, h, d):
        return [c for c in HOUSE[h] if self.live(c, d)]

    def union(self, cs):
        return set().union(*(set(digits(self.dom[c])) for c in cs))

    def fish(self, d, b, co):
        bs = [set(self.loc(h, d)) for h in b]
        cov = [set(self.loc(h, d)) for h in co]
        if any(not s for s in bs + cov) or any(
            a & z for hs in (bs, cov) for i, a in enumerate(hs) for z in hs[:i]
        ):
            return None
        base = set.union(*bs)
        cover = set.union(*cov)
        return {
            "digit": d,
            "bases": list(b),
            "covers": list(co),
            "fins": sorted(base - cover),
        }

    def targets(self, f):
        base = set().union(*(set(self.loc(h, f["digit"])) for h in f["bases"]))
        cover = set().union(*(set(self.loc(h, f["digit"])) for h in f["covers"]))
        return [c for c in sorted(cover - base) if all(peers(c, z) for z in f["fins"])]

    def basic_fish(self, d, n):
        for row in (True, False):
            bp = range(9) if row else range(9, 18)
            cp = range(9, 18) if row else range(9)
            for b in it.combinations([h for h in bp if len(self.loc(h, d)) >= 2], n):
                for co in it.combinations(cp, n):
                    f = self.fish(d, b, co)
                    if f and f["fins"] and self.targets(f):
                        yield f


def hint(t, c, ds, cells, placement=False):
    return {
        "technique": t,
        "hint_type": (
            {"SetValue": {"pos": {"row": c // 9, "col": c % 9}, "value": ds[0]}}
            if placement
            else {
                "EliminateCandidates": {
                    "pos": {"row": c // 9, "col": c % 9},
                    "values": ds,
                }
            }
        ),
        "explanation": f"{t}: independently specified pattern witness.",
        "involved_cells": [{"row": o // 9, "col": o % 9} for o in cells],
    }


def subset(s, h):
    c, ds, p = ep(h)
    t = h["technique"]
    n = 1 if "Single" in t else 2 if "Pair" in t else 3 if "Triple" in t else 4
    hidden = t.startswith("Hidden")
    for ho in range(27):
        if c not in HOUSE[ho]:
            continue
        if n == 1 and (
            (hidden and s.loc(ho, ds[0]) == [c])
            or (not hidden and digits(s.dom[c]) == ds)
        ):
            return {"Subset": {"house": ho, "cells": [c], "digits": ds}}
        if n == 1:
            continue
        if hidden:
            for dd in it.combinations(range(1, 10), n):
                ls = sorted(set().union(*(set(s.loc(ho, d)) for d in dd)))
                if (
                    len(ls) == n
                    and c in ls
                    and all(s.loc(ho, d) for d in dd)
                    and not set(dd) & set(ds)
                ):
                    return {"Subset": {"house": ho, "cells": ls, "digits": list(dd)}}
        else:
            for cs in it.combinations(
                [o for o in HOUSE[ho] if o in s.empty and o != c], n
            ):
                dd = s.union(cs)
                if len(dd) == n and set(ds) <= dd:
                    return {
                        "Subset": {"house": ho, "cells": list(cs), "digits": sorted(dd)}
                    }
    raise ValueError(t)


def wing(s, t):
    for p in s.empty:
        pd = digits(s.dom[p])
        if len(pd) != (2 if t == "XYWing" else 3):
            continue
        for a, b in it.combinations(
            [o for o in s.empty if peers(p, o) and s.dom[o].bit_count() == 2], 2
        ):
            aa = set(digits(s.dom[a]))
            bb = set(digits(s.dom[b]))
            common = aa & bb
            if len(common) != 1:
                continue
            z = common.pop()
            d = next(iter(aa - {z}))
            e = next(iter(bb - {z}))
            if d == e or set(pd) != {d, e} | ({z} if t == "XYZWing" else set()):
                continue
            for c in s.empty:
                if (
                    c not in (p, a, b)
                    and s.live(c, z)
                    and peers(c, a)
                    and peers(c, b)
                    and (t == "XYWing" or peers(c, p))
                ):
                    return hint(t, c, [z], [p, a, b, c]), {
                        "Wing": {"cells": [p, a, b], "digits": [d, e, z]}
                    }


def literal_conflict(a, b):
    return a != b and (a // 9 == b // 9 or a % 9 == b % 9 and peers(a // 9, b // 9))


def strong(s, a, b):
    if not literal_conflict(a, b):
        return False
    return (
        s.dom[a // 9].bit_count() == 2
        if a // 9 == b // 9
        else any(s.loc(h, a % 9 + 1) == sorted([a // 9, b // 9]) for h in range(27))
    )


def propagate(s, assumption, endpoint, digit_only=None):
    assigned = [0] * 729
    for c in range(81):
        for d in range(1, 10):
            n = c * 9 + d - 1
            if not s.dom[c] >> (d - 1) & 1:
                assigned[n] = -1
            elif s.val[c]:
                assigned[n] = 1
    assigned[assumption] = 1
    steps = []
    req = [[c * 9 + d - 1 for d in digits(s.dom[c])] for c in range(81)] + [
        [c * 9 + d - 1 for c in HOUSE[h] if s.dom[c] >> (d - 1) & 1]
        for h in range(27)
        for d in range(1, 10)
    ]

    def result(conflict):
        return {"assumption": assumption, "steps": steps, "conflict": conflict}

    for _ in range(730):
        changed = False
        ons = [n for n, v in enumerate(assigned) if v == 1]
        for n in range(729):
            if assigned[n] or digit_only and n % 9 + 1 != digit_only:
                continue
            cause = next((o for o in ons if literal_conflict(n, o)), None)
            if cause is not None:
                assigned[n] = -1
                steps.append({"candidate": n, "on": False, "reason": cause})
                changed = True
        for r, ns in enumerate(req):
            if digit_only and (r < 81 or (r - 81) % 9 + 1 != digit_only):
                continue
            if any(assigned[n] == 1 for n in ns):
                continue
            ls = [n for n in ns if assigned[n] != -1]
            if not ls:
                return result(r)
            if len(ls) == 1 and not assigned[ls[0]]:
                assigned[ls[0]] = 1
                steps.append({"candidate": ls[0], "on": True, "reason": r})
                changed = True
        if endpoint:
            c, ds, p = endpoint
            if all(assigned[c * 9 + d - 1] == (1 if p else -1) for d in ds):
                return result(None)
        if not changed:
            return None
    return None


def build(record):
    h = record["hint"]
    t = h["technique"]
    s = State(record["source"]["state"])
    p = record["pattern"]
    c, ds, place = ep(h)
    if p and "Basic" in p:
        return subset(s, h)
    if t in ("XYWing", "XYZWing"):
        found = wing(s, t)
        if found:
            record["hint"], rule = found
            record["pattern"] = None
            return rule
        return None
    if t == "WXYZWing":
        cs = sorted(o for a in p["Als"]["als_chain"] for o in a["cells"])
        z = p["Als"]["z_value"]
        dd = sorted(s.union(cs) - {z}) + [z]
        return {"Wing": {"cells": cs, "digits": dd}}
    if p and "Fish" in p:
        f = p["Fish"]
        part = {
            "digit": f["digit"],
            "bases": f["base_sectors"],
            "covers": f["cover_sectors"],
            "fins": f["fins"],
        }
        if t == "MutantFish":
            d = f["digit"]
            houses = [ho for ho in range(27) if s.loc(ho, d)]
            for n in (2, 3):
                for b in it.combinations(houses, n):
                    if not all(len(s.loc(h, d)) >= 2 for h in b) or not (
                        any(z < 9 for z in b) and any(9 <= z < 18 for z in b)
                    ):
                        continue
                    for co in it.combinations([z for z in houses if z not in b], n):
                        ff = s.fish(d, b, co)
                        if (
                            ff
                            and (
                                not ff["fins"]
                                or len({z // 27 * 3 + z % 9 // 3 for z in ff["fins"]})
                                == 1
                            )
                            and all(
                                set(s.loc(h, d))
                                & set().union(*(set(s.loc(z, d)) for z in b))
                                for h in co
                            )
                            and s.targets(ff)
                        ):
                            target = s.targets(ff)[0]
                            record["hint"] = hint(
                                t,
                                target,
                                [d],
                                sorted(
                                    set().union(*(set(s.loc(z, d)) for z in b + co))
                                ),
                            )
                            record["pattern"] = {
                                "Fish": {
                                    "digit": d,
                                    "base_sectors": list(b),
                                    "cover_sectors": list(co),
                                    "fins": ff["fins"],
                                }
                            }
                            return {"Fish": {"parts": [ff]}}
            return None
        if t == "SiameseFish":
            for d in (f["digit"],):
                for n in (2, 3, 4):
                    groups = {}
                    for ff in s.basic_fish(d, n):
                        key = tuple(ff["bases"])
                        for gg in groups.get(key, []):
                            if len(set(ff["covers"]) - set(gg["covers"])) == 1 and any(
                                o not in s.targets(ff) for o in s.targets(gg)
                            ):
                                target = s.targets(ff)[0]
                                record["hint"] = hint(
                                    t,
                                    target,
                                    [d],
                                    sorted(
                                        set().union(
                                            *(
                                                set(s.loc(z, d))
                                                for z in ff["bases"]
                                                + ff["covers"]
                                                + gg["covers"]
                                            )
                                        )
                                    ),
                                )
                                record["pattern"] = {
                                    "Fish": {
                                        "digit": d,
                                        "base_sectors": ff["bases"],
                                        "cover_sectors": ff["covers"],
                                        "fins": ff["fins"],
                                    }
                                }
                                return {"Fish": {"parts": [ff, gg]}}
                        groups.setdefault(key, []).append(ff)
            return None
        return {"Fish": {"parts": [part]}}
    if p and "Als" in p:
        a = p["Als"]
        sets = [
            {"cells": z["cells"], "digits": z["candidates"], "house": z["sector"]}
            for z in a["als_chain"]
        ]
        if t.startswith("Aligned"):
            cells = [c] + [z["cells"][0] for z in a["als_chain"][:-1]]
            return {
                "Aligned": {"cells": cells, "exclusions": [a["als_chain"][-1]["cells"]]}
            }
        return {"Als": {"sets": sets, "links": a["rcc_values"], "z": a["z_value"]}}
    if p and "Aic" in p:
        a = p["Aic"]
        ns = [z[0] * 9 + z[1] - 1 for z in a["chain"]]
        med = t == "ThreeDMedusa"
        return {
            "Chain": {
                "nodes": ns,
                "strong": [] if med else [z == "Strong" for z in a["link_types"]],
                "colors": [z[2] == "On" for z in a["chain"]] if med else [],
            }
        }
    if t == "BivalueUniversalGrave":
        return "Bug"
    if p and "Uniqueness" in p:
        u = p["Uniqueness"]
        cs = sorted(set(u["floor_cells"] + u["roof_cells"]))
        common = set.intersection(*(set(digits(s.dom[o])) for o in cs))
        dd = sorted(common)
        if t == "AvoidableRectangle":
            dd = sorted(set(s.val[o] for o in cs if s.val[o]))
        if t == "ExtendedUniqueRectangle":
            dd = sorted(s.union(u["floor_cells"]))
        return {"Trade": {"cells": cs, "digits": dd}}
    if p and any(z in p for z in ("Arithmetic", "ThreeHouse", "GlobalColoring")):
        return {"Native": {"proof": p}}
    if t == "Backtracking":
        return "Search"
    if p and "Forcing" in p:
        src = p["Forcing"]["source"]
        fish = None
        if t == "KrakenFish":
            for r1, r2 in it.combinations(range(9), 2):
                common = set(o % 9 for o in s.loc(r1, ds[0])) & set(
                    o % 9 for o in s.loc(r2, ds[0])
                )
                if len(common) != 2:
                    continue
                f = s.fish(ds[0], [r1, r2], [z + 9 for z in sorted(common)])
                if not f or not f["fins"] or c // 9 in (r1, r2) or c % 9 not in common:
                    continue
                branches = [propagate(s, z * 9 + ds[0] - 1, ep(h)) for z in f["fins"]]
                if all(branches):
                    return {"Forcing": {"source": src, "fish": f, "branches": branches}}
            return None
        if "Nishio" in src:
            assumptions = [src["Nishio"]["cell"] * 9 + src["Nishio"]["digit"] - 1]
        elif "Cell" in src:
            assumptions = [src["Cell"] * 9 + d - 1 for d in digits(s.dom[src["Cell"]])]
        else:
            a = src["Region"]
            assumptions = [
                o * 9 + a["digit"] - 1 for o in s.loc(a["sector"], a["digit"])
            ]
        branches = [
            propagate(s, n, ep(h), ds[0] if t == "NishioForcingChain" else None)
            for n in assumptions
        ]
        if all(branches):
            return {"Forcing": {"source": src, "fish": fish, "branches": branches}}
        return None
    raise ValueError(t)


def describe(record, rule):
    """Explain the checked local rule in grid terms."""
    h = record["hint"]
    t = h["technique"]
    c, removed, _ = ep(h)
    cell = lambda o: f"r{o//9+1}c{o%9+1}"
    houses = lambda hs: ", ".join(
        (f"row {h+1}" if h < 9 else f"column {h-8}" if h < 18 else f"box {h-17}")
        for h in hs
    )
    if t in ("XYWing", "XYZWing"):
        w = rule["Wing"]
        p, a, b = w["cells"]
        x, y, z = w["digits"]
        text = f"The pivot {cell(p)} sees pincers {cell(a)} and {cell(b)}. If the pivot is {x}, {cell(a)} must be {z}. If the pivot is {y}, {cell(b)} must be {z}. "
        if t == "XYZWing":
            text += f"The pivot can also be {z}. The target sees all three possible {z} locations. "
        else:
            text += f"The target sees both possible {z} locations. "
        text += f"Remove {z} from {cell(c)}."
    elif t == "SiameseFish":
        a, b = rule["Fish"]["parts"]
        text = f"Two finned fish on {a['digit']} share bases {houses(a['bases'])}. Their covers are {houses(a['covers'])} and {houses(b['covers'])}. Each fish has a valid elimination that sees all its fins. The first fish removes {a['digit']} from {cell(c)}."
    elif t == "MutantFish":
        f = rule["Fish"]["parts"][0]
        text = f"The bases {houses(f['bases'])} mix a row with a column. Without a fin, digit {f['digit']} fills the covers {houses(f['covers'])}. The target sees every fin, so remove {f['digit']} from {cell(c)}."
    elif t == "ExtendedUniqueRectangle":
        text = f"Six cells form a three-box rectangle on digits 1, 3 and 5. If {cell(c)} used 1 or 5, both selected columns would contain these three digits. Swapping the columns within the rectangle would give a second solution. Remove 1 and 5 from {cell(c)}."
    elif t == "NishioForcingChain":
        b = rule["Forcing"]["branches"][0]
        r = b["conflict"]
        text = f"Assume {cell(c)} is {removed[0]}. Follow only this digit's strong and weak implications. The assumption leaves {houses([(r-81)//9])} without a position for {removed[0]}. Remove {removed[0]} from {cell(c)}."
    else:
        return
    h["explanation"] = text


def main():
    corpus = json.loads((ROOT / "tests/fixtures/technique_examples.json").read_text())
    out = ROOT / "target/named-patterns"
    out.mkdir(exist_ok=True)
    state_path = ROOT / "target/hodoku-states.json"
    allstates = json.loads(state_path.read_text()) if state_path.exists() else []
    for e in corpus["examples"]:
        t = e["hint"]["technique"]
        puz = corpus["puzzles"][e["puzzle"]]
        source = {**puz, "puzzle": puz["givens"], "state": e["state"]}
        record = {
            "source": source,
            "hint": e["hint"],
            "pattern": e["certificate"]["Certified"]["pattern"],
        }
        rule = build(record)
        if rule is None and t in ("XYWing", "XYZWing"):
            for r in allstates:
                found = wing(State(r["state"]), t)
                if found:
                    record = {"source": r, "hint": found[0], "pattern": None}
                    rule = found[1]
                    break
        if rule is None:
            raise ValueError("No named witness for " + t)
        describe(record, rule)
        record["named_pattern"] = rule
        (out / f"{t}.json").write_text(json.dumps(record, indent=2) + "\n")
        print(t, "OK" if rule else "MISSING", flush=True)


if __name__ == "__main__":
    main()
