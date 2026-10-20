"""Independent declarative pattern checker. No detector, Rust call or global solver."""

import hashlib
import itertools as it
import math
from verify_technique_examples import require, peers, from_givens, check as exact_check

HOUSES = [
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


def ds(mask):
    return [d for d in range(1, 10) if mask >> (d - 1) & 1]


def bit(d):
    require(type(d) is int and 1 <= d <= 9, "Invalid pattern digit")
    return 1 << (d - 1)


def distinct(xs, limit):
    return (
        bool(xs)
        and all(type(c) is int and 0 <= c < limit for c in xs)
        and len(set(xs)) == len(xs)
    )


def digit_mask(xs):
    require(distinct([d - 1 for d in xs], 9), "Invalid pattern digit set")
    return sum(bit(d) for d in xs)


def digest(s, prefix=b"sudoku-exact-cover-v1-native"):
    h = hashlib.sha256(prefix)
    for v, d in zip(s["values"], s["domains"]):
        h.update(bytes([v]))
        h.update(d.to_bytes(2, "little"))
    return h.hexdigest()


def conflict(a, b):
    return a != b and (a // 9 == b // 9 or a % 9 == b % 9 and peers(a // 9, b // 9))


def verify_named(proof, state, hint, givens, deduction):
    method = hint["technique"]
    action = hint["hint_type"]
    placed = "SetValue" in action
    a = action["SetValue" if placed else "EliminateCandidates"]
    target = a["pos"]["row"] * 9 + a["pos"]["col"]
    digits = [a["value"]] if placed else a["values"]
    dom = state["domains"]
    values = state["values"]
    require(
        set(proof)
        == {
            "version",
            "state_hash",
            "technique",
            "conclusion",
            "witness",
            "uniqueness",
        },
        "Unknown named proof field",
    )
    require(
        proof["version"] == 1
        and proof["technique"] == method
        and proof["state_hash"] == digest(state)
        and proof["conclusion"] == deduction["purpose"],
        "Wrong named proof binding",
    )
    require(
        proof["conclusion"]
        == {"Deduction": {"cell": target, "digits": digits, "placement": placed}},
        "Wrong named endpoint",
    )
    require(
        (proof["uniqueness"] is not None)
        == (
            method
            in (
                "UniqueRectangle",
                "HiddenRectangle",
                "AvoidableRectangle",
                "ExtendedUniqueRectangle",
                "BivalueUniversalGrave",
            )
        ),
        "Wrong uniqueness premise",
    )
    require(
        all(
            not values[c] or not values[o] or not peers(c, o) or values[c] != values[o]
            for c in range(81)
            for o in range(c)
        ),
        "Conflicting placed values",
    )
    require(
        all(
            values[c]
            or all(
                not v or not peers(c, o) or not dom[c] & bit(v)
                for o, v in enumerate(values)
            )
            for c in range(81)
        ),
        "Pattern candidate conflicts with a placement",
    )
    digit_mask(digits)
    require(
        0 <= target < 81
        and not values[target]
        and all(dom[target] & bit(d) for d in digits),
        "Invalid named endpoint",
    )

    def live(c, d):
        return 0 <= c < 81 and not values[c] and bool(dom[c] & bit(d))

    def loc(h, d):
        require(0 <= h < 27, "Invalid house")
        return [c for c in HOUSES[h] if live(c, d)]

    def union(cs):
        return set().union(*(set(ds(dom[c])) for c in cs))

    def sees(cs, d):
        ls = [c for c in cs if live(c, d)]
        return bool(ls) and all(peers(target, c) for c in ls)

    def elim(d):
        return not placed and digits == [d]

    def link(a, b, strong):
        if (
            not distinct([a, b], 729)
            or not live(a // 9, a % 9 + 1)
            or not live(b // 9, b % 9 + 1)
            or not conflict(a, b)
        ):
            return False
        return not strong or (
            dom[a // 9].bit_count() == 2
            if a // 9 == b // 9
            else any(loc(h, a % 9 + 1) == sorted([a // 9, b // 9]) for h in range(27))
        )

    def fish(f, endpoint=True):
        require(set(f) == {"digit", "bases", "covers", "fins"}, "Invalid fish schema")
        d, b, co, fin = f["digit"], f["bases"], f["covers"], f["fins"]
        bit(d)
        if (
            not distinct(b, 27)
            or not distinct(co, 27)
            or len(b) != len(co)
            or len(b) > 4
            or set(b) & set(co)
        ):
            return False
        groups = [[set(loc(h, d)) for h in hs] for hs in (b, co)]
        if any(not s for g in groups for s in g) or any(
            a & b for g in groups for i, a in enumerate(g) for b in g[:i]
        ):
            return False
        bases = set.union(*groups[0])
        covers = set.union(*groups[1])
        fins = bases - covers
        return (
            len(set(fin)) == len(fin)
            and set(fin) == fins
            and (
                not endpoint
                or elim(d)
                and target in covers - bases
                and all(peers(target, c) for c in fins)
            )
        )

    def als(a):
        return (
            set(a) == {"cells", "digits", "house"}
            and distinct(a["cells"], 81)
            and len(a["cells"]) <= 8
            and 0 <= a["house"] < 27
            and all(not values[c] and c in HOUSES[a["house"]] for c in a["cells"])
            and len(a["digits"]) == len(a["cells"]) + 1
            and set(a["digits"]) == union(a["cells"])
            and len(set(a["digits"])) == len(a["digits"])
        )

    def rcc(a, b, d):
        return (
            d in a["digits"]
            and d in b["digits"]
            and all(
                peers(c, o)
                for c in a["cells"]
                if live(c, d)
                for o in b["cells"]
                if live(o, d)
            )
        )

    def assignments(cs):
        def visit(vs):
            if len(vs) == len(cs):
                yield vs
                return
            i = len(vs)
            for d in ds(dom[cs[i]]):
                if all(d != v or not peers(cs[i], cs[j]) for j, v in enumerate(vs)):
                    yield from visit(vs + [d])

        yield from visit([])

    def opposite(cs, vs):
        return target not in cs or (vs[cs.index(target)] in digits) != placed

    def unique():
        u = proof["uniqueness"]
        require(
            u
            and set(u) == {"givens", "proof"}
            and u["givens"] == givens
            and u["proof"]["purpose"] == "Unique",
            "Missing original uniqueness",
        )
        exact_check(u["proof"], from_givens(givens))
        sol = u["proof"]["witness"]
        start = from_givens(givens)
        require(
            all(
                (not start["values"][c] or values[c] == start["values"][c])
                and values[c] in (0, sol[c])
                and dom[c] & bit(sol[c])
                for c in range(81)
            ),
            "Unique premise mismatch",
        )

    def requirement(r):
        require(type(r) is int and 0 <= r < 324, "Invalid forcing requirement")
        return (
            [r * 9 + d - 1 for d in ds(dom[r])]
            if r < 81
            else [
                c * 9 + (r - 81) % 9
                for c in HOUSES[(r - 81) // 9]
                if dom[c] >> ((r - 81) % 9) & 1
            ]
        )

    def branch(b, digit_only=None):
        require(set(b) == {"assumption", "steps", "conflict"}, "Invalid branch schema")
        n = b["assumption"]
        require(
            type(n) is int
            and 0 <= n < 729
            and live(n // 9, n % 9 + 1)
            and len(b["steps"]) <= 729,
            "Invalid assumption",
        )
        status = [
            -1 if not dom[c] >> (d - 1) & 1 else 1 if values[c] else 0
            for c in range(81)
            for d in range(1, 10)
        ]
        status[n] = 1
        for step in b["steps"]:
            require(
                set(step) == {"candidate", "on", "reason"},
                "Invalid forcing step schema",
            )
            n, on, r = step["candidate"], step["on"], step["reason"]
            require(
                type(n) is int
                and 0 <= n < 729
                and status[n] == 0
                and type(on) is bool
                and (not digit_only or n % 9 + 1 == digit_only),
                "Invalid forcing literal",
            )
            if on:
                ns = requirement(r)
                require(
                    (not digit_only or r >= 81)
                    and n in ns
                    and all(status[o] == -1 for o in ns if o != n),
                    "False strong inference",
                )
            else:
                require(
                    type(r) is int
                    and 0 <= r < 729
                    and status[r] == 1
                    and conflict(n, r),
                    "False weak inference",
                )
            status[n] = 1 if on else -1
        if b["conflict"] is not None:
            r = b["conflict"]
            require(
                (not digit_only or r >= 81)
                and all(status[n] == -1 for n in requirement(r)),
                "False forcing conflict",
            )
        else:
            require(
                all(
                    status[target * 9 + d - 1] == (1 if placed else -1) for d in digits
                ),
                "Branch lacks conclusion",
            )

    w = proof["witness"]
    kind = w if isinstance(w, str) else next(iter(w))
    p = {} if isinstance(w, str) else w[kind]
    require(isinstance(w, str) or len(w) == 1, "Invalid witness variant")
    ok = False
    if kind == "Subset":
        require(set(p) == {"house", "cells", "digits"}, "Invalid subset schema")
        h, cs, dd = p["house"], p["cells"], p["digits"]
        m = digit_mask(dd)
        require(
            0 <= h < 27
            and distinct(cs, 81)
            and all(not values[c] and c in HOUSES[h] for c in cs),
            "Invalid subset cells",
        )
        allowed = {
            "NakedSingle": (1, False),
            "HiddenSingle": (1, True),
            "NakedPair": (2, False),
            "HiddenPair": (2, True),
            "NakedTriple": (3, False),
            "HiddenTriple": (3, True),
            "NakedQuad": (4, False),
            "HiddenQuad": (4, True),
        }
        require(method in allowed, "Subset label mismatch")
        n, hidden = allowed[method]
        ok = len(cs) == len(dd) == n
        if n == 1:
            ok = (
                ok
                and placed
                and cs == [target]
                and dd == digits
                and (loc(h, dd[0]) == cs if hidden else dom[target] == m)
            )
        elif hidden:
            ok = (
                ok
                and not placed
                and target in cs
                and not set(digits) & set(dd)
                and all(loc(h, d) and set(loc(h, d)) <= set(cs) for d in dd)
            )
        else:
            ok = (
                ok
                and not placed
                and target not in cs
                and target in HOUSES[h]
                and union(cs) == set(dd)
                and set(digits) <= set(dd)
            )
    elif kind == "Fish":
        require(set(p) == {"parts"} and 1 <= len(p["parts"]) <= 2, "Invalid fish parts")
        parts = p["parts"]
        f = parts[0]
        b, co, fin = f["bases"], f["covers"], f["fins"]
        n = len(b)
        basic = (all(h < 9 for h in b) and all(9 <= h < 18 for h in co)) or (
            all(9 <= h < 18 for h in b) and all(h < 9 for h in co)
        )
        mixed = lambda hs: any(h < 9 for h in hs) and any(9 <= h < 18 for h in hs)
        ok = fish(f)
        if method == "SiameseFish":
            require(len(parts) == 2, "Missing Siamese constituent")
            g = parts[1]
            base = set().union(*(set(loc(h, g["digit"])) for h in g["bases"]))
            cover = set().union(*(set(loc(h, g["digit"])) for h in g["covers"]))
            ok = (
                ok
                and fish(g, False)
                and f["digit"] == g["digit"]
                and b == g["bases"]
                and len(co) == len(g["covers"])
                and len(set(co) - set(g["covers"])) == 1
                and bool(fin)
                and bool(g["fins"])
                and any(
                    c != target and all(peers(c, z) for z in g["fins"])
                    for c in cover - base
                )
            )
        else:
            require(len(parts) == 1, "Unexpected constituent")
            if method == "PointingPair":
                ok = (
                    ok
                    and n == 1
                    and b[0] >= 18
                    and co[0] < 18
                    and not fin
                    and 2 <= len(loc(b[0], f["digit"])) <= 3
                )
            elif method == "BoxLineReduction":
                ok = ok and n == 1 and b[0] < 18 and co[0] >= 18 and not fin
            elif method in (
                "XWing",
                "Swordfish",
                "Jellyfish",
                "FinnedXWing",
                "FinnedSwordfish",
                "FinnedJellyfish",
            ):
                name = method.removeprefix("Finned")
                ok = (
                    ok
                    and basic
                    and n == {"XWing": 2, "Swordfish": 3, "Jellyfish": 4}[name]
                    and bool(fin) == method.startswith("Finned")
                    and (not fin or len({c // 27 * 3 + c % 9 // 3 for c in fin}) == 1)
                )
            elif method == "FrankenFish":
                ok = (
                    ok
                    and n >= 2
                    and any(h >= 18 for h in b + co)
                    and not mixed(b)
                    and not mixed(co)
                )
            elif method == "MutantFish":
                ok = ok and n >= 2 and (mixed(b) or mixed(co))
            elif method == "EmptyRectangle":
                boxes = [h for h in b if h >= 18]
                ok = (
                    ok
                    and n == 2
                    and len(boxes) == 1
                    and bool(fin)
                    and all(h < 18 for h in co)
                )
                if ok:
                    ls = loc(boxes[0], f["digit"])
                    ok = any(
                        all(o // 9 == r or o % 9 == c for o in ls)
                        and any(o // 9 == r and o % 9 != c for o in ls)
                        and any(o % 9 == c and o // 9 != r for o in ls)
                        for r in range(9)
                        for c in range(9)
                    )
            else:
                ok = False
    elif kind == "Wing":
        require(set(p) == {"cells", "digits"}, "Invalid wing schema")
        cs, dd = p["cells"], p["digits"]
        digit_mask(dd)
        require(
            distinct(cs, 81) and all(not values[c] for c in cs), "Invalid wing cells"
        )
        z = dd[-1]
        ok = elim(z) and target not in cs
        if method == "WXYZWing":
            ok = (
                ok
                and len(cs) == len(dd) == 4
                and union(cs) == set(dd)
                and sees(cs, z)
                and all(
                    all(
                        a == b or peers(a, b)
                        for a in cs
                        if live(a, d)
                        for b in cs
                        if live(b, d)
                    )
                    for d in dd
                    if d != z
                )
                and any(
                    a != b and not peers(a, b)
                    for a in cs
                    if live(a, z)
                    for b in cs
                    if live(b, z)
                )
            )
        else:
            require(
                method in ("XYWing", "XYZWing") and len(cs) == len(dd) == 3,
                "Wing label mismatch",
            )
            pivot, a, b = cs
            d, e, z = dd
            ok = (
                ok
                and set(ds(dom[a])) == {d, z}
                and set(ds(dom[b])) == {e, z}
                and set(ds(dom[pivot]))
                == {d, e} | ({z} if method == "XYZWing" else set())
                and peers(pivot, a)
                and peers(pivot, b)
                and peers(target, a)
                and peers(target, b)
                and (method == "XYWing" or peers(target, pivot))
            )
    elif kind == "Als":
        require(set(p) == {"sets", "links", "z"}, "Invalid ALS schema")
        sets, links, z = p["sets"], p["links"], p["z"]
        allcells = [c for a in sets for c in a["cells"]]
        require(1 <= len(sets) <= 8 and distinct(allcells, 81), "Overlapping ALSs")
        if method == "SueDeCoq":
            require(len(sets) == 3 and z is None and len(links) == 4, "SDC shape")
            i, b, l = sets
            ok = (
                len(i["cells"]) == 2
                and len(i["digits"]) == 4
                and 0 <= i["house"] < 18
                and 18 <= b["house"] < 27
                and l["house"] == i["house"]
                and als(b)
                and als(l)
                and set(i["digits"]) == union(i["cells"])
                and all(
                    not values[c]
                    and c in HOUSES[i["house"]]
                    and c in HOUSES[b["house"]]
                    for c in i["cells"]
                )
                and all(c not in HOUSES[i["house"]] for c in b["cells"])
                and all(c not in HOUSES[b["house"]] for c in l["cells"])
                and not set(b["digits"]) & set(l["digits"])
                and set(links[:2]) == set(b["digits"])
                and set(links[2:]) == set(l["digits"])
                and union(allcells) == set(i["digits"])
                and not placed
                and len(digits) == 1
                and target not in allcells
                and not any(
                    all(
                        not peers(target, c) or v != digits[0]
                        for c, v in zip(allcells, vs)
                    )
                    for vs in assignments(allcells)
                )
            )
        elif method == "DeathBlossom":
            require(
                len(sets) >= 3
                and len(sets[0]["cells"]) == 1
                and len(links) + 1 == len(sets),
                "Death Blossom shape",
            )
            stem = sets[0]["cells"][0]
            ok = (
                digit_mask(links) == dom[stem]
                and digit_mask(sets[0]["digits"]) == dom[stem]
                and 0 <= sets[0]["house"] < 27
                and stem in HOUSES[sets[0]["house"]]
                and z is not None
                and elim(z)
                and target not in allcells
                and all(
                    als(a)
                    and d != z
                    and d in a["digits"]
                    and z in a["digits"]
                    and all(peers(c, stem) for c in a["cells"] if live(c, d))
                    and sees(a["cells"], z)
                    for a, d in zip(sets[1:], links)
                )
            )
        else:
            count = {
                "AlsXz": len(sets) == 2,
                "AlsXyWing": len(sets) == 3,
                "AlsChain": len(sets) >= 4,
            }
            ok = (
                count.get(method, False)
                and len(links) + 1 == len(sets)
                and z is not None
                and elim(z)
                and target not in allcells
                and all(als(a) for a in sets)
                and all(
                    d != z
                    and (i == 0 or links[i - 1] != d)
                    and rcc(sets[i], sets[i + 1], d)
                    for i, d in enumerate(links)
                )
                and sees(sets[0]["cells"], z)
                and sees(sets[-1]["cells"], z)
            )
    elif kind == "Chain":
        require(set(p) == {"nodes", "strong", "colors"}, "Invalid chain schema")
        ns, ss, colors = p["nodes"], p["strong"], p["colors"]
        require(
            2 <= len(ns) <= 729
            and distinct(ns, 729)
            and all(live(n // 9, n % 9 + 1) for n in ns),
            "Invalid chain nodes",
        )
        if method == "ThreeDMedusa":
            require(
                len({n % 9 for n in ns}) >= 2
                and not ss
                and len(colors) == len(ns)
                and all(type(c) is bool for c in colors),
                "Invalid colors",
            )
            seen = {0}
            edges = [
                (i, j)
                for i in range(len(ns))
                for j in range(i)
                if link(ns[i], ns[j], True)
            ]
            require(
                all(colors[i] != colors[j] for i, j in edges),
                "Strong edge color conflict",
            )
            for _ in ns:
                for i, j in edges:
                    if i in seen or j in seen:
                        seen |= {i, j}
            valid = [
                col
                for col in (False, True)
                if not any(
                    colors[i] == colors[j] == col and conflict(ns[i], ns[j])
                    for i in range(len(ns))
                    for j in range(i)
                )
            ]
            ok = (
                len(seen) == len(ns)
                and bool(valid)
                and all(
                    any(
                        colors[i] == col
                        and (
                            n == target * 9 + d - 1
                            if placed
                            else conflict(n, target * 9 + d - 1)
                        )
                        for i, n in enumerate(ns)
                    )
                    for col in valid
                    for d in digits
                )
            )
        else:
            ok = (
                method in ("XChain", "AIC", "WWing")
                and not colors
                and len(ss) + 1 == len(ns)
                and len(ss) % 2 == 1
                and not placed
                and all(
                    type(v) is bool and v == (i % 2 == 0) and link(ns[i], ns[i + 1], v)
                    for i, v in enumerate(ss)
                )
            )
            if method == "XChain":
                ok = ok and len({n % 9 for n in ns}) == 1
            if method == "WWing":
                ok = (
                    ok
                    and len(ns) == 6
                    and ns[0] // 9 == ns[1] // 9
                    and ns[4] // 9 == ns[5] // 9
                    and ns[0] % 9 == ns[5] % 9
                    and len({n % 9 for n in ns[1:5]}) == 1
                    and dom[ns[0] // 9] == dom[ns[5] // 9]
                )
            ok = ok and all(
                conflict(target * 9 + d - 1, ns[0])
                and conflict(target * 9 + d - 1, ns[-1])
                for d in digits
            )
    elif kind == "Trade":
        require(set(p) == {"cells", "digits"}, "Invalid trade schema")
        cs, dd = p["cells"], p["digits"]
        m = digit_mask(dd)
        require(
            method
            in (
                "UniqueRectangle",
                "HiddenRectangle",
                "AvoidableRectangle",
                "ExtendedUniqueRectangle",
            )
            and distinct(cs, 81)
            and target in cs
            and not placed,
            "Trade label mismatch",
        )
        unique()
        start = from_givens(givens)
        require(all(not start["values"][c] for c in cs), "Trade changes a given")
        rows = {c // 9 for c in cs}
        cols = {c % 9 for c in cs}
        boxes = {c // 27 * 3 + c % 9 // 3 for c in cs}
        if method == "ExtendedUniqueRectangle":
            require(
                len(cs) == 6
                and len(dd) == 3
                and len(boxes) == 3
                and (len(rows), len(cols)) in ((2, 3), (3, 2))
                and sum(bool(dom[c] & ~m) for c in cs) == 1,
                "Extended rectangle shape",
            )
        else:
            require(
                len(cs) == 4
                and len(dd) == 2
                and len(rows) == len(cols) == len(boxes) == 2,
                "Rectangle shape",
            )
        require(
            sum(bool(values[c]) for c in cs)
            == (3 if method == "AvoidableRectangle" else 0),
            "Rectangle placements",
        )
        require(
            method != "HiddenRectangle" or any(dom[c] == m for c in cs),
            "No hidden rectangle floor",
        )
        require(
            (
                all(not values[c] or values[c] in dd for c in cs)
                if method == "AvoidableRectangle"
                else True
            ),
            "Avoidable core mismatch",
        )
        require(
            method in ("AvoidableRectangle", "ExtendedUniqueRectangle")
            or (all(dom[c] & m == m for c in cs) and any(dom[c] == m for c in cs)),
            "Rectangle core mismatch",
        )
        any_case = False
        ok = True
        for vs in assignments(cs):
            if not opposite(cs, vs):
                continue
            if not all(
                any(values[c] == d for c in HOUSES[h])
                or not set(loc(h, d)) <= set(cs)
                or any(c in loc(h, d) and v == d for c, v in zip(cs, vs))
                for h in range(27)
                for d in range(1, 10)
            ):
                continue
            any_case = True
            swaps = []
            for row in (True, False):
                if len(rows if row else cols) != 2:
                    continue
                other = [
                    vs[
                        next(
                            i
                            for i, o in enumerate(cs)
                            if (
                                o % 9 == c % 9 and o // 9 != c // 9
                                if row
                                else o // 9 == c // 9 and o % 9 != c % 9
                            )
                        )
                    ]
                    for c in cs
                ]
                swaps.append(
                    other != vs
                    and all(
                        sorted(v for c, v in zip(cs, vs) if c in hs)
                        == sorted(v for c, v in zip(cs, other) if c in hs)
                        for hs in HOUSES
                    )
                )
            ok = ok and any(swaps)
        ok = ok and any_case
    elif kind == "Bug":
        require(method == "BivalueUniversalGrave" and placed, "BUG label mismatch")
        unique()
        domains = dom.copy()
        domains[target] &= ~bit(digits[0])
        ok = all(values[c] or domains[c].bit_count() == 2 for c in range(81)) and all(
            sum(not values[c] and bool(domains[c] & bit(d)) for c in hs)
            == (0 if any(values[c] == d for c in hs) else 2)
            for hs in HOUSES
            for d in range(1, 10)
        )
    elif kind == "Aligned":
        require(set(p) == {"cells", "exclusions"}, "Invalid aligned schema")
        cs, ex = p["cells"], p["exclusions"]
        require(
            method in ("AlignedPairExclusion", "AlignedTripletExclusion")
            and len(cs) == (2 if method == "AlignedPairExclusion" else 3)
            and distinct(cs, 81)
            and target in cs
            and not placed
            and all(
                not values[c] and all(c == o or peers(c, o) for o in cs) for c in cs
            ),
            "Aligned base shape",
        )
        require(
            1 <= len(ex) <= 32
            and all(
                distinct(es, 81)
                and len(es) <= 4
                and all(
                    not values[c]
                    and c not in cs
                    and all(c == o or peers(c, o) for o in es)
                    for c in es
                )
                for es in ex
            ),
            "Invalid aligned excluder",
        )
        ok = True
        for vs in assignments(cs):
            if not opposite(cs, vs):
                continue
            excluded = False
            for es in ex:
                remaining = set()
                for c in es:
                    remaining |= set(ds(dom[c])) - {
                        v for o, v in zip(cs, vs) if peers(c, o)
                    }
                if len(remaining) < len(es):
                    excluded = True
            ok = ok and excluded
    elif kind == "Forcing":
        require(set(p) == {"source", "fish", "branches"}, "Invalid forcing schema")
        src, f, bs = p["source"], p["fish"], p["branches"]
        require(1 <= len(bs) <= 9, "Invalid branch count")
        expected = None
        if method in ("NishioForcingChain", "DynamicForcingChain") and "Nishio" in src:
            a = src["Nishio"]
            require(
                a["cell"] == target and elim(a["digit"]), "Nishio endpoint mismatch"
            )
            expected = [target * 9 + a["digit"] - 1]
        elif method in ("CellForcingChain", "DynamicForcingChain") and "Cell" in src:
            c = src["Cell"]
            require(
                0 <= c < 81 and not values[c] and dom[c].bit_count() >= 2,
                "Invalid forcing cell",
            )
            expected = requirement(c)
        elif method == "RegionForcingChain" and "Region" in src:
            a = src["Region"]
            bit(a["digit"])
            expected = requirement(81 + a["sector"] * 9 + a["digit"] - 1)
        elif method == "KrakenFish" and "Region" in src:
            require(
                f
                and fish(f, False)
                and len(f["bases"]) >= 2
                and src["Region"] == {"sector": f["bases"][0], "digit": f["digit"]}
                and f["fins"]
                and elim(f["digit"])
                and not any(target in HOUSES[h] for h in f["bases"])
                and any(target in HOUSES[h] for h in f["covers"]),
                "Invalid Kraken fish",
            )
            expected = [c * 9 + f["digit"] - 1 for c in f["fins"]]
        require(
            expected is not None
            and (method == "KrakenFish" or f is None)
            and sorted(expected) == sorted(b["assumption"] for b in bs)
            and len({b["assumption"] for b in bs}) == len(bs),
            "Incomplete forcing source cover",
        )
        for b in bs:
            branch(b, digits[0] if method == "NishioForcingChain" else None)
        ok = True
    elif kind == "Native":
        require(set(p) == {"proof"}, "Invalid native witness")
        native = p["proof"]
        require(len(native) == 1, "Invalid native proof")
        family = next(iter(native))
        require(
            {
                "Arithmetic": "ArithmeticCounting",
                "ThreeHouse": "ThreeHouseMatching",
                "GlobalColoring": "GlobalColoring",
            }.get(family)
            == method,
            "Native label mismatch",
        )
        q = native[family]
        require(
            q["cell"] == target
            and digits == [q["digit"]]
            and placed == (False if family == "ThreeHouse" else q["value"]),
            "Native endpoint mismatch",
        )
        verify_native(family, q, state)
        ok = True
    elif kind == "Search":
        require(method == "Backtracking", "Search label mismatch")
        require(
            set(p) == {"proof"} and p["proof"] == deduction, "Search trace mismatch"
        )
        exact_check(p["proof"], state)
        ok = True
    require(ok, "Named " + method + " pattern rejected")


def verify_native(family, p, s):
    dom, val = s["domains"], s["values"]
    c, d = p["cell"], p["digit"]
    if family == "Arithmetic":
        h = hashlib.sha256(b"sudoku-arithmetic-v1-classic-324")
        for v, m in zip(val, dom):
            h.update(m.to_bytes(2, "little"))
            h.update(bytes([v]))
        require(p["state_hash"] == h.hexdigest(), "Arithmetic state mismatch")
        terminal = p.get("terminal", "IntervalGcd")
        limit = (
            2
            if p["version"] == 1 and terminal == "IntervalGcd"
            else (
                8
                if p["version"] == 2
                and isinstance(terminal, dict)
                and "Residue" in terminal
                else 0
            )
        )
        require(limit and 1 <= len(p["terms"]) <= 6, "Arithmetic format")
        coeff = [0] * 729
        rhs = 0
        used = set()
        for t in p["terms"]:
            r = t["requirement"]
            w = t["weight"]
            require(type(w) is int and 0 < abs(w) <= limit, "Arithmetic weight")
            rhs += w
            if "Cell" in r:
                cell = r["Cell"]["cell"]
                require(0 <= cell < 81, "Arithmetic cell")
                key = cell
                ns = [cell * 9 + z - 1 for z in ds(dom[cell])]
            else:
                a = r["SectorDigit"]
                sector, digit = a["sector"], a["digit"]
                bit(digit)
                require(0 <= sector < 27, "Arithmetic house")
                key = 81 + sector * 9 + digit - 1
                ns = [o * 9 + digit - 1 for o in HOUSES[sector] if dom[o] & bit(digit)]
            require(key not in used, "Duplicate arithmetic term")
            used.add(key)
            for n in ns:
                coeff[n] += w
        n = c * 9 + d - 1
        residual = rhs - (0 if p["value"] else coeff[n])
        coeff[n] = 0
        lower = sum(min(0, k) for k in coeff)
        upper = sum(max(0, k) for k in coeff)
        if terminal == "IntervalGcd":
            gcd = math.gcd(*coeff)
            require(
                residual < lower
                or residual > upper
                or (residual != 0 if gcd == 0 else residual % gcd != 0),
                "False arithmetic contradiction",
            )
        else:
            modulus = terminal["Residue"]["modulus"]
            require(2 <= modulus <= 16, "Invalid modulus")
            reachable = {0}
            for k in coeff:
                if k:
                    reachable |= {(r + k) % modulus for r in reachable}
            require(residual % modulus not in reachable, "False arithmetic residue")
    elif family == "ThreeHouse":
        require(
            p["version"] == 1
            and p["state_hash"] == digest(s, b"sudoku-three-house-v1-native")
            and 0 <= p["center"] < 81,
            "Three-house binding",
        )
        center = p["center"]
        r, co = center // 9, center % 9
        hs = [r, 9 + co, 18 + r // 3 * 3 + co // 3]
        cells = sorted(set().union(*(set(HOUSES[h]) for h in hs)))
        atoms = [[] for _ in range(6)]
        for o in cells:
            rr, cc, bx = (
                o // 9 == r,
                o % 9 == co,
                o // 27 == center // 27 and o % 9 // 3 == co // 3,
            )
            a = (
                0
                if rr and cc
                else 1 if rr and bx else 2 if cc and bx else 3 if rr else 4 if cc else 5
            )
            atoms[a].append(o)
        witness = p["witness"]
        require(
            c in cells
            and len(witness) == 21
            and all(dom[o] & bit(v) for o, v in zip(cells, witness)),
            "Three-house witness",
        )
        vv = dict(zip(cells, witness))
        require(
            vv[c] != d
            and all({vv[o] for o in HOUSES[h]} == set(range(1, 10)) for h in hs),
            "Invalid three-house feasible witness",
        )
        domains = dom.copy()
        domains[c] = bit(d)
        blocked = set()
        require(len(p["blockers"]) <= 366, "Too many blockers")
        for b in p["blockers"]:
            a, m, sub = b["atom"], b["digits"], b["cells"]
            require(
                type(a) is int
                and 1 <= a <= 5
                and 0 <= m < 512
                and m.bit_count() == len(atoms[a])
                and 0 < sub < 1 << len(atoms[a])
                and (a, m) not in blocked,
                "Invalid Hall blocker",
            )
            union = 0
            for i, o in enumerate(atoms[a]):
                if sub >> i & 1:
                    union |= domains[o] & m
            require(union.bit_count() < sub.bit_count(), "False Hall obstruction")
            blocked.add((a, m))
        for z in ds(domains[center]):
            z = bit(z)
            for aa in it.combinations(
                [bit(i) for i in range(1, 10) if not bit(i) & z], 2
            ):
                a = sum(aa)
                for bb in it.combinations(
                    [bit(i) for i in range(1, 10) if not bit(i) & (z | a)], 2
                ):
                    b = sum(bb)
                    sets = [z, a, b, 511 ^ (z | a), 511 ^ (z | b), 511 ^ (z | a | b)]
                    require(
                        any((i, sets[i]) in blocked for i in range(1, 6)),
                        "Uncovered three-house signature",
                    )
    elif family == "GlobalColoring":
        require(
            p["version"] == 1
            and p["state_hash"] == digest(s, b"sudoku-global-coloring-v1-native"),
            "Coloring binding",
        )
        cs, wit, nodes = p["cells"], p["witness"], p["nodes"]
        require(
            distinct(cs, 81)
            and cs == sorted(cs)
            and 1 <= len(cs) <= 32
            and c in cs
            and len(wit) == len(cs)
            and 1 <= len(nodes) <= 32768,
            "Coloring shape",
        )
        require(
            all(
                dom[o] & bit(v)
                and all(not peers(o, cs[j]) or wit[j] != v for j in range(i))
                for i, (o, v) in enumerate(zip(cs, wit))
            )
            and (wit[cs.index(c)] == d) == p["value"],
            "Coloring witness",
        )
        domains = [dom[o] for o in cs]
        domains[cs.index(c)] &= ~bit(d) if p["value"] else bit(d)
        seen = set()

        def visit(n, domains, assigned):
            require(
                type(n) is int and 0 <= n < len(nodes) and n not in seen,
                "Coloring tree edge",
            )
            seen.add(n)
            node = nodes[n]
            if "Conflict" in node:
                require(
                    domains[cs.index(node["Conflict"]["cell"])] == 0,
                    "False coloring conflict",
                )
                return
            a = node["Split"]
            i = cs.index(a["cell"])
            require(
                i not in assigned
                and domains[i]
                and sorted(b["digit"] for b in a["branches"]) == ds(domains[i]),
                "Incomplete coloring split",
            )
            for b in a["branches"]:
                require(b["child"] < n, "Coloring postorder")
                next_domains = domains.copy()
                next_domains[i] = bit(b["digit"])
                for j, o in enumerate(cs):
                    if peers(o, cs[i]):
                        next_domains[j] &= ~bit(b["digit"])
                visit(b["child"], next_domains, assigned | {i})

        visit(len(nodes) - 1, domains, set())
        require(len(seen) == len(nodes), "Unused coloring node")
