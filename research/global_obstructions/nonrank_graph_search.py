"""Bounded search for native non-q-colorability with exact fractional q cover.

The rectangle mode selects q cells in each of four native boxes. Every
independent set in a fractional cover of total q must meet each box once.
Consequently at most q**4 templates suffice, and an integral q-coloring is
exactly a partition into q disjoint templates. Floating LP only proposes a
witness; every positive is checked with Fraction before being retained.
"""
import os
for _variable in ("OMP_NUM_THREADS","OPENBLAS_NUM_THREADS","MKL_NUM_THREADS"):
    os.environ[_variable] = "1"
import argparse
from collections import Counter
from fractions import Fraction
import itertools
import json
from pathlib import Path
import random
import time

import numpy as np
from scipy.optimize import linprog

HERE=Path(__file__).resolve().parent
ROWS=[set(range(9*r,9*r+9)) for r in range(9)]
COLS=[set(range(c,81,9)) for c in range(9)]
BOXES=[set((3*br+r)*9+3*bc+c for r in range(3) for c in range(3)) for br in range(3) for bc in range(3)]
HOUSES=ROWS+COLS+BOXES


def sees(a,b):
    return a//9 == b//9 or a%9 == b%9 or (a//27,a%9//3) == (b//27,b%9//3)


def templates(groups):
    return [tuple(cells) for cells in itertools.product(*groups)
            if all(not sees(a,b) for a,b in itertools.combinations(cells,2))]


def integral_partition(cells, independent_sets):
    cells=tuple(sorted(cells));index={c:i for i,c in enumerate(cells)}
    masks=[sum(1<<index[c] for c in independent) for independent in independent_sets]
    support=[[m for m in masks if m&(1<<i)] for i in range(len(cells))]
    calls=0;dead=set()
    def search(remaining):
        nonlocal calls
        calls+=1
        if not remaining:return []
        if remaining in dead:return None
        options=min(([m for m in support[i] if m&remaining==m] for i in range(len(cells)) if remaining&(1<<i)),key=len)
        for mask in options:
            result=search(remaining^mask)
            if result is not None:return [mask]+result
        dead.add(remaining);return None
    partition=search((1<<len(cells))-1)
    return partition,calls


def exact_fractional_cover(cells, independent_sets, q):
    if not independent_sets:return None
    matrix=np.array([[int(c in st) for st in independent_sets] for c in cells],dtype=np.int8)
    result=linprog(np.zeros(len(independent_sets)),A_eq=matrix,b_eq=np.ones(len(cells)),bounds=(0,None),method="highs-ds")
    if not result.success:
        if result.status!=2:raise RuntimeError(result.message)
        return None
    weights=[Fraction(float(w)).limit_denominator(1000000) for w in result.x]
    if not all(w>=0 for w in weights):raise RuntimeError("negative rationalization")
    if not all(sum(w for st,w in zip(independent_sets,weights) if c in st)==1 for c in cells):
        raise RuntimeError("rationalization failed exact coverage")
    if sum(weights)!=q:raise RuntimeError("fractional total differs from q")
    return [{"cells":list(st),"weight":[w.numerator,w.denominator]} for st,w in zip(independent_sets,weights) if w]


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--q",type=int,default=4,choices=[3,4,5,6])
    parser.add_argument("--samples",type=int,default=1000000)
    parser.add_argument("--seconds",type=int,default=600)
    parser.add_argument("--seed",type=int,default=219404)
    parser.add_argument("--output",type=Path,default=HERE/"nonrank_graph_rectangle_results.json")
    parser.add_argument("--max-positives",type=int,default=3)
    parser.add_argument("--tridagon-source",type=Path,help="Enumerate q-3 extra cells per box around this known12-cell control")
    args=parser.parse_args();args.output.parent.mkdir(parents=True,exist_ok=True)
    q=args.q;rng=random.Random(args.seed);stats=Counter();started=time.monotonic()
    # All rectangles are equivalent under native Sudoku row/column/band/stack
    # symmetries, so one representative suffices for these abstract supports.
    selected_boxes=[BOXES[i] for i in [0,1,3,4]]
    patterns=[list(itertools.combinations(sorted(box),q)) for box in selected_boxes]
    enumeration=None
    if args.tridagon_source:
        source=json.loads(args.tridagon_source.read_text())
        core=set(source["graph_cells"])
        selected_boxes=[box for box in BOXES if core & box]
        assert len(selected_boxes)==4 and all(len(box&core)==3 for box in selected_boxes)
        patterns=[list(tuple(sorted((box&core)|set(extra))) for extra in
                       itertools.combinations(sorted(box-core),q-3)) for box in selected_boxes]
        enumeration=iter(itertools.product(*patterns))
    positives=[];controls=[];seen=set();last_checkpoint=started
    def checkpoint():
        args.output.write_text(json.dumps({"schema_version":1,"seed":args.seed,"q":q,
          "arguments":{k:str(v) if isinstance(v,Path) else v for k,v in vars(args).items()},
          "seconds":time.monotonic()-started,"counts":dict(stats),"positives":positives,"rank_rejected_controls":controls,
          "scope":"Random q cells per box in boxes1,2,4,5; repeated supports skipped; floating negatives are discovery-only"},indent=2)+"\n")
    for _ in range(args.samples):
        now=time.monotonic()
        if now-started>=args.seconds:break
        if enumeration is None:
            groups=[rng.choice(group) for group in patterns]
        else:
            groups=next(enumeration,None)
            if groups is None:break
        cells=tuple(sorted(c for group in groups for c in group))
        stats["presentations"]+=1
        if cells in seen:stats["duplicates"]+=1;continue
        seen.add(cells)
        if any(len(set(cells)&house)>q for house in HOUSES):stats["house_capacity_rejected"]+=1;continue
        stats["house_capacity_survivors"]+=1
        independent_sets=templates(groups)
        partition,calls=integral_partition(cells,independent_sets)
        stats["partition_search_nodes"]+=calls
        if partition is not None:stats["q_colorable"]+=1;continue
        stats["non_q_colorable"]+=1
        cover=exact_fractional_cover(cells,independent_sets,q)
        if cover is None:
            stats["fractional_q_rejected"]+=1
            if len(controls)<8:controls.append({"cells":list(cells),"full_box_templates":len(independent_sets)})
        else:
            stats["exact_fractional_q_positive"]+=1
            positive={"cells":list(cells),"groups":[list(g) for g in groups],"q":q,
                      "all_full_box_independent_sets":[list(st) for st in independent_sets],
                      "fractional_cover":cover,"integral_refutation_nodes":calls,
                      "status":"Exact graph separation; respecting native completion and givens tests still required"}
            positives.append(positive);checkpoint();print("EXACT POSITIVE",json.dumps(positive),flush=True)
            if len(positives)>=args.max_positives:break
        if now-last_checkpoint>=15:
            checkpoint();print("PROGRESS",json.dumps(dict(stats)),flush=True);last_checkpoint=now
    checkpoint();print(json.dumps({"counts":dict(stats),"seconds":time.monotonic()-started,"positives":len(positives)}),flush=True)


if __name__=="__main__":main()
