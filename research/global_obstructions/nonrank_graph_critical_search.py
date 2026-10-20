"""Randomize critical subgraphs before rejecting a weighted-rank explanation.

A rank-explained large support can contain a different critical subgraph that
has a genuine coloring/fractional-coloring gap. Therefore one deterministic
vertex-deletion order is insufficient discovery coverage. This script tries
many orders in supports anchored to different respecting Sudoku completions.
Positive fractional covers pass exact rational checks; floating rejections are
only exploratory until their rational dual has been independently checked.
"""
import os
for _name in ("OMP_NUM_THREADS","OPENBLAS_NUM_THREADS","MKL_NUM_THREADS"):
    os.environ[_name]="1"
import argparse
from collections import Counter
from fractions import Fraction
from functools import lru_cache
import itertools
import json
from pathlib import Path
import random
import time

import networkx as nx
import numpy as np
from scipy.optimize import linprog
from structural_discover import HOUSES, PEERS, four_colorable
from real_verify import exact_completions

HERE=Path(__file__).resolve().parent
ROOT=HERE.parent.parent


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--seconds",type=int,default=420)
    parser.add_argument("--orders",type=int,default=30)
    parser.add_argument("--seed",type=int,default=219440)
    parser.add_argument("--output",type=Path,default=HERE/"nonrank_graph_critical_results.json")
    parser.add_argument("--require-native-lp",action="store_true",help="Continue past pure graph gaps whose symmetric native board LP is infeasible")
    args=parser.parse_args();rng=random.Random(args.seed);started=time.monotonic();stats=Counter();classified={};positives=[];pure_gaps=[]
    coefficients=np.zeros((27,81),dtype=np.int8)
    for i,h in enumerate(HOUSES):
        for c in h:coefficients[i,c]=1
    def expired():return time.monotonic()-started>=args.seconds
    def checkpoint():
        args.output.write_text(json.dumps({"schema_version":1,"seed":args.seed,"seconds":time.monotonic()-started,
          "arguments":{k:str(v) if isinstance(v,Path) else v for k,v in vars(args).items()},"counts":dict(stats),
          "classified_graphs":list(classified.values()),"positives":positives,"pure_graph_gaps":pure_gaps,
          "scope":"Known and new completion+guardian supports, randomized minimal critical subgraphs; discovery negatives are not a classification"},indent=2)+"\n")
    def inspect(vertices,context):
        cache={}
        def colorable(cells):
            key=tuple(sorted(cells))
            if key not in cache:
                cache[key]=four_colorable(key,node_limit=1000000)[0]
                stats["exact_coloring_queries"]+=1
            return cache[key]
        if colorable(vertices):return
        stats["non_four_colorable_supports"]+=1
        for _ in range(args.orders):
            if expired():return
            remaining=set(vertices);order=sorted(vertices);rng.shuffle(order)
            for c in order:
                if not colorable(remaining-{c}):remaining.remove(c)
            stats["critical_subgraphs"]+=1
            cells=sorted(remaining);graph=nx.Graph();graph.add_nodes_from(cells)
            graph.add_edges_from((a,b) for a,b in itertools.combinations(cells,2) if b in PEERS[a])
            fingerprint=(len(cells),graph.number_of_edges(),nx.weisfeiler_lehman_graph_hash(graph,iterations=5))
            if fingerprint in classified:
                stats["duplicate_graph_fingerprints"]+=1;continue
            sets=sorted(sorted(st) for st in nx.find_cliques(nx.complement(graph)))
            matrix=np.array([[int(c in st) for st in sets] for c in cells],dtype=np.int8)
            answer=linprog(np.ones(len(sets)),A_ub=-matrix,b_ub=-np.ones(len(cells)),bounds=(0,None),method="highs-ds")
            assert answer.success
            dual=[Fraction(float(-v)).limit_denominator(1000000) for v in answer.ineqlin.marginals]
            assert all(v>=0 for v in dual)
            assert all(sum(dual[i] for i,c in enumerate(cells) if c in st)<=1 for st in sets)
            primal=[Fraction(float(v)).limit_denominator(1000000) for v in answer.x]
            assert all(v>=0 for v in primal)
            assert all(sum(v for st,v in zip(sets,primal) if c in st)>=1 for c in cells)
            assert sum(primal)==sum(dual)
            record={"cells":cells,"native_edges":graph.number_of_edges(),"alpha":max(map(len,sets)),
              "wl_hash":fingerprint[-1],"exact_fractional_chromatic":[sum(primal).numerator,sum(primal).denominator],
              "dual_vertex_weights":[[v.numerator,v.denominator] for v in dual],"context":context}
            classified[fingerprint]=record
            if sum(primal)>4:
                stats["exact_weighted_rank_explained"]+=1
            else:
                stats["exact_nonrank_positive"]+=1
                record["fractional_cover"]=[{"cells":st,"weight":[v.numerator,v.denominator]} for st,v in zip(sets,primal) if v]
                record["status"]="Exact pure native graph gap; domain closure/full template LP/givens still to verify"
                role=linprog(np.zeros(81),A_eq=coefficients,b_eq=np.full(27,4),bounds=[(1,1) if c in cells else (0,1) for c in range(81)],method="highs-ds")
                record["floating_native_role_lp_feasible"]=bool(role.success)
                pure_gaps.append(record)
                if not role.success:
                    stats["native_role_lp_rejected_pure_gap"]+=1
                    if args.require_native_lp:
                        checkpoint();print("PURE GAP REJECTED BY NATIVE LP",len(cells),flush=True)
                        continue
                else:
                    point=[Fraction(float(v)).limit_denominator(1000000) for v in role.x]
                    assert all(0<=v<=1 for v in point) and all(point[c]==1 for c in cells)
                    assert all(sum(point[c] for c in h)==4 for h in HOUSES)
                    record["exact_native_role_witness"]=[[v.numerator,v.denominator] for v in point]
                    stats["exact_native_role_positive_pure_gap"]+=1
                positives.append(record);checkpoint();print("EXACT NONRANK POSITIVE",json.dumps(record),flush=True)
                return
        checkpoint()
    # Revisit large original supports, not just the single old minimized motif.
    for path in [HERE/"expanded_motif_results.json",HERE/"expanded_motif_mixed_results.json"]:
        for old in json.loads(path.read_text())["motifs"]:
            if expired() or positives:break
            completion=list(map(int,old["completion"]));palette=old["palette"];guardian=old["guardian"]
            vertices={c for c,d in enumerate(completion) if d in palette}-set(old["removed_completion_cells"])|{guardian}
            stats["seeded_supports"]+=1
            inspect(vertices,{"source_id":old["source_id"],"completion":old["completion"],"palette":palette,"guardian":guardian,"extra":completion[guardian]})
        if expired() or positives:break
    corpus=json.loads((ROOT/"research/nonmodular/real_game_corpus.json").read_text());rng.shuffle(corpus)
    for source in corpus:
        if expired() or positives:break
        completion=exact_completions(source["puzzle"],1)[0];digits=list(map(int,completion));stats["new_completed_grids"]+=1
        palettes=list(itertools.combinations(range(1,10),4));rng.shuffle(palettes)
        print("GRID",source["id"],"stats",json.dumps(dict(stats)),flush=True)
        for palette in palettes:
            base={c for c,d in enumerate(digits) if d in palette};guardians=[c for c,d in enumerate(digits) if d not in palette];rng.shuffle(guardians)
            for guardian in guardians[:10]:
                choices=[sorted(base&h) for h in HOUSES if guardian in h];removals=list(itertools.product(*choices));rng.shuffle(removals)
                for removed in removals:
                    if expired() or positives:break
                    vertices=base-set(removed)|{guardian};stats["new_support_presentations"]+=1
                    inspect(vertices,{"source_id":source["id"],"completion":completion,"palette":list(palette),"guardian":guardian,"extra":digits[guardian]})
                if expired() or positives:break
            if expired() or positives:break
        checkpoint()
    checkpoint();print(json.dumps({"counts":dict(stats),"classified_graphs":len(classified),"positives":len(positives),"seconds":time.monotonic()-started}),flush=True)


if __name__=="__main__":main()
