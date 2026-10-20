"""Exact independent audit of the pure graph gap and its local rejection.

This verifies a structural separation, not a new Sudoku rule or a unique
givens/full-old-engine-stall benchmark. The 27-cell state is rejected because
existing local propagation already places its guardian.
"""
from fractions import Fraction
import hashlib
import itertools
import json
from pathlib import Path

from expanded_motif_verify import HOUSES, maximal_independent_sets, sees
from real_verify import replay_trace, triple_feasible
from nonrank_graph_certificate import check as check_research_tree

HERE=Path(__file__).resolve().parent


def solve(cells, domains, dynamic=False):
    cells=tuple(sorted(cells));adj={c:[v for v in cells if v!=c and sees(c,v)] for c in cells}
    assigned={};nodes=0;leaves=0
    def recurse():
        nonlocal nodes,leaves
        nodes+=1
        available={c:domains[c]&~sum(set(assigned.get(v,0) for v in adj[c]))
                   for c in cells if c not in assigned}
        if not available:return dict(assigned)
        def degree(c):return sum(v not in assigned for v in adj[c]) if dynamic else len(adj[c])
        c=min(available,key=lambda c:(bin(available[c]).count("1"),-degree(c),c))
        choices=available[c]
        if not choices:leaves+=1
        while choices:
            bit=choices&-choices;choices-=bit;assigned[c]=bit
            answer=recurse()
            if answer is not None:return answer
            del assigned[c]
        return None
    answer=recurse()
    return answer,nodes,leaves


def verify():
    if not __debug__:raise SystemExit("Run without -O; this exact verifier uses assertions")
    proof=json.loads((HERE/"nonrank_graph_27_proof.json").read_text())["pure_graph_gap"]
    cells=proof["cells"];context=proof["context"];completion=list(map(int,context["completion"]))
    cover=[(entry["cells"],Fraction(*entry["weight"])) for entry in proof["fractional_cover"]]
    assert all(w>=0 and set(st)<=set(cells) and all(not sees(a,b) for a,b in itertools.combinations(st,2)) for st,w in cover)
    assert sum(w for st,w in cover)==4
    assert all(sum(w for st,w in cover if c in st)==1 for c in cells)
    clique=[c for c,w in zip(cells,proof["dual_vertex_weights"]) if Fraction(*w)==1]
    assert len(clique)==4 and all(sees(a,b) for a,b in itertools.combinations(clique,2))
    assert all({completion[c] for c in house}==set(range(1,10)) for house in HOUSES)
    assert all(completion[c] in context["palette"] for c in cells if c!=context["guardian"])
    assert completion[context["guardian"]]==context["extra"] not in context["palette"]
    answer,nodes,leaves=solve(cells,{c:15 for c in cells},True)
    assert answer is None
    assert nodes==proof["independent_full_branch_refutation"]["nodes"]==5081
    assert leaves==proof["independent_full_branch_refutation"]["empty_domain_leaves"]==504
    assert all(solve(set(cells)-{c},{v:15 for v in cells},True)[0] is not None for c in cells)
    candidate=json.loads((HERE/"nonrank_graph_27_candidate.json").read_text())
    base=replay_trace(candidate["initial_domains"],candidate["base_trace"],completion)
    assert base==candidate["base_domains"]
    assert base[context["guardian"]]==1<<(context["extra"]-1)
    controls=json.loads((HERE/"nonrank_graph_rejected_controls.json").read_text())["controls"]
    for control in controls:
        vertices=control["cells"];weights=[Fraction(*pair) for pair in control["dual_vertex_weights"]]
        independent=maximal_independent_sets(vertices)
        assert len(independent)==control["maximal_independent_sets_count"]
        assert all(w>=0 for w in weights)
        by_cell=dict(zip(vertices,weights))
        assert all(sum(by_cell[c] for c in st)<=1 for st in independent)
        assert sum(weights)==Fraction(*control["exact_weighted_capacity_lower_bound"])>4
    modern=json.loads((HERE/"nonrank_graph_modern_classification.json").read_text())
    coloring={int(c):d for c,d in modern["uniform_four_coloring"].items()}
    assert set(coloring)==set(modern["selected_cells"]) and all(1<=v<=4 for v in coloring.values())
    assert all(coloring[a]!=coloring[b] for a,b in itertools.combinations(coloring,2) if sees(a,b))
    closed_core=modern["closed_tridagon_core"]
    assert solve(closed_core,dict(zip(closed_core,modern["closed_core_domains"])))[0] is None
    mapping={int(c):v for c,v in modern["known_tridagon_mapping"].items()}
    assert set(mapping)==set(closed_core) and len(set(mapping.values()))==12
    assert all(sees(a,b)==sees(mapping[a],mapping[b]) for a,b in itertools.combinations(closed_core,2))
    masks=dict(zip(modern["selected_cells"],modern["base_core_opposite_domains"]))
    minimal=modern["minimum_found_base_opposite_cell_cores"][0]
    assert len(minimal)==14 and solve(minimal,masks)[0] is None
    digest=hashlib.sha256();count=0
    for subset in itertools.combinations(modern["selected_cells"],13):
        assignment,_,_=solve(subset,masks);assert assignment is not None
        digest.update(json.dumps([[c,assignment[c].bit_length()] for c in subset],separators=(",",":")).encode());count+=1
    assert count==560 and digest.hexdigest()==modern["thirteen_subset_witness_digest"]
    second=json.loads((HERE/"nonrank_graph_native27_proof.json").read_text())["pure_graph_gap"]
    vertices=second["cells"]
    second_cover=[(entry["cells"],Fraction(*entry["weight"])) for entry in second["fractional_cover"]]
    assert sum(w for st,w in second_cover)==4
    assert all(w>=0 and set(st)<=set(vertices) and all(not sees(a,b) for a,b in itertools.combinations(st,2)) for st,w in second_cover)
    assert all(sum(w for st,w in second_cover if c in st)==1 for c in vertices)
    second_clique=[c for c,w in zip(vertices,second["dual_vertex_weights"]) if Fraction(*w)==1]
    assert len(second_clique)==4 and all(sees(a,b) for a,b in itertools.combinations(second_clique,2))
    answer,second_nodes,second_leaves=solve(vertices,{c:15 for c in vertices},True)
    assert answer is None and second_nodes==16745 and second_leaves==2016
    assert all(solve(set(vertices)-{c},{v:15 for v in vertices},True)[0] is not None for c in vertices)
    second_case=json.loads((HERE/"nonrank_graph_native27_candidate.json").read_text())
    sol=list(map(int,second_case["completion"]))
    assert all({sol[c] for c in h}==set(range(1,10)) for h in HOUSES)
    second_base=replay_trace(second_case["initial_domains"],second_case["base_trace"],sol)
    assert second_base==second_case["base_domains"]
    assumed=second_base.copy();assumed[second_case["target"]["cell"]]&=~(1<<(second_case["target"]["digit"]-1))
    second_branch=replay_trace(assumed,second_case["assumption_trace"])
    assert second_branch==second_case["assumed_domains"]
    assert all(triple_feasible(second_branch,houses) for houses in itertools.combinations(range(27),3))
    role=second_case["role_witness"];palette=second_case["palette"]
    point=[[Fraction(role[c],4) if d in palette else Fraction(1-role[c],5) for d in range(1,10)] for c in range(81)]
    assert all(0<=v<=1 for row in point for v in row)
    assert all(sum(row)==1 for row in point)
    assert all(sum(point[c][d] for c in h)==1 for h in HOUSES for d in range(9))
    assert all(second_branch[c]&(1<<d) or point[c][d]==0 for c in range(81) for d in range(9))
    assert second_base[second_case["target"]["cell"]]&(second_base[second_case["target"]["cell"]]-1)
    tree_check=check_research_tree()
    print(json.dumps({"pure_graph_chromatic_number":5,"pure_graph_fractional_chromatic_number":4,
      "full_branch_nodes":nodes,"locally_rejected_guardian_placement":True,
      "exact_weighted_rejection_controls":len(controls),"modern_known_tridagon_minimum_within_selected_cells":14,
      "second_27_native_equations_verified":324,"second_27_refutation_nodes":second_nodes,
      "second_27_hall81_traces_verified":True,"second_27_all_house_triples_feasible":2925,
      "second_27_full_template_status":"49term primal verified by separate nonrank_audit_template checker",
      "research_finite_tree_nodes":tree_check["nodes"],"second_27_vertex_critical":True,
      "unique_puzzle_new_family_acceptance":False}))


if __name__=="__main__":verify()
