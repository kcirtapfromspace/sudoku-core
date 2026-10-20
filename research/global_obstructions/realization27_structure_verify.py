"""Standalone exact universal trade proof for the native27 house incidence.

The search enumerates ALL labelled colorings without symmetry pruning, then
independently compares them with the55 stored normalized representatives and
their24 global color permutations. No solver, numerical library or engine is
trusted. This certifies an unavoidable support, not a unique puzzle.
"""
import argparse
from collections import Counter,defaultdict
import copy
import itertools
import json
from pathlib import Path

HERE=Path(__file__).resolve().parent
PROOF=HERE/"realization27_structure_proof.json"
HOUSES=[tuple(range(9*r,9*r+9)) for r in range(9)]
HOUSES += [tuple(range(c,81,9)) for c in range(9)]
HOUSES += [tuple((3*br+r)*9+3*bc+c for r in range(3) for c in range(3)) for br in range(3) for bc in range(3)]


def sees(a,b):
    return a//9==b//9 or a%9==b%9 or (a//27,a%9//3)==(b//27,b%9//3)


def enumerate_labelled(cells):
    """Full labelled DFS: no fixed anchor, memoization or color symmetry."""
    cells=tuple(cells);n=len(cells);adj=[[j for j,b in enumerate(cells) if b!=a and sees(a,b)] for a in cells]
    assignment=[-1]*n;answers=[];nodes=0
    def walk():
        nonlocal nodes
        nodes+=1;available=[]
        for i in range(n):
            if assignment[i]!=-1:continue
            blocked={assignment[j] for j in adj[i] if assignment[j]!=-1}
            choices=[d for d in range(4) if d not in blocked]
            if not choices:return
            available.append((len(choices),-len(adj[i]),i,choices))
        if not available:answers.append(tuple(assignment));return
        _,_,i,choices=min(available)
        for d in choices:
            assignment[i]=d;walk();assignment[i]=-1
    walk()
    return set(answers),nodes


def valid_board(board):
    return len(board)==81 and all(set(board[c] for c in h)==set("123456789") for h in HOUSES)


def verify(proof=None):
    if not __debug__:raise SystemExit("Run without -O")
    p=json.loads(PROOF.read_text()) if proof is None else proof
    assert p["schema_version"]==1 and p["kind"]=="native27-universal-unavoidable-support"
    fixture=json.loads((HERE/"nonrank_graph_native27_candidate.json").read_text())
    cells=p["cells"];guardian=p["guardian"];palette=p["palette"]
    assert cells==fixture["graph_cells"] and guardian==fixture["target"]["cell"]==41
    assert palette==fixture["palette"]==[2,3,5,9] and p["guardian_extra"]==1
    colored=[c for c in cells if c!=guardian]
    assert p["colored_cells"]==colored and cells==sorted(set(cells))
    anchor=p["anchor"];assert len(anchor)==4 and set(anchor)<=set(colored)
    assert all(sees(a,b) for a,b in itertools.combinations(anchor,2))
    assert p["anchor_colors"]==[0,1,2,3]
    # The full27graph has no4coloring; hence any respecting completion must
    # give its sole extra candidate to the guardian.
    impossible,_=enumerate_labelled(cells);assert not impossible
    labelled,search_nodes=enumerate_labelled(colored)
    assert len(labelled)==p["labelled_coloring_count"]==1320
    normalized={tuple(v) for v in p["normalized_colorings"]}
    assert len(normalized)==len(p["normalized_colorings"])==p["normalized_coloring_count"]==55
    for coloring in normalized:
        assert len(coloring)==len(colored) and all(type(v) is int and 0<=v<4 for v in coloring)
        assert [coloring[colored.index(c)] for c in anchor]==[0,1,2,3]
    expanded={tuple(perm[v] for v in coloring) for coloring in normalized for perm in itertools.permutations(range(4))}
    assert expanded==labelled
    signatures=defaultdict(list)
    for coloring in labelled:
        signature=tuple(sum(1<<coloring[i] for i,c in enumerate(colored) if c in house) for house in HOUSES)
        signatures[signature].append(coloring)
    histogram=Counter(len(bucket) for bucket in signatures.values())
    assert len(signatures)==p["house_signature_count"]==408
    assert dict(sorted(histogram.items()))=={int(k):v for k,v in p["house_signature_multiplicity_histogram"].items()}=={2:192,3:96,5:72,6:48}
    assert min(histogram)==p["minimum_house_signature_multiplicity"]==2
    # Every labelled coloring therefore has a distinct same-signature partner.
    # Pasting the partner into any full completion preserves every house's
    # selected-cell digit multiset, and therefore the whole Sudoku.
    source=p["source_completion"];alternate=p["alternate_completion"]
    assert source==fixture["completion"] and valid_board(source) and valid_board(alternate)
    assert source!=alternate and source[guardian]==alternate[guardian]=="1"
    assert all(source[c]==alternate[c] for c in range(81) if c not in cells)
    differences=[c for c in range(81) if source[c]!=alternate[c]]
    assert differences==p["source_trade_cells"] and len(differences)==6
    source_coloring=tuple(palette.index(int(source[c])) for c in colored)
    alternate_coloring=tuple(palette.index(int(alternate[c])) for c in colored)
    assert source_coloring in labelled and alternate_coloring in labelled
    assert any(source_coloring in bucket and alternate_coloring in bucket for bucket in signatures.values())
    maximal="".join("." if c in cells else d for c,d in enumerate(source))
    assert maximal==p["maximal_outside_puzzle"]
    assert all(ch=="." or (ch==source[c]==alternate[c]) for c,ch in enumerate(maximal))
    source_bucket=next(bucket for bucket in signatures.values() if source_coloring in bucket)
    assert len(source_bucket)==p["maximal_outside_solution_count"]==6
    # All nonguardian initial domains of this maximal puzzle lie in P;
    # the guardian permits only P or its sole extra. Non4colorability forces
    # that extra. Thus the source signature bucket describes ALL completions.
    for cell in cells:
        allowed=set("123456789")-{source[v] for v in range(81) if v not in cells and sees(cell,v)}
        permitted=set(map(str,palette))|({str(p["guardian_extra"])} if cell==guardian else set())
        assert allowed<=permitted
    for coloring in source_bucket:
        board=list(source)
        for cell,color in zip(colored,coloring):board[cell]=str(palette[color])
        assert valid_board("".join(board))
    boxes={3*(c//27)+c%9//3 for c in colored}
    assert sorted(boxes)==p["native_boxes_touched"]==list(range(9))
    # Every possible palette clue outside support sees a full-palette cell.
    assert all(any(sees(c,s) for s in colored) for c in range(81) if c not in cells)
    return {"normalized_colorings":55,"labelled_colorings":1320,"house_signatures":408,
      "signature_multiplicities":dict(sorted(histogram.items())),"minimum_multiplicity":2,
      "universal_outside_clue_nonuniqueness":True,"stored_trade_cells":differences,
      "native_boxes_touched":9,"unpruned_labelled_search_nodes":search_nodes}


def self_test():
    source=json.loads(PROOF.read_text());edits=[lambda p:p.update(schema_version=2),
      lambda p:p["normalized_colorings"].pop(),lambda p:p["normalized_colorings"][0].__setitem__(0,4),
      lambda p:p.update(labelled_coloring_count=1321),lambda p:p.update(house_signature_count=409),
      lambda p:p["house_signature_multiplicity_histogram"].update({"1":1}),
      lambda p:p.update(alternate_completion=p["source_completion"]),lambda p:p["cells"].pop()]
    for edit in edits:
        mutated=copy.deepcopy(source);edit(mutated)
        try:verify(mutated)
        except (AssertionError,IndexError,KeyError,ValueError,TypeError):pass
        else:raise AssertionError("Malformed proof was accepted")
    return len(edits)


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument("--self-test",action="store_true");args=parser.parse_args()
    result=verify()
    if args.self_test:result["mutations_rejected"]=self_test()
    print(json.dumps(result))
