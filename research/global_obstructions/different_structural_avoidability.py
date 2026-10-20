"""Exact early trade screen for a proposed support and respecting completion.

Only palette-valued support cells are changed; all guardians keep their actual
completion values. We enumerate the source27-house signature, stopping after
two colorings. Two are a decisive internal trade; one is merely a necessary
screen, not a universal classification of other colorings or embeddings.
"""
import argparse
import json
from pathlib import Path

HOUSES=[set(range(9*r,9*r+9)) for r in range(9)]
HOUSES += [set(range(c,81,9)) for c in range(9)]
HOUSES += [set((3*br+r)*9+3*bc+c for r in range(3) for c in range(3)) for br in range(3) for bc in range(3)]


def sees(a,b):
    return a//9==b//9 or a%9==b%9 or (a//27,a%9//3)==(b//27,b%9//3)


def screen(cells,palette,completion,node_limit=1000000):
    assert len(completion)==81 and all({int(completion[c]) for c in h}==set(range(1,10)) for h in HOUSES)
    assert cells==sorted(set(cells)) and all(0<=c<81 for c in cells)
    assert palette==sorted(set(palette)) and all(1<=d<=9 for d in palette)
    colors=[c for c in cells if int(completion[c]) in palette]
    house_masks=[sum(1<<(int(completion[c])-1) for c in colors if c in h) for h in HOUSES]
    allowed={c:(1<<9)-1 for c in colors}
    for c in colors:
        for h,mask in zip(HOUSES,house_masks):
            if c in h:allowed[c]&=mask
    adjacency={c:[v for v in colors if v!=c and sees(c,v)] for c in colors}
    assigned={};answers=[];nodes=0
    def search():
        nonlocal nodes
        nodes+=1
        if nodes>node_limit:raise RuntimeError("Node cutoff: no unique-signature conclusion")
        if len(answers)>=2:return
        candidates={c:allowed[c]&~sum(set(assigned.get(v,0) for v in adjacency[c])) for c in colors if c not in assigned}
        if not candidates:answers.append(dict(assigned));return
        c=min(candidates,key=lambda c:(bin(candidates[c]).count("1"),-len(adjacency[c]),c));choices=candidates[c]
        while choices and len(answers)<2:
            bit=choices&-choices;choices-=bit;assigned[c]=bit;search();del assigned[c]
    search();assert answers
    partner=None
    for answer in answers:
        board=list(completion)
        for c,bit in answer.items():board[c]=str(bit.bit_length())
        board="".join(board)
        assert all({int(board[c]) for c in h}==set(range(1,10)) for h in HOUSES)
        assert all(board[c]==completion[c] for c in range(81) if c not in cells)
        if board!=completion:partner=board
    if len(answers)==2:assert partner is not None
    return {"source_house_signature_solutions_capped_2":len(answers),"nodes":nodes,
      "outside_clue_nonuniqueness_witness":partner,
      "source_completion_trade_found":partner is not None,
      "necessary_avoidability_screen_passed":len(answers)==1,
      "scope":"Source completion signature only; no claim for all compatible completions"}


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument("fixture",type=Path);parser.add_argument("--node-limit",type=int,default=1000000);args=parser.parse_args()
    case=json.loads(args.fixture.read_text());cells=case.get("graph_cells",case.get("cells"));completion=case.get("completion",case.get("solution"));palette=case["palette"]
    print(json.dumps(screen(cells,palette,completion,args.node_limit)))
