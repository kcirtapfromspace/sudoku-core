"""Build/check the research-only 27-cell exhaustive coloring refutation tree.

Every allowed branch is present. Checker geometry is reconstructed from native
coordinates; it does not trust the search heuristic, claimed chromatic number,
fractional optimum, conflict label, or any omitted branch. Production limits
are intentionally unchanged.
"""
import argparse
import json
from pathlib import Path
import random
import copy

def sees(a,b):
    return a//9==b//9 or a%9==b%9 or (a//27,a%9//3)==(b//27,b%9//3)

HERE=Path(__file__).resolve().parent
PATH=HERE/"nonrank_graph_native27_tree.json"


def build(trial=0):
    case=json.loads((HERE/"nonrank_graph_native27_candidate.json").read_text())
    cells=case["graph_cells"];n=len(cells);adj=[sum(1<<j for j,b in enumerate(cells) if a!=b and sees(a,b)) for a in cells]
    ds=[case["base_domains"][c] for c in cells];target=case["target"];ds[cells.index(target["cell"])]&=~(1<<(target["digit"]-1))
    rng=random.Random(991+trial);tie=list(range(n));rng.shuffle(tie);rank={j:i for i,j in enumerate(tie)};nodes=[]
    def search(domains,assigned):
        if assigned==(1<<n)-1:raise RuntimeError("Satisfying coloring: no refutation certificate")
        def key(j):
            d=bin(domains[j]).count("1");ad=adj[j]&~assigned
            if trial%3==0:s=-bin(ad).count("1")
            elif trial%3==1:s=-bin(adj[j]).count("1")
            else:s=-sum(1/(bin(domains[k]).count("1") or .1) for k in range(n) if ad>>k&1)
            return d,s,rank[j]
        j=min((i for i in range(n) if not assigned>>i&1),key=key);choices=domains[j]
        if not choices:
            nodes.append(["conflict",cells[j]]);return len(nodes)-1
        branches=[]
        while choices:
            bit=choices&-choices;choices-=bit;updated=domains.copy();updated[j]=bit
            for k in range(n):
                if adj[j]>>k&1 and not assigned>>k&1:updated[k]&=~bit
            child=search(updated,assigned|1<<j);branches.append([bit.bit_length(),child])
        nodes.append(["split",cells[j],branches]);return len(nodes)-1
    root=search(ds,0)
    certificate={"schema_version":1,"kind":"research-only-native-list-coloring-refutation",
       "cells":cells,"base_domains":[case["base_domains"][c] for c in cells],
       "base_witness":[int(case["completion"][c]) for c in cells],
       "target":{"cell":target["cell"],"digit":target["digit"],"value":True},
       "root":root,"nodes":nodes,"search_trial":trial,
       "scope":"27cells exceeds production24celllimit; no uniqueness or historical novelty claim"}
    PATH.write_text(json.dumps(certificate,separators=(",",":"))+"\n")
    return certificate


def check(certificate=None):
    if not __debug__:raise SystemExit("Run without -O")
    c=json.loads(PATH.read_text()) if certificate is None else certificate
    assert isinstance(c,dict) and c["schema_version"]==1 and c["kind"]=="research-only-native-list-coloring-refutation"
    cells=c["cells"];assert isinstance(cells,list) and cells==sorted(set(cells))
    assert len(cells)==27 and all(type(x) is int and 0<=x<81 for x in cells)
    nodes=c["nodes"];assert isinstance(nodes,list) and 0<len(nodes)<=100000 and type(c["root"]) is int and c["root"]==len(nodes)-1
    assert len(c["base_domains"])==len(c["base_witness"])==len(cells)
    domains=dict(zip(cells,c["base_domains"]));assert all(0<m<512 for m in domains.values())
    target=c["target"];assert target["value"] is True and target["cell"] in cells
    case=json.loads((HERE/"nonrank_graph_native27_candidate.json").read_text())
    assert cells==case["graph_cells"]
    assert target=={"cell":case["target"]["cell"],"digit":case["target"]["digit"],"value":True}
    assert type(target["cell"]) is int and type(target["digit"]) is int
    assert c["base_domains"]==[case["base_domains"][cell] for cell in cells]
    assert all(type(d) is int and 1<=d<=9 and domains[cell]&(1<<(d-1)) for cell,d in zip(cells,c["base_witness"]))
    witness=dict(zip(cells,c["base_witness"]))
    assert witness[target["cell"]]==target["digit"]
    assert all(witness[a]!=witness[b] for a in cells for b in cells if a<b and sees(a,b))
    domains[target["cell"]]&=~(1<<(target["digit"]-1));visited=set();conflicts=0
    def walk(index,ds,assigned):
        nonlocal conflicts
        assert index not in visited;visited.add(index);node=nodes[index];cell=node[1]
        assert cell in ds and cell not in assigned
        if node[0]=="conflict":
            assert len(node)==2 and ds[cell]==0;conflicts+=1;return
        assert node[0]=="split" and len(node)==3
        branches=node[2];assert len(branches)==bin(ds[cell]).count("1")
        assert {d for d,child in branches}=={d for d in range(1,10) if ds[cell]&(1<<(d-1))}
        for digit,child in branches:
            assert 0<=child<index;updated=ds.copy();bit=1<<(digit-1);updated[cell]=bit
            for peer in cells:
                if peer!=cell and peer not in assigned and sees(cell,peer):updated[peer]&=~bit
            walk(child,updated,assigned|{cell})
    walk(c["root"],domains,set());assert visited==set(range(len(nodes)))
    return {"cells":len(cells),"nodes":len(nodes),"conflicts":conflicts,"complete_branch_coverage":True,"production_limits_unchanged":True}


def self_test():
    original=json.loads(PATH.read_text());mutations=[]
    for label,edit in [
        ("schema",lambda c:c.update(schema_version=2)),
        ("kind",lambda c:c.update(kind="untrusted")),
        ("target-digit",lambda c:c["target"].update(digit=2)),
        ("target-polarity",lambda c:c["target"].update(value=False)),
        ("short-witness",lambda c:c["base_witness"].pop()),
        ("witness-target",lambda c:c["base_witness"].__setitem__(c["cells"].index(c["target"]["cell"]),2)),
        ("duplicate-cell",lambda c:c["cells"].__setitem__(0,c["cells"][1])),
        ("root",lambda c:c.update(root=0)),
        ("missing-branch",lambda c:c["nodes"][c["root"]][2].pop()),
        ("cycle",lambda c:c["nodes"][c["root"]][2][0].__setitem__(1,c["root"])),
        ("illegal-digit",lambda c:c["nodes"][c["root"]][2][0].__setitem__(0,0)),
        ("premise",lambda c:c["base_domains"].__setitem__(0,511)),
    ]:
        mutated=copy.deepcopy(original);edit(mutated)
        try:check(mutated)
        except (AssertionError,KeyError,IndexError,TypeError,ValueError):mutations.append(label)
        else:raise AssertionError("Accepted malformed certificate: "+label)
    return len(mutations)


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument("--build",action="store_true");parser.add_argument("--trial",type=int,default=0);parser.add_argument("--self-test",action="store_true");args=parser.parse_args()
    result=check(build(args.trial) if args.build else None)
    if args.self_test:result["malformed_mutations_rejected"]=self_test()
    print(json.dumps(result))
