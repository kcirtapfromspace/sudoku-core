"""Persistent, provenance-bound guardian clauses for offline Sudoku discovery."""
from collections import Counter,deque
from copy import deepcopy
import hashlib
import itertools
import json
from pathlib import Path
import time

from catalog630_match import relation,proper_coloring
from catalog630_verify import native_map
from real_verify import HOUSES,sees

HERE=Path(__file__).resolve().parent
CATALOGUE_PATH=HERE/'catalog630_catalogue.json'
CATALOGUE=json.loads(CATALOGUE_PATH.read_text())
CATALOGUE_HASH=hashlib.sha256(CATALOGUE_PATH.read_bytes()).hexdigest()
PATTERNS={r['id']:r for r in CATALOGUE['patterns']}


def digest(value):
    return hashlib.sha256(json.dumps(value,sort_keys=True,separators=(',',':')).encode()).hexdigest()


def valid_domains(domains):
    if type(domains)is not list or len(domains)!=81 or any(type(m)is not int or not 0<m<=511 for m in domains):
        raise ValueError('Expected81 nonempty integer domains')


class Cutoff(Exception):pass


class Budget:
    def __init__(self,seconds=30,nodes=100000):
        self.started=time.monotonic();self.deadline=self.started+seconds
        self.limit=nodes;self.nodes=0
    def tick(self,n=1):
        self.nodes+=n
        if self.nodes>self.limit or time.monotonic()>self.deadline:
            raise Cutoff('Explicit time or node budget exhausted')


def state(domains,inherited=None,context_id=None):
    valid_domains(domains)
    if inherited is None:
        return {'version':1,'context_id':context_id or digest(domains),'catalogue_sha256':CATALOGUE_HASH,
                'origins':{},'clauses':[],'coverage':[]}
    result=deepcopy(inherited)
    if result['version']!=1 or result['catalogue_sha256']!=CATALOGUE_HASH:
        raise ValueError('Incompatible relation provenance')
    if context_id is not None and context_id!=result['context_id']:
        raise ValueError('A new puzzle must reset inherited relations')
    for key,origin in result['origins'].items():
        valid_domains(origin)
        if digest(origin)!=key or any(m&~old for m,old in zip(domains,origin)):
            raise ValueError('Inherited relation state is not a monotone ancestor')
    for clause in result['clauses']:
        validate_clause(clause,result,domains)
    return result


def validate_clause(clause,relations,domains=None):
    original=relations['origins'][clause['origin']]
    occurrence={k:clause[k]for k in ['palette','transpose','mapping','guardians']}
    native_map(PATTERNS[clause['catalogue_id']]['cells'],original,occurrence)
    if proper_coloring(PATTERNS[clause['catalogue_id']]['cells'],3)is not None:
        raise ValueError('Catalogue lemma graph is three-colorable')
    expected=digest({k:clause[k]for k in clause if k!='id'})
    if clause['id']!=expected:raise ValueError('Clause identity mismatch')
    if domains is not None and any(m&~old for m,old in zip(domains,original)):
        raise ValueError('Clause premises are not inherited monotonically')


def bounded_embeddings(pattern,host,transpose,budget):
    if len(pattern)>len(host):return
    pr={(a,b):relation(a,b,transpose)for a in pattern for b in pattern}
    hr={(a,b):relation(a,b)for a in host for b in host}
    ps={a:Counter(pr[a,b]for b in pattern if b!=a)for a in pattern}
    hs={a:Counter(hr[a,b]for b in host if b!=a)for a in host}
    available={a:[b for b in host if all(hs[b][r]>=n for r,n in ps[a].items())]for a in pattern}
    mapping,used={},set()
    def visit():
        budget.tick()
        if len(mapping)==len(pattern):yield dict(mapping);return
        choices={a:[b for b in available[a]if b not in used and all(pr[a,x]==hr[b,y]for x,y in mapping.items())]
                 for a in pattern if a not in mapping}
        a=min(choices,key=lambda c:(len(choices[c]),-sum(pr[c,z]!=0 for z in mapping),c))
        for b in choices[a]:
            mapping[a]=b;used.add(b);yield from visit();used.remove(b);del mapping[a]
    yield from visit()


def covered(relations,domains):
    return any(row['complete'] and row['catalogue_sha256']==CATALOGUE_HASH
               and all(not m&~old for m,old in zip(domains,relations['origins'][row['origin']]))
               for row in relations['coverage'])


def discover(domains,relations,budget,max_clauses=20000):
    """Strict scan BEFORE local deletion; partial discoveries stay sound at cutoff."""
    valid_domains(domains)
    if covered(relations,domains):return {'status':'complete','inherited_complete_coverage':True,'new_clauses':0}
    key=digest(domains);relations['origins'][key]=domains.copy()
    ids={c['id']for c in relations['clauses']};seen=set();before=len(ids);queries=0
    minimum=min(len(p['cells'])for p in PATTERNS.values())
    status='complete'
    try:
        for digits in itertools.combinations(range(1,10),3):
            mask=sum(1<<(d-1)for d in digits)
            host=[c for c,m in enumerate(domains)if m&mask==mask]
            if len(host)<minimum:continue
            for pattern in PATTERNS.values():
                for transpose in (False,True):
                    budget.tick();queries+=1
                    for mapping in bounded_embeddings(pattern['cells'],host,transpose,budget):
                        signature=(pattern['id'],digits,tuple(sorted(mapping.values())))
                        if signature in seen:continue
                        seen.add(signature)
                        guardians=sorted([c,d]for c in mapping.values()for d in range(1,10)
                                         if domains[c]&(1<<(d-1)) and d not in digits)
                        clause={'origin':key,'catalogue_id':pattern['id'],'palette':list(digits),
                                'transpose':transpose,'mapping':sorted([a,b]for a,b in mapping.items()),'guardians':guardians}
                        clause['id']=digest(clause)
                        if clause['id']not in ids:
                            relations['clauses'].append(clause);ids.add(clause['id'])
                        if len(ids)>max_clauses:raise Cutoff('Guardian-clause storage budget exhausted')
    except Cutoff:
        status='cutoff'
    row={'origin':key,'catalogue_sha256':CATALOGUE_HASH,'complete':status=='complete',
         'orientation_queries':queries,'embedding_nodes':budget.nodes,'new_clauses':len(ids)-before}
    relations['coverage'].append(row)
    return {'status':status,**row}


def clause_status(clause,domains):
    remaining=[]
    for cell,digit in clause['guardians']:
        bit=1<<(digit-1)
        if domains[cell]==bit:return 'satisfied',[]
        if domains[cell]&bit:remaining.append([cell,digit])
    return ('false'if not remaining else 'unit'if len(remaining)==1 else 'open'),remaining


def unit_step(domains,relations):
    for clause in relations['clauses']:
        status,remaining=clause_status(clause,domains)
        if status=='false':return {'kind':'guardian-contradiction','clause':clause['id']}
        if status=='unit':
            cell,digit=remaining[0]
            return {'kind':'guardian-unit','clause':clause['id'],'cell':cell,'digit':digit,'value':True}
    return None


def implication_graph(domains,relations):
    """Native conflict edges and binary positive ORs; endpoints encode9*cell+d (1..729)."""
    graph={}
    def edge(a,b,reason):graph.setdefault(a,[]).append((b,reason))
    def positive(literals,reason):
        if len(literals)==2:
            a,b=literals;edge(-a,b,reason);edge(-b,a,reason)
    available=[(c,d,9*c+d)for c,m in enumerate(domains)for d in range(1,10)if m&(1<<(d-1))]
    for c,m in enumerate(domains):
        positive([9*c+d for d in range(1,10)if m&(1<<(d-1))],{'kind':'cell-or','cell':c})
    for h,cells in enumerate(HOUSES):
        for d in range(1,10):positive([9*c+d for c in cells if domains[c]&(1<<(d-1))],{'kind':'house-or','house':h,'digit':d})
    for i,(a,d,x)in enumerate(available):
        for b,e,y in available[i+1:]:
            if a==b or d==e and sees(a,b):
                reason={'kind':'native-conflict','a':[a,d],'b':[b,e]}
                edge(x,-y,reason);edge(y,-x,reason)
    for clause in relations['clauses']:
        status,left=clause_status(clause,domains)
        if status=='open'and len(left)==2:positive([9*c+d for c,d in left],{'kind':'guardian-or','clause':clause['id']})
    return graph,available


def binary_step(domains,relations,budget):
    graph,available=implication_graph(domains,relations)
    for cell,digit,x in available:
        if domains[cell]&(domains[cell]-1)==0:continue
        for start,goal,value in [(x,-x,False),(-x,x,True)]:
            parent={start:None};queue=deque([start])
            while queue and goal not in parent:
                budget.tick();v=queue.popleft()
                for w,reason in graph.get(v,[]):
                    if w not in parent:parent[w]=(v,reason);queue.append(w)
            if goal in parent:
                path=[];v=goal
                while v!=start:
                    u,reason=parent[v];path.append({'from':u,'to':v,'reason':reason});v=u
                return {'kind':'binary-forcing','cell':cell,'digit':digit,'value':value,'path':list(reversed(path))}
    return None
