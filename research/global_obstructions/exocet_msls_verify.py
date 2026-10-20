"""Standard-library acceptance for the scoped Mith515 Exocet/MSLS audit."""
import hashlib
import itertools
import json
from fractions import Fraction
from pathlib import Path

from nonrank_audit_template import mutations, verify as verify_template
from real_four_digit_mine import initial_domains
from real_verify import exact_completions,replay_trace

HERE=Path(__file__).resolve().parent
HOUSES=[set(range(9*r,9*r+9)) for r in range(9)]+[set(range(c,81,9)) for c in range(9)]
HOUSES += [set((3*br+r)*9+3*bc+c for r in range(3) for c in range(3)) for br in range(3) for bc in range(3)]

def peer(a,b):return a!=b and any(a in h and b in h for h in HOUSES)


def structural_candidates(domains,relaxed=False):
    """Different iterator: native cell pairs, then unordered cross-line triples."""
    result={}
    for a,b in itertools.combinations(range(81),2):
        if any(not domains[c]&(domains[c]-1) for c in (a,b)):continue
        ar,ac=divmod(a,9);br,bc=divmod(b,9)
        if (ar//3,ac//3)!=(br//3,bc//3):continue
        horizontal=ar==br;vertical=ac==bc
        if not(horizontal or vertical):continue
        mask=domains[a]|domains[b]
        if bin(mask).count('1')<2 or not relaxed and bin(mask).count('1') not in (3,4):continue
        index=(lambda c:c%9) if horizontal else (lambda c:c//9)
        band=(lambda c:c//27==ar//3) if horizontal else (lambda c:c%9//3==ac//3)
        for cross in itertools.combinations(range(9),3):
            if {i//3 for i in cross}!={0,1,2} or index(a) in cross or index(b) in cross:continue
            crossing={c for c in range(81) if index(c) in cross}
            targets=sorted(c for c in crossing if band(c) and domains[c]&mask and not(peer(c,a) and peer(c,b)))
            if len(targets)!=2 or not relaxed and peer(*targets) or any(not domains[c]&(domains[c]-1) for c in targets):continue
            s_cells=sorted(c for c in crossing if not band(c))
            result[((a,b),tuple(targets),tuple(s_cells))]=mask
    return result


def main():
    if not __debug__:raise SystemExit('Run without -O')
    path=HERE/'different_lead.json';case=json.loads(path.read_text());result=json.loads((HERE/'exocet_msls_results.json').read_text())
    assert result['schema_version']==1 and result['fixture_sha256']==hashlib.sha256(path.read_bytes()).hexdigest()
    assert case['target']==result['target_assumption']=={'cell':7,'digit':3,'value':True}
    assert case['conclusion']=={'cell':7,'digit':3,'value':False,'kind':'elimination'}
    assert result['selected_cells']==case['graph_cells'] and len(case['graph_cells'])==20
    initial=initial_domains([int(c) if c in '123456789' else 0 for c in case['puzzle']])
    assert initial==case['initial_domains'] and exact_completions(case['puzzle'],2)==[case['solution']]
    assert replay_trace(initial,case['base_trace'],list(map(int,case['solution'])))==case['base_domains']
    raw=case['base_domains'].copy();raw[7]=4
    assert replay_trace(raw,case['assumption_trace'])==case['assumed_domains']
    witness=json.loads((HERE/'different_lead_template_primal.json').read_text());template=verify_template(case,witness)
    rejected_mutations=mutations(case,witness)
    point=[[Fraction(0) for _ in range(9)] for _ in range(81)]
    for term in witness['terms']:
        weight=Fraction(term['numerator'],term['denominator'])
        for c in term['cells']:point[c][term['digit']-1]+=weight
    assert point[7][2]==1
    assert all(sum(v)==1 for v in point)
    assert all(sum(point[c][d] for c in h)==1 for h in HOUSES for d in range(9))
    assert all(m>0 and m&~b==0 for m,b in zip(case['assumed_domains'],case['base_domains']))
    assert all(not point[c][d] or case['base_domains'][c]>>d&1 for c in range(81) for d in range(9))
    assert sum(sum(point[c]) for c in case['graph_cells'])==20
    rejected=0
    for name,domains in [('unforced',case['base_domains']),('closed_opposite',case['assumed_domains'])]:
        expected=structural_candidates(domains,True)
        saved=result['records'][name+'_full_board']['relaxed_three_crossline_capacity']
        rejected_saved={}
        for item in saved['capacity_rejections']:
            a,b=item['bases'];targets=item['targets'];s_cells=item['s_cells'];d=item['violating_digit'];cells=item['three_native_nonpeer_cells']
            key=((a,b),tuple(targets),tuple(s_cells));assert key in expected
            assert expected[key]>>(d-1)&1 and len(cells)==len(set(cells))==3 and set(cells)<=set(s_cells)
            assert all(domains[c]>>(d-1)&1 for c in cells)
            assert all(not peer(a,b) for a,b in itertools.combinations(cells,2))
            assert key not in rejected_saved;rejected_saved[key]=item;rejected+=1
        assert set(expected)==set(rejected_saved) and saved['patterns']==[]
        assert saved['counts']['capacity_valid_patterns']==0
    # Reproduce the explicitly bounded rectangle catalogue and selected-scope scan.
    from exocet_msls_scan import junior,msls
    for name,domains in [('unforced',case['base_domains']),('closed_opposite',case['assumed_domains'])]:
        for scoped,scope in [('full_board',None),('selected20',case['graph_cells'])]:
            saved=result['records'][name+'_'+scoped]
            assert junior(domains,scope)==saved['basic_junior_exocet']
            assert junior(domains,scope,True)==saved['relaxed_three_crossline_capacity']
            assert msls(domains,scope)==saved['rectangular_msls']
            assert saved['rectangular_msls']['counts']['rank_zero_configurations']==0
            assert saved['rectangular_msls']['findings']==[]
    assert rejected==58
    control=json.loads((HERE/'exocet_msls_controls.json').read_text())
    assert control['schema_version']==1 and control['kind']=='artificial-basic-junior-positive-control'
    board=list(map(int,control['completion']));masks=control['domains']
    assert len(board)==len(masks)==81 and all({board[c] for c in h}==set(range(1,10)) for h in HOUSES)
    assert all(m>>(d-1)&1 for m,d in zip(masks,board))
    positive=junior(masks);assert len(positive['patterns'])==1
    pattern=positive['patterns'][0];assert pattern['bases']==control['bases'] and pattern['targets']==control['targets']
    assert pattern['base_digits']==control['base_digits'] and len(pattern['eliminations'])==12
    assert all(board[c]!=d for c,d in pattern['eliminations'])
    relaxed_control=masks.copy();relaxed_control[38]|=1
    assert junior(relaxed_control)['patterns']==[]
    print(json.dumps({'fixture':case['name'],'opposite_target_mass':str(point[7][2]),
      'exact_template_relaxation':template,'native_equations':324,'selected_cell_truths':20,'template_mutations_rejected':rejected_mutations,
      'relaxed_junior_geometries_rejected_by_explicit_three_packings':rejected,
      'rectangles_per_snapshot':74601,'rank_zero_rectangle_hits':0,'junior_positive_and_broken_capacity_controls':2,
      'scope':'Native rank0 cover theorem and specified basic Junior/rectangle scans only; no Senior/Weak/Complex/conditional exclusion.'}))

if __name__=='__main__':main()
