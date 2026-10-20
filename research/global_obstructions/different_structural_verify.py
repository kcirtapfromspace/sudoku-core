"""Standard-library checks for the bounded different-geometry negative controls.

The search uses SciPy and NetworkX. This checker independently establishes each
retained rank rejection by maximizing its rational weight over ALL native
independent sets. It does not rely on a floating solver or the search's list of
maximal independent sets.
"""
import argparse
from collections import Counter
from fractions import Fraction
from functools import lru_cache
from math import gcd
import json
from pathlib import Path

HERE=Path(__file__).resolve().parent


def sees(a,b):
    return a//9==b//9 or a%9==b%9 or (a//27,a%9//3)==(b//27,b%9//3)


def rank_check(record):
    cells=record['cells']; assert cells==sorted(set(cells))
    assert all(isinstance(c,int) and 0<=c<81 for c in cells)
    weights=[Fraction(*v) for v in record['dual_vertex_weights']]
    assert len(weights)==len(cells) and all(v>=0 for v in weights)
    denominator=1
    for v in weights:denominator=denominator*v.denominator//gcd(denominator,v.denominator)
    integers=[int(v*denominator) for v in weights]
    adjacency=[sum(1<<j for j,b in enumerate(cells) if i!=j and sees(a,b)) for i,a in enumerate(cells)]
    @lru_cache(None)
    def maximum(available):
        if not available:return 0
        vertices=[i for i in range(len(cells)) if available>>i&1]
        i=max(vertices,key=lambda j:((adjacency[j]&available).bit_count() if hasattr(int,'bit_count') else bin(adjacency[j]&available).count('1'),integers[j],-j))
        remaining=available&~(1<<i)
        return max(maximum(remaining),integers[i]+maximum(remaining&~adjacency[i]))
    optimum=maximum((1<<len(cells))-1)
    assert optimum<=denominator, 'Invalid weighted independence bound'
    assert sum(integers)>4*denominator, 'Not a four-color rank rejection'
    assert len({3*(c//27)+c%9//3 for c in cells})<=8
    return maximum.cache_info().currsize


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--summary',type=Path,default=HERE/'different_structural_summary.json')
    args=parser.parse_args();summary=json.loads(args.summary.read_text())
    assert summary['schema_version']==1
    assert summary['kind']=='bounded-omitted-box-geometry-search'
    counts=summary['counts'];controls=summary['rank_controls']
    assert len(controls)==counts['exact_weighted_rank_explained']==177
    assert len({tuple(q['cells']) for q in controls})==len(controls)
    assert counts['critical_subgraphs']==len(controls)+counts['duplicate_graph_fingerprints']
    assert counts['critical_subgraphs']==6*counts['non_four_colorable_supports']
    assert counts['forced_box_omissions']==counts['new_support_presentations']
    assert summary['pure_graph_nonrank_hits']==0
    histogram={str(n):v for n,v in Counter(len(q['cells']) for q in controls).items()}
    assert histogram==summary['critical_order_histogram']
    nodes=sum(rank_check(record) for record in controls)
    for control in summary['coordinate_rejections']:
        clique=control['native_five_clique']
        assert len(clique)==len(set(clique))==5
        assert all(sees(a,b) for i,a in enumerate(clique) for b in clique[i+1:])
        assert set(clique)<=set(control['cells'])
    from different_structural_avoidability import screen
    fixture=json.loads((HERE/'nonrank_graph_native27_candidate.json').read_text())
    import itertools
    coordinate_expected={}
    for permutation in itertools.permutations(range(4)):
        image=[]
        for c in fixture['graph_cells']:
            row,column=divmod(c,9);digits=[row//3,row%3,column//3,column%3]
            a,b,d,e=[digits[i] for i in permutation];image.append(9*(3*a+b)+3*d+e)
        image=sorted(image)
        if len({3*(c//27)+c%9//3 for c in image})<=8:coordinate_expected[permutation]=image
    assert len(coordinate_expected)==8
    assert {tuple(q['coordinate_permutation']):q['cells'] for q in summary['coordinate_rejections']}==coordinate_expected
    result=screen(fixture['graph_cells'],fixture['palette'],fixture['completion'])
    assert result['source_completion_trade_found'] and not result['necessary_avoidability_screen_passed']
    assert result['outside_clue_nonuniqueness_witness']!=fixture['completion']
    single=screen([0],[int(fixture['completion'][0])],fixture['completion'])
    assert single['necessary_avoidability_screen_passed'] and not single['source_completion_trade_found']
    malformed=json.loads(json.dumps(summary['rank_controls'][0]));malformed['dual_vertex_weights'][0]=[100,1]
    try:rank_check(malformed)
    except AssertionError:pass
    else:raise AssertionError('Malformed rank witness accepted')
    print(json.dumps({'rank_controls':len(summary['rank_controls']),'independent_set_dp_states':nodes,
      'native_five_clique_controls':len(summary['coordinate_rejections']),
      'source_signature_trade_and_rigid_controls':2,'malformed_rank_rejected':True}))


if __name__=='__main__':main()
