"""Bounded native MSLS and basic Junior Exocet audit of the fixed Mith515 lead.

Independent implementation of explicitly specified variants, not a run of the
external .NET library. No absence result covers Senior/Weak/Complex Exocets.
"""
import argparse
from collections import Counter
from functools import lru_cache
import hashlib
import itertools
import json
from pathlib import Path

HERE=Path(__file__).resolve().parent
HOUSES=[set(range(9*r,9*r+9)) for r in range(9)]+[set(range(c,81,9)) for c in range(9)]
HOUSES += [set((3*br+r)*9+3*bc+c for r in range(3) for c in range(3)) for br in range(3) for bc in range(3)]
def peer(a,b):return a!=b and (a//9==b//9 or a%9==b%9 or (a//27,a%9//3)==(b//27,b%9//3))
def count(m):return bin(m).count('1')
def digits(m):return [d for d in range(1,10) if m>>(d-1)&1]


def msls(domains,scope=None):
    """Author-library rectangle sizes, vacancy caps and homogeneous-family covers.

Each digit independently chooses rows OR columns OR boxes, whichever has the
fewest occupied houses. Ties are all inspected for eliminations. Native common
peer intersections (beyond the chosen house itself) are allowed, as in source.
    """
    stats=Counter();findings=[];allowed=set(range(81)) if scope is None else set(scope)
    for nr,nc in [(3,3),(3,4),(4,3),(4,4),(4,5),(5,4)]:
        for rows in itertools.combinations(range(9),nr):
            if len({r//3 for r in rows})==1:continue
            for columns in itertools.combinations(range(9),nc):
                if len({c//3 for c in columns})==1:continue
                stats['rectangles_generated']+=1
                rectangle={9*r+c for r in rows for c in columns}
                cells={c for c in rectangle if domains[c]&(domains[c]-1)}
                if not cells<=allowed:continue
                stats['rectangles_in_scope']+=1
                if (len(rectangle)<12 and len(rectangle)-len(cells)>1) or len(rectangle)-len(cells)>2:continue
                stats['vacancy_eligible']+=1
                positions=[[c for c in cells if domains[c]>>(d-1)&1] for d in range(1,10)]
                housegroups=[];capacity=0
                for ds in positions:
                    families=[[i for i,h in enumerate(HOUSES[first:first+9],first) if any(c in h for c in ds)] for first in (0,9,18)]
                    minimum=min(map(len,families));capacity+=minimum
                    housegroups.append([family for family in families if len(family)==minimum])
                if capacity!=len(cells):continue
                stats['rank_zero_configurations']+=1;eliminations=[]
                for d,(ds,families) in enumerate(zip(positions,housegroups),1):
                    targets=set()
                    for family in families:
                        for h in family:
                            inside=set(ds)&HOUSES[h]
                            targets|={c for c in range(81) if domains[c]>>(d-1)&1 and all(peer(c,v) for v in inside)}
                    for c in sorted(targets):
                        if domains[c]&(domains[c]-1):eliminations.append([c,d])
                if eliminations:
                    stats['configurations_with_eliminations']+=1
                    findings.append({'cells':sorted(cells),'rectangle_rows':rows,'rectangle_columns':columns,
                     'cover_choices':housegroups,'eliminations':eliminations})
    for name in ('rectangles_generated','rectangles_in_scope','vacancy_eligible','rank_zero_configurations','configurations_with_eliminations'):stats[name]+=0
    return {'counts':dict(stats),'findings':findings}


@lru_cache(None)
def independent_capacity(cells):
    """Exact native packing number, upper bounded by tiny selected cross-lines."""
    cells=tuple(cells)
    @lru_cache(None)
    def walk(available):
        if not available:return 0
        v=available[0];rest=available[1:]
        return max(walk(rest),1+walk(tuple(w for w in rest if not peer(v,w))))
    return walk(cells)


def junior(domains,scope=None,relaxed=False):
    """Two-base, three cross-lines (one in each stack), two nonseeing targets.

Base union has three/four digits; companions contain none. S-cell capacity is
checked exactly (not guessed from a drawing), with solved values included.
All relevant base/targets/S support must lie in scope if supplied.
    """
    stats=Counter();patterns=[];rejections=[];allowed=set(range(81)) if scope is None else set(scope)
    for transpose in (False,True):
        def cell(r,c):return 9*c+r if transpose else 9*r+c
        for band in range(3):
            bandcells={cell(r,c) for r in range(3*band,3*band+3) for c in range(9)}
            for row in range(3*band,3*band+3):
                for stack in range(3):
                    for bp in itertools.combinations(range(3*stack,3*stack+3),2):
                        bases=[cell(row,c) for c in bp]
                        if any(not domains[c]&(domains[c]-1) for c in bases):continue
                        mask=domains[bases[0]]|domains[bases[1]]
                        if count(mask)<2 or not relaxed and count(mask) not in (3,4):continue
                        stats['eligible_base_pairs']+=1
                        remaining=next(c for c in range(3*stack,3*stack+3) if c not in bp)
                        other=[s for s in range(3) if s!=stack]
                        for crosschoices in itertools.product(range(3),repeat=2):
                            columns=sorted([remaining]+[3*s+k for s,k in zip(other,crosschoices)])
                            crossing={cell(r,c) for r in range(9) for c in columns}
                            inner=crossing&bandcells;s_cells=sorted(crossing-bandcells)
                            targets=sorted(c for c in inner if domains[c]&mask and not all(peer(c,b) for b in bases))
                            stats['crossline_geometries']+=1
                            if len(targets)!=2 or not relaxed and peer(*targets):continue
                            if any(not domains[c]&(domains[c]-1) for c in targets):continue
                            stats['two_live_targets' if relaxed else 'two_live_nonseeing_targets']+=1
                            relevant=set(bases+targets+[c for c in s_cells if domains[c]&mask])
                            if not relevant<=allowed:continue
                            stats['geometries_in_scope']+=1
                            caps={d:independent_capacity(tuple(c for c in s_cells if domains[c]>>(d-1)&1)) for d in digits(mask)}
                            if any(k>2 for k in caps.values()):
                                d=next(d for d,k in caps.items() if k>2)
                                witness=next(list(cs) for cs in itertools.combinations([c for c in s_cells if domains[c]>>(d-1)&1],3) if all(not peer(a,b) for a,b in itertools.combinations(cs,2)))
                                rejections.append({'bases':bases,'targets':targets,'base_digits':digits(mask),'s_cells':s_cells,'violating_digit':d,'three_native_nonpeer_cells':witness})
                                continue
                            stats['capacity_valid_patterns']+=1
                            eliminations=[[c,d] for c in targets for d in digits(domains[c]&~mask)]
                            patterns.append({'orientation':'column-base' if transpose else 'row-base','bases':bases,'targets':targets,
                             'base_digits':digits(mask),'crossline_indices':columns,'s_cells':s_cells,
                             'digit_capacities':caps,'eliminations':eliminations,
                             'digit_s_positions':{d:[c for c in s_cells if domains[c]>>(d-1)&1] for d in digits(mask)}})
    for name in ('eligible_base_pairs','crossline_geometries','two_live_targets' if relaxed else 'two_live_nonseeing_targets','geometries_in_scope','capacity_valid_patterns'):stats[name]+=0
    return {'counts':dict(stats),'patterns':patterns,'capacity_rejections':rejections}


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--fixture',type=Path,default=HERE/'different_lead.json');p.add_argument('--output',type=Path,default=HERE/'exocet_msls_results.json');a=p.parse_args()
    case=json.loads(a.fixture.read_text());records={}
    for name,domains in [('unforced',case['base_domains']),('closed_opposite',case['assumed_domains'])]:
        for scoped,scope in [('full_board',None),('selected20',case['graph_cells'])]:
            records[name+'_'+scoped]={'basic_junior_exocet':junior(domains,scope),'relaxed_three_crossline_capacity':junior(domains,scope,True),'rectangular_msls':msls(domains,scope)}
    result={'schema_version':1,'fixture_sha256':hashlib.sha256(a.fixture.read_bytes()).hexdigest(),
      'target_assumption':case['target'],'selected_cells':case['graph_cells'],'records':records,
      'limits':'Exactly the variants documented by this script. No Senior/Weak/Double/Complex Exocet or arbitrary derived-set/conditional MSLS exclusion.'}
    a.output.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:{t:v['counts'] for t,v in q.items()} for k,q in records.items()}))

if __name__=='__main__':main()
