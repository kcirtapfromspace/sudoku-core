#!/usr/bin/env python3
"""Direct exact theorem checks for the declared Junior/rectangular-MSLS steps.

No pattern searcher is imported. Acceptance proves the supplied instance; it
does not claim an exhaustive search for other Exocet or MSLS variants.
"""
import itertools
import argparse
import copy
import json
from pathlib import Path


def need(ok,message):
    if not ok:raise ValueError(message)


def peer(a,b):return a!=b and (a//9==b//9 or a%9==b%9 or (a//27,a%9//3)==(b//27,b%9//3))
def digits(mask):return [d for d in range(1,10)if mask&(1<<(d-1))]
def live(mask):return mask&(mask-1)!=0


def houses():
    return [set(9*r+c for c in range(9))for r in range(9)]+[set(9*r+c for r in range(9))for c in range(9)]+[
        set(9*(3*b+r)+3*s+c for r in range(3)for c in range(3))for b in range(3)for s in range(3)]


def junior(domains,pattern,relaxed):
    orientation=pattern['orientation'];need(orientation in ('row-base','column-base'),'base orientation')
    transpose=orientation=='column-base'
    def canonical(c):
        r,k=divmod(c,9);return (k,r)if transpose else(r,k)
    def cell(r,c):return 9*c+r if transpose else 9*r+c
    bases=pattern['bases'];need(len(bases)==2 and len(set(bases))==2,'two distinct bases')
    need(all(type(c)is int and 0<=c<81 and live(domains[c])for c in bases),'live native bases')
    (r,a),(s,b)=map(canonical,bases)
    need(r==s and a//3==b//3,'same-row same-box base geometry')
    mask=domains[bases[0]]|domains[bases[1]]
    need(pattern['base_digits']==digits(mask),'actual base union')
    need(len(digits(mask))>=2 and (relaxed or len(digits(mask))in (3,4)),'declared base width')
    cross=pattern['crossline_indices']
    need(cross==sorted(set(cross))and len(cross)==3 and all(type(c)is int and 0<=c<9 for c in cross),'three crosslines')
    need({c//3 for c in cross}=={0,1,2},'one crossline per stack')
    remainder=next(c for c in range(3*(a//3),3*(a//3)+3)if c not in (a,b))
    need(remainder in cross,'base companion crossline')
    band={cell(row,col)for row in range(3*(r//3),3*(r//3)+3)for col in range(9)}
    crossing={cell(row,col)for row in range(9)for col in cross}
    expected_s=sorted(crossing-band);need(pattern['s_cells']==expected_s,'complete S cells')
    targets=sorted(c for c in crossing&band if domains[c]&mask and not all(peer(c,b)for b in bases))
    need(pattern['targets']==targets and len(targets)==2 and all(live(domains[c])for c in targets),'exact two live targets')
    need(relaxed or not peer(*targets),'ordinary Junior target geometry')
    positions={str(d):[c for c in expected_s if domains[c]&(1<<(d-1))]for d in digits(mask)}
    need({str(k):v for k,v in pattern['digit_s_positions'].items()}==positions,'actual S support lists')
    caps={}
    for key,cells in positions.items():
        need(not any(all(not peer(a,b)for a,b in itertools.combinations(cs,2))for cs in itertools.combinations(cells,3)),'S support permits three placements')
        caps[key]=0 if not cells else 2 if any(not peer(a,b)for a,b in itertools.combinations(cells,2))else 1
    need({str(k):v for k,v in pattern['digit_capacities'].items()}==caps,'exact S capacity')
    # Each actual base digit needs3 occurrences in the3 crosslines, at most2
    # in S; cells seeing both bases cannot host it. The two targets therefore
    # contain the two distinct actual base digits and no outside digit.
    return [[c,d]for c in targets for d in digits(domains[c]&~mask)]


def msls(domains,pattern):
    rows=list(pattern['rectangle_rows']);columns=list(pattern['rectangle_columns'])
    need(rows==sorted(set(rows))and columns==sorted(set(columns)),'rectangle coordinate sets')
    need(all(type(x)is int and 0<=x<9 for x in rows+columns),'native rectangle range')
    need((len(rows),len(columns))in [(3,3),(3,4),(4,3),(4,4),(4,5),(5,4)],'declared rectangle sizes')
    need(len({r//3 for r in rows})>1 and len({c//3 for c in columns})>1,'declared rectangle scope')
    rectangle={9*r+c for r in rows for c in columns};cells={c for c in rectangle if live(domains[c])}
    need(pattern['cells']==sorted(cells),'complete live rectangle support')
    vacancies=len(rectangle)-len(cells)
    need(vacancies<=2 and (len(rectangle)>=12 or vacancies<=1),'declared vacancy scope')
    hs=houses();cover_choices=[];positions=[];capacity=0
    for d in range(1,10):
        ds={c for c in cells if domains[c]&(1<<(d-1))};positions.append(ds)
        families=[[i for i in range(start,start+9)if hs[i]&ds]for start in (0,9,18)]
        minimum=min(map(len,families));capacity+=minimum
        cover_choices.append([f for f in families if len(f)==minimum])
    need(pattern['cover_choices']==cover_choices and capacity==len(cells),'exact rank-zero cover equality')
    # Homogeneous covers partition each digit's selected positions. Sum of
    # capacities equals |S|, hence every minimum cover house is saturated.
    eliminations=[]
    for d,(ds,families)in enumerate(zip(positions,cover_choices),1):
        targets=set()
        for family in families:
            for h in family:
                inside=hs[h]&ds;need(inside,'nonempty saturated group')
                targets|={c for c in range(81)if domains[c]&(1<<(d-1))and all(peer(c,v)for v in inside)}
        eliminations.extend([c,d]for c in sorted(targets)if live(domains[c]))
    return eliminations


def pattern_step(domains,event):
    need(__debug__,'Run without -O')
    need(len(domains)==81 and all(type(m)is int and 0<m<=511 for m in domains),'native domains')
    need(event['kind']=='pattern-elimination'and event['value']is False,'pattern elimination polarity')
    need(type(event['cell'])is int and 0<=event['cell']<81 and type(event['digit'])is int and 1<=event['digit']<=9,'native elimination endpoint')
    rule=event['rule'];pattern=event['pattern']
    if rule in ('junior','relaxed-junior'):eliminations=junior(domains,pattern,rule=='relaxed-junior')
    elif rule=='rectangle-msls':eliminations=msls(domains,pattern)
    else:raise ValueError('Outside declared pattern scope')
    need(pattern['eliminations']==eliminations,'exact claimed pattern eliminations')
    need([event['cell'],event['digit']]in eliminations,'requested elimination established')
    return True


def self_test():
    records=json.loads((Path(__file__).resolve().parent/'joint_audit_pattern_controls.json').read_text())['records']
    need(len(records)==2,'both nonvacuous pattern controls')
    rejected=0
    for record in records:
        domains=record['domains'];event=record['event'];completion=list(map(int,record['completion']))
        need(len(completion)==81 and all({completion[c]for c in h}==set(range(1,10))for h in houses()),'respecting completed grid')
        need(all(domains[c]&(1<<(d-1))for c,d in enumerate(completion)),'completion respects all lists')
        pattern_step(domains,event)
        need(all(completion[c]!=d for c,d in event['pattern']['eliminations']),'all actual eliminations respect completion')
        changes=[lambda e:e.__setitem__('value',True),lambda e:e.__setitem__('digit',0),
                 lambda e:e['pattern']['eliminations'].clear()]
        if event['rule']=='junior':
            changes.extend([lambda e:e['pattern'].__setitem__('base_digits',[1,2]),
                            lambda e:e['pattern']['s_cells'].pop(),
                            lambda e:e['pattern']['digit_capacities'].__setitem__('1',3)])
        else:
            changes.extend([lambda e:e['pattern']['cells'].pop(),
                            lambda e:e['pattern']['cover_choices'][0].clear(),
                            lambda e:e['pattern']['rectangle_columns'].pop()])
        for change in changes:
            bad=copy.deepcopy(event);change(bad)
            try:pattern_step(domains,bad)
            except (ValueError,KeyError,TypeError,IndexError):rejected+=1
            else:raise ValueError('Invalid pattern control accepted')
    return {'accepted_nonvacuous_patterns':len(records),'corruptions_rejected':rejected,
            'scope':'Direct basic-Junior and declared rectangular-MSLS theorem checks only.'}


def main():
    need(__debug__,'Run without -O');p=argparse.ArgumentParser(description=__doc__);p.add_argument('--self-test',action='store_true');args=p.parse_args()
    need(args.self_test,'Use --self-test or import pattern_step');print(json.dumps(self_test(),sort_keys=True))


if __name__=='__main__':main()
