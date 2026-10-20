"""Candidate-preserving, explicitly bounded OFFLINE Sudoku discovery baseline.

Binary reachability and guardian unit propagation are only a subset of
CSP-Rules ORk-/forcing-whips. Every enabled profile and every cutoff is exposed.
"""
import argparse
from copy import deepcopy
from dataclasses import dataclass,asdict,replace
import json
from pathlib import Path
import time

from strong_baseline_relations import (Budget,Cutoff,state,discover,covered,digest,
                                      valid_domains,unit_step,binary_step,clause_status)
from deeper_forcing_search import propagate
from deeper_forcing_verify import blockers
from exocet_msls_scan import junior,msls


@dataclass
class Options:
    context_id:str=None
    seconds:float=60
    max_rounds:int=1000
    catalogue:bool=True
    catalogue_seconds:float=10
    catalogue_nodes:int=100000
    max_clauses:int=20000
    unit_binary:bool=True
    junior:bool=True
    relaxed_junior:bool=True
    msls:bool=True
    forcing_depth:int=2
    forcing_nodes:int=20000
    forcing_fresh_catalogue:bool=True
    forcing_inner_je_msls:bool=False
    forcing_guardian_width:int=6
    forcing_targets:tuple=()
    history_snapshots:tuple=()


def apply_step(domains,step):
    cell,digit=step['cell'],step['digit'];bit=1<<(digit-1)
    if not domains[cell]&bit:raise ValueError('Inference endpoint absent')
    if step.get('value',False):domains[cell]=bit
    else:domains[cell]&=~bit


def terminal_node(closed):
    result=deepcopy(closed)
    if result['kind']=='three-house-contradiction':
        result['cover']=blockers(result['domains'],result['center'])
    return result


def nonforcing(initial,relations,options,budget,branch=False):
    domains=initial.copy();trace=[];phases=[]
    try:
        for iteration in range(options.max_rounds):
            budget.tick(0);before=domains.copy();round_phases=[]
            do_catalogue=options.catalogue and (not branch or options.forcing_fresh_catalogue)
            if do_catalogue:
                cb=Budget(min(options.catalogue_seconds,max(0,budget.deadline-time.monotonic())),options.catalogue_nodes)
                old=len(relations['clauses']);scan=discover(domains,relations,cb,options.max_clauses)
                new=deepcopy(relations['clauses'][old:])
                trace.append({'kind':'catalogue-scan','origin_domains':domains.copy(),'new_clauses':new,'scan':scan})
                round_phases.append({'phase':'catalogue','status':scan['status']})
                if scan['status']=='cutoff':return domains,trace,relations,'cutoff',phases+round_phases
            elif options.catalogue and branch:
                round_phases.append({'phase':'catalogue','status':'omitted-in-declared-inner-profile'})
            closed=propagate(tuple(domains))
            trace.append({'kind':'local','certificate':terminal_node(closed)if closed['kind']!='fixedpoint'else deepcopy(closed)})
            domains=closed['domains'].copy();round_phases.append({'phase':'Hall81','status':'complete'})
            if closed['kind']!='fixedpoint':return domains,trace,relations,'refuted',phases+round_phases
            budget.tick(0)
            if options.unit_binary:
                found=unit_step(domains,relations)
                if found and found['kind']=='guardian-contradiction':
                    trace.append(found);return domains,trace,relations,'refuted',phases+round_phases
                if found is None:found=binary_step(domains,relations,budget)
                if found:
                    trace.append(found);apply_step(domains,found)
                    phases+=round_phases+[{'phase':'unit/binary','status':'progress'}]
                    continue
                round_phases.append({'phase':'unit/binary','status':'complete'})
            do_extended=not branch or options.forcing_inner_je_msls
            changed=False
            if do_extended:
                scans=[]
                if options.junior:scans.append(('junior',lambda:junior(domains)))
                if options.relaxed_junior:scans.append(('relaxed-junior',lambda:junior(domains,relaxed=True)))
                if options.msls:scans.append(('rectangle-msls',lambda:msls(domains)))
                for name,scanfun in scans:
                    budget.tick(0);result=scanfun();budget.tick(0)
                    patterns=result.get('patterns',result.get('findings',[]))
                    finding=next((p for p in patterns if p['eliminations']),None)
                    round_phases.append({'phase':name,'status':'complete','counts':result['counts']})
                    if finding:
                        cell,digit=finding['eliminations'][0]
                        step={'kind':'pattern-elimination','rule':name,'pattern':deepcopy(finding),'cell':cell,'digit':digit,'value':False}
                        trace.append(step);apply_step(domains,step);changed=True;break
            elif options.junior or options.relaxed_junior or options.msls:
                round_phases.append({'phase':'JE/MSLS','status':'omitted-in-declared-inner-profile'})
            phases+=round_phases
            if changed or domains!=before:continue
            return domains,trace,relations,'complete',phases
        raise Cutoff('Monotone nonforcing round budget exhausted')
    except Cutoff as exc:
        return domains,trace,relations,'cutoff',phases+[{'phase':'nonforcing-incomplete','status':'cutoff','reason':str(exc)}]



class Forcing:
    def __init__(self,options,budget):
        self.options=options;self.budget=budget;self.calls=0
    def refute(self,initial,relations,depth):
        self.calls+=1;self.budget.tick()
        if self.calls>self.options.forcing_nodes:raise Cutoff('Forcing node budget exhausted')
        domains,trace,local,status,phases=nonforcing(initial,deepcopy(relations),self.options,self.budget,True)
        if status=='cutoff':raise Cutoff('Assumption phase incomplete')
        node={'trace':trace,'domains':domains,'phase_results':phases}
        if status=='refuted':return {**node,'kind':'refuted'}
        if depth==0:return None
        choices=[]
        if self.options.unit_binary:
            for clause in local['clauses']:
                st,left=clause_status(clause,domains)
                if st=='open'and 2<=len(left)<=self.options.forcing_guardian_width:
                    choices.append(('guardian',clause['id'],left))
        for c in sorted((c for c,m in enumerate(domains)if m&(m-1)),key=lambda c:(bin(domains[c]).count('1'),c)):
            choices.append(('cell',c,[[c,d]for d in range(1,10)if domains[c]&(1<<(d-1))]))
        for kind,pivot,literals in choices:
            branches=[]
            for cell,digit in literals:
                assumed=domains.copy();assumed[cell]=1<<(digit-1)
                child=self.refute(assumed,local,depth-1)
                if child is None:break
                branches.append({'literal':[cell,digit],'proof':child})
            else:return {**node,'kind':'split','split_kind':kind,'pivot':pivot,'branches':branches}
        return None


def close_state(domains,inherited_relations=None,options=None):
    options=options or Options();options=Options(**options)if type(options)is dict else options
    valid_domains(domains);initial=domains.copy();working=domains.copy()
    relations=state(working,inherited_relations,options.context_id)
    input_relations=deepcopy(relations)
    budget=Budget(options.seconds,max(1000000,options.forcing_nodes*100));trace=[];phases=[];status='cutoff';error=None
    try:
        # History scope is exactly supplied snapshots plus this starting state.
        # No claim that an unavailable givens-to-state history was scanned.
        for snapshot in options.history_snapshots:
            valid_domains(list(snapshot))
            if any(m&~old for m,old in zip(working,snapshot)):raise ValueError('History snapshot not an ancestor')
            cb=Budget(min(options.catalogue_seconds,max(0,budget.deadline-time.monotonic())),options.catalogue_nodes)
            old=len(relations['clauses']);scan=discover(list(snapshot),relations,cb,options.max_clauses)
            trace.append({'kind':'history-catalogue-scan','origin_domains':list(snapshot),
                          'new_clauses':deepcopy(relations['clauses'][old:]),'scan':scan})
            phases.append({'phase':'history-catalogue','status':scan['status']})
            if scan['status']=='cutoff':raise Cutoff('History catalogue scope incomplete')
        for _ in range(options.max_rounds):
            working,part,relations,status,phase=nonforcing(working,relations,options,budget)
            trace+=part;phases+=phase
            if status in ('cutoff','refuted'):break
            if all(m&(m-1)==0 for m in working):status='complete';break
            if options.forcing_depth<0:status='fixedpoint';break
            forcing=Forcing(options,budget);progress=False
            targets=list(options.forcing_targets)if options.forcing_targets else [(c,d)for c,m in enumerate(working)if m&(m-1)for d in range(1,10)if m&(1<<(d-1))]
            for cell,digit in targets:
                if not working[cell]&(1<<(digit-1))or working[cell]&(working[cell]-1)==0:continue
                assumed=working.copy();assumed[cell]=1<<(digit-1)
                proof=forcing.refute(assumed,relations,options.forcing_depth)
                if proof:
                    step={'kind':'forcing-elimination','cell':cell,'digit':digit,'value':False,'proof':proof}
                    trace.append(step);apply_step(working,step);progress=True;break
            phases.append({'phase':'forcing','status':'progress'if progress else'complete',
                           'calls':forcing.calls,'candidate_scope':'explicit focus list'if options.forcing_targets else'all live candidates'})
            if not progress:status='fixedpoint';break
        else:raise Cutoff('Top-level round budget exhausted')
    except Cutoff as exc:status='cutoff';error=str(exc)
    return {'schema_version':1,'input_domains':initial,'domains':working,'trace':trace,'relations':relations,'input_relations':input_relations,
            'status':status,'error':error,'phase_results':phases,'options':asdict(options),
            'history_scope':{'starting_state_hash':digest(initial),'provided_ancestor_snapshots':len(options.history_snapshots),
                             'givens_history_claimed':False},
            'forcing_profile':{'additional_split_depth':options.forcing_depth,'fresh_catalogue':options.forcing_fresh_catalogue,
                               'inner_JE_MSLS':options.forcing_inner_je_msls,'guardian_split_width':options.forcing_guardian_width},
            'scope':'Offline configured baseline only. Unit/binary guardian reasoning is a subset of ORk-whips/forcing-whips, not their full implementation. No fixedpoint claim on cutoff.',
            'seconds':time.monotonic()-budget.started}


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('fixture',type=Path);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--opposite',action='store_true');p.add_argument('--seconds',type=float,default=60)
    p.add_argument('--no-forcing',action='store_true');p.add_argument('--inner-je-msls',action='store_true');p.add_argument('--history',action='store_true')
    a=p.parse_args();case=json.loads(a.fixture.read_text());domains=case['assumed_domains']if a.opposite else case.get('base_domains',case['state']['domains'])
    options=Options(seconds=a.seconds,forcing_depth=-1 if a.no_forcing else 2,forcing_inner_je_msls=a.inner_je_msls,
                    history_snapshots=(case['initial_domains'],)if a.history else())
    result=close_state(domains,options=options);a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps(result,separators=(',',':'))+'\n');print(json.dumps({k:result[k]for k in ['status','error','seconds','history_scope','forcing_profile']}))


if __name__=='__main__':main()
