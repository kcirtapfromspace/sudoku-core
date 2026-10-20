"""Bounded exploratory clue augmentation, with reproducible negative results.

SciPy LP is a floating discovery screen only. These outputs are not certificates
or committed fixtures. Positive acceptance would require real_verify plus an
exact linear_witness artifact. The baseline is the frozen pre-global engine.
"""
import argparse,json,sys,random,time,subprocess,tempfile
from pathlib import Path
import numpy as np
from scipy.optimize import linprog
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE))
from real_verify import close,exact_completions
from real_existing_check import build_engine
args=argparse.ArgumentParser(description=__doc__)
args.add_argument('--seeds',type=int,default=30)
args.add_argument('--binary',type=Path,help='optional already-built frozen baseline probe')
args.add_argument('--output-dir',type=Path,default=Path('/tmp/sudoku-global-q4-augmentation'))
args=args.parse_args()
args.output_dir.mkdir(parents=True,exist_ok=True)
temporary=None
if args.binary is None:
 temporary=tempfile.TemporaryDirectory(prefix='global-q4-augmentation-engine-')
 binary=build_engine(temporary.name)
else:binary=args.binary.resolve()
from real_four_digit_mine import initial_domains
r=json.loads((HERE/'real_four_digit_candidate.json').read_text());sol=list(map(int,r['completion']));core=r['graph_cells'];guardian=r['target']['cell'];extra=r['target']['digit'];palette=r['palette']
A=np.zeros((324,729),dtype=np.int8)
for c in range(81):
 rr,cc=divmod(c,9);b=rr//3*3+cc//3
 for d in range(9):A[[c,81+rr*9+d,162+cc*9+d,243+b*9+d],c*9+d]=1

def lp(ds):
 vs=[c*9+d for c,m in enumerate(ds) for d in range(9) if m&(1<<d)]
 return linprog(np.zeros(len(vs)),A_eq=A[:,vs],b_eq=np.ones(324),bounds=(0,None),method='highs-ds').success

def engine_close(puzzle,base,unique):
 alltrace=[]
 for i in range(50):
  state={'version':1,'domains':base,'values':[m.bit_length() if m&(m-1)==0 else 0 for m in base]}
  e=json.loads(subprocess.check_output([binary],input=json.dumps({'id':'augmentation','puzzle':puzzle,'state':state,'unique':unique,'full':False})+'\n',text=True))
  if e['state']['domains'][guardian]==1:return None,None
  if not e['trace']:return base,alltrace
  alltrace.append({'kind':'existing-engine','steps':e['trace'],'state_after':e['state']})
  base,trace=close(e['state']['domains']);alltrace+=trace
  if base[guardian]==1:return None,None
 raise Exception('engine local closure limit')

seedgivens=[int(c) if c!='.' else 0 for c in r['puzzle']]
start=time.time();attempts=[];best=None;totaltested=0
for seed in range(args.seeds):
 rng=random.Random(21410+seed);g=seedgivens.copy();available=[c for c in range(81) if c not in core and not g[c]];rng.shuffle(available);available.sort(key=lambda c:sol[c] not in palette)
 accepted=[];tested=0
 for c in available:
  tested+=1;totaltested+=1;n=g.copy();n[c]=sol[c];ds=initial_domains(n);a=ds.copy();a[guardian]&=~(1<<(extra-1))
  if not lp(a):continue
  puzzle=''.join(str(d) if d else '.' for d in n)
  try:
   base,bt=close(ds);base,et=engine_close(puzzle,base,False)
   if base is None:continue
   bt+=et;a=base.copy();a[guardian]&=~(1<<(extra-1));branch,at=close(a)
  except ValueError:continue
  if not lp(branch):continue
  s=exact_completions(puzzle)
  if len(s)==1:
   base,et=engine_close(puzzle,base,True)
   if base is None:continue
   bt+=et;a=base.copy();a[guardian]&=~(1<<(extra-1))
   try:branch,at=close(a)
   except ValueError:continue
   if not lp(branch):continue
  g=n;accepted.append(c)
  out={'schema_version':1,'name':'generated-four-digit-hard-unique' if len(s)==1 else 'generated-four-digit-hard-augmentation','puzzle':puzzle,'initial_domains':ds,'base_domains':base,'base_trace':bt,'assumed_domains':branch,'assumption_trace':at,'completion':r['completion'],'solution':r['completion'],'graph_cells':core,'palette':palette,'target':r['target'],'solutions_capped_at_2':len(s),'seed':21410+seed,'additional_givens':accepted,'seconds':time.time()-start,'trials':totaltested}
  for key,dom in [('state',base),('assumed_state',branch)]:out[key]={'version':1,'domains':dom,'values':[m.bit_length() if m&(m-1)==0 else 0 for m in dom]}
  if best is None or sum(bool(d) for d in g)>sum(ch!='.' for ch in best['puzzle']):
   best=out;(args.output_dir/'best.json').write_text(json.dumps(out,indent=2)+'\n')
  print('seed',seed,'givens',sum(bool(d) for d in g),'added',c,'solutions',len(s),'seconds',time.time()-start,flush=True)
  if len(s)==1:
   (args.output_dir/'unique.json').write_text(json.dumps(out,indent=2)+'\n');print('UNIQUE HARD',puzzle,flush=True);raise SystemExit
 attempts.append({'seed':seed,'tested':tested,'accepted':accepted,'givens':sum(bool(d) for d in g)})
 print('finished seed',seed,'givens',sum(bool(d) for d in g),flush=True)
(args.output_dir/'summary.json').write_text(json.dumps({'attempts':attempts,'seconds':time.time()-start,'best':best},indent=2)+'\n')

if temporary is not None:temporary.cleanup()
