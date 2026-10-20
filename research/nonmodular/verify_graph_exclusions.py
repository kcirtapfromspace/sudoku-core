"""Finite native-geometry exclusions; historical novelty is unverified.
Counts refer to diamond-pair configurations at one representative cell.
Cell transitivity covers every shared-tip position.
"""
if not __debug__:
    raise SystemExit("Run without -O: this research checker uses assertions.")


from itertools import combinations
N=81
P=[set() for _ in range(N)]
def peer(a,b):
 r,c=divmod(a,9);s,d=divmod(b,9)
 return a!=b and (r==s or c==d or (r//3,c//3)==(s//3,d//3))
for a in range(N):P[a]={b for b in range(N) if peer(a,b)}
S=0
diamonds=[]
for a,b in combinations(sorted(P[S]),2):
 if b not in P[a]:continue
 for t in (P[a]&P[b])-P[S]-{S}:
  diamonds.append((a,b,t))
print('Diamonds with common tip r1c1:',len(diamonds))
count=0;hits=[]
for d1,d2 in combinations(diamonds,2):
 if set(d1)&set(d2) or d2[2] not in P[d1[2]]:continue
 V=[S]+list(d1)+list(d2)
 if any(all(peer(a,b) for a,b in combinations(Q,2)) for Q in combinations(V,4)):continue
 expected={tuple(sorted(e)) for e in [(S,d1[0]),(S,d1[1]),(d1[0],d1[1]),(d1[0],d1[2]),(d1[1],d1[2]),(S,d2[0]),(S,d2[1]),(d2[0],d2[1]),(d2[0],d2[2]),(d2[1],d2[2]),(d1[2],d2[2])]}
 actual={tuple(sorted((a,b))) for a,b in combinations(V,2) if peer(a,b)}
 count+=1
 hits.append((V,actual==expected,len(actual)))
print('Native K4-free spindle supergraph embeddings:',count)
if hits:
 V,induced,edges=hits[0]
 print('First:',[(v//9+1,v%9+1) for v in V],'induced',induced,'edges',edges)
else: print('No such embedding at any common tip, by cell transitivity.')
# Any K4-free odd wheel implies the hub neighborhood contains an induced odd cycle.
# Search induced cycles of length 5,7,... up to possible local-graph max hole size.
hood=P[S];found=[]
for size in range(5, min(len(hood),9)+1,2):
 def dfs(path):
  if len(path)==size:
   if path[0] in P[path[-1]]:found.append(path);return True
   return False
  for b in sorted(hood&P[path[-1]]):
   if b<=path[0] or b in path:continue
   if any(b in P[a] for a in path[:-1] if not(len(path)==size-1 and a==path[0])):continue
   if dfs(path+[b]):return True
  return False
 for start in sorted(hood):
  if dfs([start]):break
 if found:break
print('Induced odd wheels with rim size 5,7,9:',found)
# All omitted extra-edge cases contain K4 already inside a diamond. For the
# induced-diamond cases above, count exactly how every spindle support collapses.
raw_count=0;k4_hist={}
for d1,d2 in combinations(diamonds,2):
 if set(d1)&set(d2) or d2[2] not in P[d1[2]]:continue
 V=[S]+list(d1)+list(d2)
 number=sum(all(peer(a,b) for a,b in combinations(Q,2)) for Q in combinations(V,4))
 raw_count+=1;k4_hist[number]=k4_hist.get(number,0)+1
assert raw_count==2376 and k4_hist=={5:648,1:1296,2:432}
print('All spindle embeddings:',raw_count,'K4-count histogram:',k4_hist)
# A perfect elimination ordering proves the entire native neighborhood chordal,
# not merely the bounded odd-cycle checks above. First remove row/column peers
# outside the center box, then the remaining eight box peers.
order=[v for v in sorted(hood) if v//9==0 and v%9>=3]
order += [v for v in sorted(hood) if v%9==0 and v//9>=3]
order += [v for v in sorted(hood) if v not in order]
assert len(order)==20 and set(order)==hood
for i,v in enumerate(order):
 later=P[v]&set(order[i+1:])
 assert all(peer(a,b) for a,b in combinations(later,2))
print('Exact perfect elimination ordering verifies native neighborhoods are chordal.')
print('Therefore a K4-free native wheel of ANY rim length is impossible:')
print('its rim support would be triangle-free and chordal, hence a forest.')
