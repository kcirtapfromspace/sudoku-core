"""Exact sanity benchmark for rational PSD Sudoku certificates.
The pattern is the KNOWN bivalue oddagon, not a novelty claim.
No repository files are modified; all arithmetic checks use fractions.Fraction.
"""
if not __debug__:
    raise SystemExit("Run without -O: this research checker uses assertions.")


from fractions import Fraction as F
from itertools import product

CYCLE = [(0, 0), (0, 3), (1, 4), (3, 4), (3, 0)]
COMPLETION = '324156789897324156156789324249618537631572498578493261963241875712865943485937612'
# A y-value 1/2 means x_(cell,1)=x_(cell,2)=1/2; y=0 means
# x_(cell,d)=1/7 for all d=3,...,9. There are two y=1/2 cells
# in every row, column, and box, including every cell in CYCLE.
HALVES = {
 (0,0),(0,3), (1,4),(1,6), (2,2),(2,7),
 (3,0),(3,4), (4,3),(4,8), (5,1),(5,8),
 (6,2),(6,5), (7,1),(7,7), (8,5),(8,6),
}

def peers(a, b):
 r,c=a; s,d=b
 return r==s or c==d or (r//3,c//3)==(s//3,d//3)

# Establish exactly a chordless 5-cycle under native Sudoku peer relations.
for i in range(5):
 for j in range(i+1,5):
  assert peers(CYCLE[i], CYCLE[j]) == ((i-j)%5 in (1,4))

units = ([[9*r+c for c in range(9)] for r in range(9)] +
         [[9*r+c for r in range(9)] for c in range(9)] +
         [[9*r+c for r in range(3*br,3*br+3)
                    for c in range(3*bc,3*bc+3)] for br in range(3) for bc in range(3)])
# Exact satisfiability witness for guardian state: A allows 1,2,3; B-E allow 1,2.
board = list(map(int, COMPLETION))
assert all({board[i] for i in unit} == set(range(1,10)) for unit in units)
assert board[0] == 3
assert all(board[9*r+c] in (1,2) for r,c in CYCLE[1:])

# Full 729-variable fractional solution after additionally deleting guardian A3.
def x(r,c,d):
 y=F(1,2) if (r,c) in HALVES else F(0)
 return y if d in (1,2) else (1-2*y)/7
assert all(sum(x(r,c,d) for d in range(1,10)) == 1 for r in range(9) for c in range(9))
assert all(sum(x(i//9,i%9,d) for i in unit) == 1 for unit in units for d in range(1,10))
assert all(x(r,c,d)==0 for r,c in CYCLE for d in range(3,10))
print('Full 324-equation native LP relaxation: exactly feasible, denominators only 2 and 7.')

# M=8Q has diagonal 10, C5 edge coefficients 3, and nonedge coefficients -8.
M=[[F(10 if i==j else 3 if (i-j)%5 in (1,4) else -8) for j in range(5)] for i in range(5)]
L=[[F(i==j) for j in range(5)] for i in range(5)]
D=[F(0)]*5
for j in range(5):
 D[j]=M[j][j]-sum(L[j][k]**2*D[k] for k in range(j))
 for i in range(j+1,5):
  residual=M[i][j]-sum(L[i][k]*L[j][k]*D[k] for k in range(j))
  assert D[j] or not residual
  if D[j]:L[i][j]=residual/D[j]
assert all(d>=0 for d in D)
assert all(M[i][j]==sum(L[i][k]*D[k]*L[j][k] for k in range(5)) for i in range(5) for j in range(5))
print('Rational PSD factorization verified. LDL pivots:', ', '.join(map(str,D)))

# Independent exhaustive combinatorial check, separate from PSD factorization.
stable = [z for z in product((0,1), repeat=5) if all(not(z[i] and z[(i+1)%5]) for i in range(5))]
assert len(stable)==11 and max(map(sum,stable))==2
assert not any(all(z[i]!=z[(i+1)%5] for i in range(5)) for z in product((0,1),repeat=5))
print('All 32 binary patterns checked: max independent set 2; no 2-coloring.')
print('Known oddagon deduction: r1c1=3.')
print('Respecting completion:', COMPLETION)
