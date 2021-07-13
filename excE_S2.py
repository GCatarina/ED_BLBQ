# packages #
import sys
import numpy as np
import time
from quspin.basis import spin_basis_1d
from quspin.operators import hamiltonian
###############################################################################

# general functions #
# Hamiltonian: BLBQ, 1D
def H_BLBQ_1D(N,J,beta,BC,basis):
	#J
	SpSm = [[J/2,i,(i+1)%N] for i in range(N-1+BC)]
	SmSp = [[J/2,i,(i+1)%N] for i in range(N-1+BC)]
	SzSz = [[J,i,(i+1)%N] for i in range(N-1+BC)]

	#β
	SzSzSzSz = [[beta*J,i,i,(i+1)%N,(i+1)%N] for i in range(N-1+BC)]
	SzSpSzSm = [[beta*J/2,i,i,(i+1)%N,(i+1)%N] for i in range(N-1+BC)]
	SzSmSzSp = [[beta*J/2,i,i,(i+1)%N,(i+1)%N] for i in range(N-1+BC)]
	SpSzSmSz = [[beta*J/2,i,i,(i+1)%N,(i+1)%N] for i in range(N-1+BC)]
	SmSzSpSz = [[beta*J/2,i,i,(i+1)%N,(i+1)%N] for i in range(N-1+BC)]
	SpSpSmSm = [[beta*J/4,i,i,(i+1)%N,(i+1)%N] for i in range(N-1+BC)]
	SpSmSmSp = [[beta*J/4,i,i,(i+1)%N,(i+1)%N] for i in range(N-1+BC)]
	SmSpSpSm = [[beta*J/4,i,i,(i+1)%N,(i+1)%N] for i in range(N-1+BC)]
	SmSmSpSp = [[beta*J/4,i,i,(i+1)%N,(i+1)%N] for i in range(N-1+BC)]

	static = [
	    ["+-",SpSm],
	    ["-+",SmSp],
	    ["zz",SzSz],
		["zzzz", SzSzSzSz],
		["z+z-", SzSpSzSm],
		["z-z+", SzSmSzSp],
		["+z-z", SpSzSmSz],
		["-z+z", SmSzSpSz],
		["++--", SpSpSmSm],
		["+--+", SpSmSmSp],
		["-++-", SmSpSpSm],
		["--++", SmSmSpSp]
	]
	dynamic = []

	no_checks = dict(check_pcon=False,check_symm=False,check_herm=False)
	H = hamiltonian(static,dynamic,basis=basis,dtype=np.float64,**no_checks)

	return H

# Sz(i) operator
def Szi_op(i,basis):
	Sz = [[1,i]]

	static = [
		['z',Sz]
	]
	dynamic = []

	no_checks = dict(check_pcon=False,check_symm=False,check_herm=False)
	Szi = hamiltonian(static,dynamic,basis=basis,dtype=np.float64,**no_checks)

	return Szi

# Sz operator
def Sz_op(Nsites,basis):
	Sz = 0
	for i in range(Nsites):
		Sz += Szi_op(i,basis)

	return Sz

# S^2 operator
def S2_op(Nsites,basis):
	Sz = Sz_op(Nsites,basis)

	#S+S-
	SpSm = [[1,i,j] for i in range(Nsites) for j in range(Nsites)]

	static = [
		['+-',SpSm]
	]
	dynamic = []

	no_checks = dict(check_pcon=False,check_symm=False,check_herm=False)
	SpSm = hamiltonian(static,dynamic,basis=basis,dtype=np.float64,**no_checks)

	S2 = SpSm + np.dot(Sz,Sz) - Sz

	return S2
###############################################################################

# main #
start_time = time.time()

# read inputs
if len(sys.argv) != 9:
	print('Error: run code as #python_<code.py>_<s>_<N>_<J(meV)>_<beta>_<BC>'
		+ '_<Sz>_<DeltaE(meV)>_<nLanczos>\n')
	sys.exit(0)

# physical parameters
## s
try:
	s = float(sys.argv[1])
except:
	print('Error: insert integer or half-integer s > 0\n')
	sys.exit(0)
if (2*s)%1 != 0 or s <= 0:
	print('Error: insert integer or half-integer s > 0\n')
	sys.exit(0)
## N
try:
	N = int(sys.argv[2])
except:
	print('Error: insert integer N >= 2\n')
	sys.exit(0)
if N%1 != 0 or N <= 1:
	print('Error: insert integer N >= 2\n')
	sys.exit(0)
## J
try:
	J = float(sys.argv[3])
except:
	print('Error: insert real J != 0\n')
	sys.exit(0)
if J==0:
	print('Error: insert real J != 0\n')
	sys.exit(0)
## beta
try:
	beta = float(sys.argv[4])
except:
	print('Error: insert real beta\n')
	sys.exit(0)
## BC
try:
	BC = int(sys.argv[5])
except:
	print('Error: insert BC=0 (1) for open (periodic) boundary conditions\n')
	sys.exit(0)
if BC!=0 and BC!=1:
	print('Error: insert BC=0 (1) for open (periodic) boundary conditions\n')
	sys.exit(0)

# other parameters
## Sz
try:
	Sz = float(sys.argv[6])
except:
	print('Error: insert Sz in [-N*s, -N*s+1, ..., N*s]\n')
	sys.exit(0)
if Sz%1 != (N*s)%1 or abs(Sz) > N*s:
	print('Error: insert Sz in [-N*s, -N*s+1, ..., N*s]\n')
	sys.exit(0)
## DeltaE
try:
	DeltaE = float(sys.argv[7])
except:
	print('Error: insert real DeltaE >= 0\n')
	sys.exit(0)
if DeltaE < 0:
	print('Error: insert real DeltaE >= 0\n')
	sys.exit(0)
## nLanczos
try:
	nLanczos = int(sys.argv[8])
except:
	print('Error: insert integer nLanczos > 0\n')
	sys.exit(0)
if nLanczos%1 != 0 or nLanczos <= 0:
	print('Error: insert integer nLanczos > 0\n')
	sys.exit(0)

# open writing file
fw = open("results_excE-S2/s" + str(s) + "_N" + str(N) + "_J" + str(J)
	+ "meV_beta" + str(beta) + "_BC" + str(BC) + "_Sz" + str(Sz) + "_DeltaE"
	+ str(DeltaE) + "meV_nLanczos" + str(nLanczos) + ".txt", "w")

# basis
if (2*s)%2 == 0:
	basis = spin_basis_1d(N, m=Sz/N, S=str(int(s)), pauli=False)
else:
	basis = spin_basis_1d(N, m=Sz/N, S=str(int(2*s)) + '/2', pauli=False)

t1 = time.time() - start_time

# Hamiltonian
H = H_BLBQ_1D(N,J,beta,BC,basis)

t2 = time.time() - start_time

# diagonalization
En,psin = H.eigsh(k=nLanczos, which='SA')

t3 = time.time() - start_time

# nmax
nmax = len(En)
for n in range(1,len(En)):
	if En[n]-En[0] > DeltaE:
		nmax = n
		break
if nmax==len(En):
	print("Warning: larger nLanczos is required\n")

# excitation energies
excEn = [En[n]-En[0] for n in range(1,nmax)]

# S^2 operator
S2op = S2_op(N,basis)

# <psin|S^2|psin>
S2n = [np.dot(psin[:,n].conj(),S2op.dot(psin[:,n])) for n in range(nmax)]

t4 = time.time() - start_time

# outputs
for n in range(nmax):
	fw.write("#E" + str(n) + " = " + str(En[n]) + " meV\n")
fw.write("--------------------\n\n")
fw.write("#List of excitation energies (meV):\n")
fw.write(str(excEn))
fw.write("\n\n")
fw.write("#List of S^2:\n")
fw.write(str(S2n))
fw.write("\n\n")
fw.write("--------------------\n")
fw.write("#time to initialize and find basis = " + str(t1) + " s\n")
fw.write("#time to build Hamiltonian = " + str(t2-t1) + " s\n")
fw.write("#time to diagonalize = " + str(t3-t2) + " s\n")
fw.write("#time to compute excitation energies and S^2 = " + str(t4-t3)
	+ " s\n")
fw.write("#total time = " + str(time.time() - start_time) + " s\n")

## close file
fw.close()
###############################################################################
