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

# Sx(i) operator
def Sxi_op(i,basis):
	Sp = [[1/2,i]]
	Sm = [[1/2,i]]

	static = [
		['+',Sp],
		['-',Sm]
	]
	dynamic = []

	no_checks = dict(check_pcon=False,check_symm=False,check_herm=False)
	Sxi = hamiltonian(static,dynamic,basis=basis,dtype=np.float64,**no_checks)

	return Sxi

# ISy(i) operator
def ISyi_op(i,basis):
	Sp = [[1/2,i]]
	Sm = [[-1/2,i]]

	static = [
		['+',Sp],
		['-',Sm]
	]
	dynamic = []

	no_checks = dict(check_pcon=False,check_symm=False,check_herm=False)
	ISyi = hamiltonian(static,dynamic,basis=basis,dtype=np.float64,**no_checks)

	return ISyi

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

# spin spectral weight GS(S=1) -> S0
def calc_ssw_S1S0(EnGS,psiGS,EnS0,psinS0,Nsites,basis):
	Sziop = [Szi_op(i,basis) for i in range(Nsites)]
	Sxiop = [Sxi_op(i,basis) for i in range(Nsites)]
	ISyiop = [ISyi_op(i,basis) for i in range(Nsites)]

	excE_nS0 = [0.0 for nS0 in range(len(EnS0))]
	sswi_nS0 = [[0.0 for i in range(Nsites)] for nS0 in range(len(EnS0))]
	for nS0 in range(len(EnS0)):
		excE_nS0[nS0] = EnS0[nS0] - EnGS

		for i in range(Nsites):
			for jGS in range(3):

				#|<GS|Sz(i)|psinS0>|^2
				sswi_nS0[nS0][i] += 1/3*abs( np.dot(psiGS[jGS].conj(),
					Sziop[i].dot(psinS0[nS0])) )**2

				#|<GS|Sx(i)|psinS0>|^2
				sswi_nS0[nS0][i] += 1/3*abs( np.dot(psiGS[jGS].conj(),
					Sxiop[i].dot(psinS0[nS0])) )**2

				#|<GS|ISy(i)|psinS0>|^2
				sswi_nS0[nS0][i] += 1/3*abs( np.dot(psiGS[jGS].conj(),
					ISyiop[i].dot(psinS0[nS0])) )**2

	return excE_nS0, sswi_nS0

# spin spectral weight GS(S=1) -> S1
def calc_ssw_S1S1(EnGS,psiGS,EnS1,psinS1,Nsites,basis):
	Sziop = [Szi_op(i,basis) for i in range(Nsites)]
	Sxiop = [Sxi_op(i,basis) for i in range(Nsites)]
	ISyiop = [ISyi_op(i,basis) for i in range(Nsites)]

	excE_nS1 = [0.0 for nS1 in range(len(EnS1)//3)]
	sswi_nS1 = [[0.0 for i in range(Nsites)] for nS1 in range(len(EnS1)//3)]
	for nS1 in range(len(EnS1)//3):
		excE_nS1[nS1] = EnS1[nS1*3] - EnGS

		for i in range(Nsites):
			for jGS in range(3):
				for j in range(3):

					#|<GS|Sz(i)|psinS1>|^2
					sswi_nS1[nS1][i] += 1/3*abs( np.dot(psiGS[jGS].conj(),
						Sziop[i].dot(psinS1[nS1*3+j])) )**2

					#|<GS|Sx(i)|psinS1>|^2
					sswi_nS1[nS1][i] += 1/3*abs( np.dot(psiGS[jGS].conj(),
						Sxiop[i].dot(psinS1[nS1*3+j])) )**2

					#|<GS|ISy(i)|psinS1>|^2
					sswi_nS1[nS1][i] += 1/3*abs( np.dot(psiGS[jGS].conj(),
						ISyiop[i].dot(psinS1[nS1*3+j])) )**2

	return excE_nS1, sswi_nS1

# spin spectral weight GS(S=1) -> S2
def calc_ssw_S1S2(EnGS,psiGS,EnS2,psinS2,Nsites,basis):
	Sziop = [Szi_op(i,basis) for i in range(Nsites)]
	Sxiop = [Sxi_op(i,basis) for i in range(Nsites)]
	ISyiop = [ISyi_op(i,basis) for i in range(Nsites)]

	excE_nS2 = [0.0 for nS2 in range(len(EnS2)//5)]
	sswi_nS2 = [[0.0 for i in range(Nsites)] for nS2 in range(len(EnS2)//5)]
	for nS2 in range(len(EnS2)//5):
		excE_nS2[nS2] = EnS2[nS2*5] - EnGS

		for i in range(Nsites):
			for jGS in range(3):
				for j in range(5):

					#|<GS|Sz(i)|psinS1>|^2
					sswi_nS2[nS2][i] += 1/3*abs( np.dot(psiGS[jGS].conj(),
						Sziop[i].dot(psinS2[nS2*5+j])) )**2

					#|<GS|Sx(i)|psinS1>|^2
					sswi_nS2[nS2][i] += 1/3*abs( np.dot(psiGS[jGS].conj(),
						Sxiop[i].dot(psinS2[nS2*5+j])) )**2

					#|<GS|ISy(i)|psinS1>|^2
					sswi_nS2[nS2][i] += 1/3*abs( np.dot(psiGS[jGS].conj(),
						ISyiop[i].dot(psinS2[nS2*5+j])) )**2

	return excE_nS2, sswi_nS2
###############################################################################

# main #
start_time = time.time()

# read inputs
if len(sys.argv) != 8:
	print('Error: run code as #python_<code.py>_<s>_<N>_<J(meV)>_<beta>_<BC>'
		+ '_<DeltaE(meV)>_<nLanczos>\n')
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
## DeltaE
try:
	DeltaE = float(sys.argv[6])
except:
	print('Error: insert real DeltaE >= 0\n')
	sys.exit(0)
if DeltaE < 0:
	print('Error: insert real DeltaE >= 0\n')
	sys.exit(0)
## nLanczos
try:
	nLanczos = int(sys.argv[7])
except:
	print('Error: insert integer nLanczos > 0\n')
	sys.exit(0)
if nLanczos%1 != 0 or nLanczos <= 0:
	print('Error: insert integer nLanczos > 0\n')
	sys.exit(0)

# open writing file
fw = open("results_ssw-S1GS/s" + str(s) + "_N" + str(N) + "_J" + str(J)
	+ "meV_beta" + str(beta) + "_BC" + str(BC) + "_DeltaE" + str(DeltaE)
	+ "meV_nLanczos" + str(nLanczos) + ".txt", "w")

# Szlist for ssw with S=1 GS
Szlist = [-2,-1,0,+1,+2]

# Nup list from Szlist
Nup = [int(Sz+N*s) for Sz in Szlist]

# basis
if (2*s)%2 == 0:
	basis = spin_basis_1d(N, Nup=Nup, S=str(int(s)), pauli=False)
else:
	basis = spin_basis_1d(N, Nup=Nup, S=str(int(2*s)) + '/2', pauli=False)

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

# S^2 operator
S2op = S2_op(N,basis)

# <psin|S^2|psin>
S2n = [np.dot(psin[:,n].conj(),S2op.dot(psin[:,n])) for n in range(nmax)]

# ssw for S=1 GS
## GS
EnGS, psiGS = En[0], [psin[:,0],psin[:,1],psin[:,2]]
### error alerts
for j in range(3):
	if abs(S2n[j]-2) > 1e-8:
		print('Error: GS triad does not have S=1\n')
		sys.exit(0)
	if abs(EnGS-En[j]) > 1e-8:
		print('Error: GS triad is not degenerate\n')
		sys.exit(0)
if abs(EnGS-En[3]) < 1e-8:
	print('Error: S=1 GS is not unique\n')
	sys.exit(0)
## relevant states (S=0, S=1, S=2)
EnS0, psinS0 = [], []
EnS1, psinS1 = [], []
EnS2, psinS2 = [], []
for n in range(nmax):
	if abs(S2n[n]-0) < 1e-8:
		psinS0.append(psin[:,n])
		EnS0.append(En[n])
	if abs(S2n[n]-2) < 1e-8:
		psinS1.append(psin[:,n])
		EnS1.append(En[n])
	if abs(S2n[n]-6) < 1e-8:
		psinS2.append(psin[:,n])
		EnS2.append(En[n])
### error alerts
if len(EnS1)%3 != 0:
	print('Error: number of S=1 states is not multiple of 3\n')
	sys.exit(0)
for nS1 in range(len(EnS1)//3):
	for j in range(1,3):
		if abs(EnS1[nS1*3]-EnS1[nS1*3+j]) > 1e-8:
			print('Error: triad nr ' + str(nS1) + ' of S=1 states is not ' +
				'degenerate\n')
			sys.exit(0)
if len(EnS2)%5 != 0:
	print('Error: number of S=2 states is not multiple of 5\n')
	sys.exit(0)
for nS2 in range(len(EnS2)//5):
	for j in range(1,5):
		if abs(EnS2[nS2*5]-EnS2[nS2*5+j]) > 1e-8:
			print('Error: quintet nr ' + str(nS2) + ' of S=2 states is not ' +
				'degenerate\n')
			sys.exit(0)
## list of te, ssw(i)
excE_nS0, sswi_nS0 = calc_ssw_S1S0(EnGS,psiGS,EnS0,psinS0,N,basis)
excE_nS1, sswi_nS1 = calc_ssw_S1S1(EnGS,psiGS,EnS1,psinS1,N,basis)
excE_nS2, sswi_nS2 = calc_ssw_S1S2(EnGS,psiGS,EnS2,psinS2,N,basis)

t4 = time.time() - start_time

# outputs
for n in range(nmax):
	fw.write("#E" + str(n) + " = " + str(En[n]) + " meV\n")
fw.write("\n")
for n in range(nmax):
	fw.write("#S^2 for state " + str(n) + " = " + str(S2n[n]) + "\n")
fw.write("--------------------\n\n")
fw.write("#List of excitation energies for GS (S=1) -> S=0 (meV):\n")
fw.write(str(excE_nS0))
fw.write("\n")
fw.write("#List of spin spectral weights for GS (S=1) -> S=0:\n")
fw.write(str(sswi_nS0))
fw.write("\n\n")
fw.write("#List of excitation energies for GS (S=1) -> S=1 (meV):\n")
fw.write(str(excE_nS1))
fw.write("\n")
fw.write("#List of spin spectral weights for GS (S=1) -> S=1:\n")
fw.write(str(sswi_nS1))
fw.write("\n\n")
fw.write("#List of excitation energies for GS (S=1) -> S=2 (meV):\n")
fw.write(str(excE_nS2))
fw.write("\n")
fw.write("#List of spin spectral weights for GS (S=1) -> S=2:\n")
fw.write(str(sswi_nS2))
fw.write("\n\n")
fw.write("--------------------\n")
fw.write("#time to initialize and find basis = " + str(t1) + " s\n")
fw.write("#time to build Hamiltonian = " + str(t2-t1) + " s\n")
fw.write("#time to diagonalize = " + str(t3-t2) + " s\n")
fw.write("#time to compute spin spectral weights = " + str(t4-t3) + " s\n")
fw.write("#total time = " + str(time.time() - start_time) + " s\n")

## close file
fw.close()
###############################################################################
