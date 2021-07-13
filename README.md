# ED_BLBQ
Python codes to solve bilinear biquadratic (BLBQ) spin models by exact diagonalization (ED).

# System requirements
These codes use the open-source QuSpin Python package (https://weinbe58.github.io/QuSpin/), whose installation guide can be found at https://weinbe58.github.io/QuSpin/Installation.html. 
The recommended installation uses Anaconda/Miniconda to install Python3 and manage all package dependencies.

# Hamiltonian
The spin Hamiltonian under consideration reads as

<a href="https://www.codecogs.com/eqnedit.php?latex=H_{\textup{BLBQ}}&space;=&space;J&space;\sum_i&space;\left[&space;\bold{S}_i&space;\cdot&space;\bold{S}_{i&plus;1}&space;&plus;&space;\beta&space;\left(&space;\bold{S}_i&space;\cdot&space;\bold{S}_{i&plus;1}&space;\right&space;)^2&space;\right&space;]," target="_blank"><img src="https://latex.codecogs.com/gif.latex?H_{\textup{BLBQ}}&space;=&space;J&space;\sum_i&space;\left[&space;\bold{S}_i&space;\cdot&space;\bold{S}_{i&plus;1}&space;&plus;&space;\beta&space;\left(&space;\bold{S}_i&space;\cdot&space;\bold{S}_{i&plus;1}&space;\right&space;)^2&space;\right&space;]," title="H_{\textup{BLBQ}} = J \sum_i \left[ \bold{S}_i \cdot \bold{S}_{i+1} + \beta \left( \bold{S}_i \cdot \bold{S}_{i+1} \right )^2 \right ]," /></a>

where J is the exchange coupling, 
<a href="https://www.codecogs.com/eqnedit.php?latex=\beta" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\beta" title="\beta" /></a>
is a parameter that determines the relative strength between the bilinear and biquadratic terms, and
<a href="https://www.codecogs.com/eqnedit.php?latex=\bold{S}_i" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\bold{S}_i" title="\bold{S}_i" /></a>
denote spin-s operators at site i.
The sum over i runs from 1 to N-1 (N mod N) for open (periodic) boundary conditions, where N is the size of the spin chain.

# Code usage

## Types of code
There are 3 types of code available, which compute different quantities.

1) The code 'excE_S2.py' computes excitation energies with respect to the ground state, and the expectation value of the operator
<a href="https://www.codecogs.com/eqnedit.php?latex=S^2&space;=&space;\bold{S}&space;\cdot&space;\bold{S}," target="_blank"><img src="https://latex.codecogs.com/gif.latex?S^2&space;=&space;\bold{S}&space;\cdot&space;\bold{S}," title="S^2 = \bold{S} \cdot \bold{S}," /></a>
with
<a href="https://www.codecogs.com/eqnedit.php?latex=\bold{S}&space;=&space;\sum_{i=1}^N&space;\bold{S}_i," target="_blank"><img src="https://latex.codecogs.com/gif.latex?\bold{S}&space;=&space;\sum_{i=1}^N&space;\bold{S}_i," title="\bold{S} = \sum_{i=1}^N \bold{S}_i," /></a>
for the corresponding eigenstates.

2) The code 'Szi.py' computes the expectation value of the operator
<a href="https://www.codecogs.com/eqnedit.php?latex=S^z_i" target="_blank"><img src="https://latex.codecogs.com/gif.latex?S^z_i" title="S^z_i" /></a>
for eigenstates.

3) The codes 'ssw_S0GS.py' and 'ssw_S1GS.py' compute the spin spectral weights of excitations from a ground state with total spin S=0 and S=1, respectively. 
The spin spectral weight at site i for the state M' is formally defined as

<a href="https://www.codecogs.com/eqnedit.php?latex=\mathcal{S}_{M'}(i)&space;=&space;\sum_M&space;P_M&space;\sum_{\alpha=x,y,z}&space;|\langle&space;M&space;|&space;S^a_i&space;|&space;M'&space;\rangle&space;|^2," target="_blank"><img src="https://latex.codecogs.com/gif.latex?\mathcal{S}_{M'}(i)&space;=&space;\sum_M&space;P_M&space;\sum_{\alpha=x,y,z}&space;|\langle&space;M&space;|&space;S^a_i&space;|&space;M'&space;\rangle&space;|^2," title="\mathcal{S}_{M'}(i) = \sum_M P_M \sum_{\alpha=x,y,z} |\langle M | S^a_i | M' \rangle |^2," /></a>

where 
<a href="https://www.codecogs.com/eqnedit.php?latex=P_M" target="_blank"><img src="https://latex.codecogs.com/gif.latex?P_M" title="P_M" /></a>
denotes the equilibrium occupation of the state M. 
This calculation is done at zero temperature, so that 
<a href="https://www.codecogs.com/eqnedit.php?latex=P_M" target="_blank"><img src="https://latex.codecogs.com/gif.latex?P_M" title="P_M" /></a>
is the inverse of the ground state degeneracy for the ground states, and zero for all other eigenstates.
It must also be noted that a selection rule imposes that the spin spectral weights can only be non-zero if the total spin of the ground state and the state M' differ from 0,+-1.

## Code inputs

### Physical parameters

All codes take as inputs the following "physical parameters": s, N, J, 
<a href="https://www.codecogs.com/eqnedit.php?latex=\beta" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\beta" title="\beta" /></a> 
and a variable BC that can be set to 0 (1) for open (periodic) boundary conditions.

### Other parameters

In addition, the codes also have as inputs some of the following "other parameters": Sz denotes the 
<a href="https://www.codecogs.com/eqnedit.php?latex=S^z" target="_blank"><img src="https://latex.codecogs.com/gif.latex?S^z" title="S^z" /></a>
symmetry sector in which the calculations are done, DeltaE sets the energy range for the excitations, nLanczos is the number of Lanczos vectors used in the diagonalization routine, and neigen is the number of eigenstates considered.

## Example of usage

1) Excitation spectrum, within the range of 50 mev, of an N=5 spin-1 open chain, with J = 18.4 meV and 
<a href="https://www.codecogs.com/eqnedit.php?latex=\beta&space;=&space;0.085" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\beta&space;=&space;0.085" title="\beta = 0.085" /></a> (note: for spin-1 chains, we know that we can obtain all energy levels in the symmetry sector with 
<a href="https://www.codecogs.com/eqnedit.php?latex=S^z&space;=&space;0" target="_blank"><img src="https://latex.codecogs.com/gif.latex?S^z&space;=&space;0" title="S^z = 0" /></a>).
```
python excE_S2.py 1 5 18.4 0.085 0 0 50 10
```

2) Same as in 1 for a cyclic chain.
```
python excE_S2.py 1 5 18.4 0.085 1 0 50 10
```

3) Average magnetization of the first 3 eigenstates with 
<a href="https://www.codecogs.com/eqnedit.php?latex=S^z&space;=&space;&plus;1," target="_blank"><img src="https://latex.codecogs.com/gif.latex?S^z&space;=&space;&plus;1," title="S^z = +1," /></a> for the same model as in 1.
```
python Szi.py 1 5 18.4 0.085 0 1 3
```

4) Spin spectral weigths of all allowed excitations within the range of 50 meV, for the same model as in 1 (note: we already know that this model gives a ground state with S=1).
```
python ssw_S1GS.py 1 5 18.4 0.085 0 50 30
```

5) Spin spectral weigths of all allowed excitations within the range of 50 meV, for the same model as in 2 (note: we already know that this model gives a ground state with S=0).
```
python ssw_S0GS.py 1 5 18.4 0.085 1 50 30
```
