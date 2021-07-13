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

3) The codes 'ssw_S0GS.py' and 'ssw_S1GS.py' compute the spin spectral weight for transitions from a ground state with total spin S=0 and S=1, respectively. 
The spin spectral weight at site i for the state M' is formally defined as

<a href="https://www.codecogs.com/eqnedit.php?latex=\mathcal{S}_{M'}(i)&space;=&space;\sum_M&space;P_M&space;\sum_{\alpha=x,y,z}&space;|\langle&space;M&space;|&space;S^a_i&space;|&space;M'&space;\rangle&space;|^2," target="_blank"><img src="https://latex.codecogs.com/gif.latex?\mathcal{S}_{M'}(i)&space;=&space;\sum_M&space;P_M&space;\sum_{\alpha=x,y,z}&space;|\langle&space;M&space;|&space;S^a_i&space;|&space;M'&space;\rangle&space;|^2," title="\mathcal{S}_{M'}(i) = \sum_M P_M \sum_{\alpha=x,y,z} |\langle M | S^a_i | M' \rangle |^2," /></a>

where 
<a href="https://www.codecogs.com/eqnedit.php?latex=P_M" target="_blank"><img src="https://latex.codecogs.com/gif.latex?P_M" title="P_M" /></a>
denotes the equilibrium occupation of the state M. 
This calculation is done at zero temperature, so that 
<a href="https://www.codecogs.com/eqnedit.php?latex=P_M" target="_blank"><img src="https://latex.codecogs.com/gif.latex?P_M" title="P_M" /></a>
is the inverse of the ground state degeneracy for the ground states, and zero for all other eigenstates.

## Code inputs
All codes take as inputs the following parameters: J, 
<a href="https://www.codecogs.com/eqnedit.php?latex=\beta" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\beta" title="\beta" /></a>, 
s, N and a variable BC that should be set to 0 (1) for open (periodic) boundary conditions.
These variables are called "physical parameters" and can be changed in the source codes.

In addition, all codes also take as inputs the number of excited states ('nexc') and the number of Lanczos vectors for the diagonalization routine ('nLanczos'). These variables are called "other parameters" and can also be changed in the source codes.

Finally, the codes "te_S2.py" and "localSz.py" also allow to choose the 
<a href="https://www.codecogs.com/eqnedit.php?latex=S^z" target="_blank"><img src="https://latex.codecogs.com/gif.latex?S^z" title="S^z" /></a>
sector in which the calculations are done. This variable ('Sz') is also in the group of "other parameters" and can be changed in the source codes.

## Example of usage
We now provide a representative example of usage of all codes.
We consider an N=9 spin-1 (s=1) chain, with J=18.4 [meV] and beta=0.085.

First, we run the "te_S2.py" code for both BC=0,1. 
For spin-1 chains, we know that all eigenstates can be obtained within the Sz=0 sector, and therefore we set Sz=0.
For BC=0 (1), we also set nexc=50, nLanczos=51 (nexc=10, nLanczos=11).

