# ED_BLBQ
Python codes to solve bilinear-biquadratic spin models.

# System requirements
These codes use the open-source QuSpin Python package (https://weinbe58.github.io/QuSpin/), whose installation guide can be found at https://weinbe58.github.io/QuSpin/Installation.html. 
The recommended installation uses Anaconda/Miniconda to install Python3 and manage all package dependencies.

# Hamiltonian
The spin Hamiltonian under consideration, referred to as bilinear-biquadratic (BLBQ) model, reads as

<a href="https://www.codecogs.com/eqnedit.php?latex=H&space;=&space;J&space;\sum_{i}&space;\left[&space;\bold{S}_i&space;\cdot&space;\bold{S}_{i&plus;1}&space;&plus;&space;\beta&space;\left(&space;\bold{S}_i&space;\cdot&space;\bold{S}_{i&plus;1}&space;\right&space;)^2&space;\right&space;]," target="_blank"><img src="https://latex.codecogs.com/gif.latex?H&space;=&space;J&space;\sum_{i}&space;\left[&space;\bold{S}_i&space;\cdot&space;\bold{S}_{i&plus;1}&space;&plus;&space;\beta&space;\left(&space;\bold{S}_i&space;\cdot&space;\bold{S}_{i&plus;1}&space;\right&space;)^2&space;\right&space;]," title="H = J \sum_{i} \left[ \bold{S}_i \cdot \bold{S}_{i+1} + \beta \left( \bold{S}_i \cdot \bold{S}_{i+1} \right )^2 \right ]," /></a>

where J is the exchange coupling, 
<a href="https://www.codecogs.com/eqnedit.php?latex=\beta" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\beta" title="\beta" /></a>
is a parameter that determines the relative strength between bilinear and biquadratic terms, and
<a href="https://www.codecogs.com/eqnedit.php?latex=\bold{S}_i" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\bold{S}_i" title="\bold{S}_i" /></a>
denotes spin-S operators at site i.
The sum over i runs from 1 to N-1 (N) for open (periodic) boundary conditions, where N is the size of the spin chain.
