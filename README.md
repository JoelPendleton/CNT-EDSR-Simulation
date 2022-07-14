# Electric Dipole Spin Resonance Simulator
A simulation of EDSR in a carbon nanotube Double Quantum Dot using kwant package.

Computations are performed in atomic units, and are only converted to SI units 
through the execution of the visualise function.

Theory largely based on '_Coherent Single Electron Spin Control in a Slanting Zeeman Field_' by Tokura et. al.

* https://kwant-project.org
* https://tkwant.kwant-project.org

## Installation 
**docker images: joelpendleton/cnt-edsr-simulation**   

I recommend using my docker image by running the following:

```
docker pull joelpendleton/cnt-edsr-simulation
docker start joelpendleton/cnt-edsr-simulation
conda activate env-tkwant
```

Alternatively, you run the following installation commands:

### Linux and MacOS

```
conda create -n env-tkwant python=3.7
conda activate env-tkwant
conda install tkwant -c conda-forge
```

### Windows
```
conda create -n env-tkwant python=3.7
conda activate env-tkwant
conda install -c intel mpi4py
conda install tkwant -c conda-forge
```

## Directories

* ```evolve-output/```: the files produced by the the time evolution of the simulator.
  * Each file is named according the time of the simulation was initiated.
  * Each file is a json object with the following keys.
    * `B_0`: the applied external static field (atomic units).
    * `lattice_points`: the number of lattice points.
    * `length`: the length of the carbon nanotube (atomic units)
    * `E_sl`: the energy associated with the slanted magnetic field (atomic units).
    * `E_x`: the energy shift associated with the apparent oscillation of the slanted magnetic field (atomic units).
    * `E_delta`: the energy difference between the lowest two eigenstates (atomic units).
    * `perturbation`: the type of time-dependent perturbation used (sin or cos).
    * `potential_type`: the type of confining potential used (infinite-well or parabolic).
    * `effective_mass`: effective mass of electrons within quantum dot (atomic units)
    * `E_1`: the lowest energy eigenstate (atomic units).
    * `E_2`: the second lowest energy eigenstate (atomic units).
    * `states`: a list containing the states of the system at each time step during the evolution. </br> Each state contains the following keys:  
      * `time`: the time associated with the evolution step (atomic units).
      * `pdf`: a list containing the values of the probability density function over the lattice point coordinates.
      * `B_x`: effective magnetic field experienced by electrons in x-direction (atomic units).
      * `B_y`: effective magnetic field experienced by electrons in y-direction (atomic units).
      * `B_z`: effective magnetic field experienced by electrons in z-direction (atomic units).
      * `rho_sx`: expectation of pauli spin x matrix.
      * `rho_sy`: expectation of pauli spin y matrix.
      * `rho_sz`: expectation of pauli spin z matrix.
  * Example: 
  ```json
    {
      "B_0": 2.127659574468085e-8,
      "lattice_points": 4,
      "length": 12472.192422511227,
      "eV_0": 0.0003674930360069677,
      "E_sl": 0.000036749303600696764,
      "E_x": 0.000002088290610937897,
      "E_z": 4.203467974491095e-8,
      "perturbation": "sin",
      "potential_type": "infinite-well",
      "effective_mass": 1,
      "E_1": 1.0018601870742315e-8,
      "E_2": 5.205328161565326e-8,
      "states": [
        {
          "time": 0,
          "pdf": [
            4.0560596446038007e-10,
            6.476893426371965e-9,
            3.268146788485447e-8,
            1.0281407512141786e-7
          ],
          "B_x": 1.599900632463728e-9,
          "B_y": 1.121243705358237e-8,
          "B_z": 7.930829912939263e-8,
          "rho_sx": 0.006612777555865185,
          "rho_sy": -0.11966290116446851,
          "rho_sz": 0.9926646474106606
        },
        {
          "time": 30391.612126166216,
          "pdf": [
            4.0212530049091855e-10,
            6.422358439157277e-9,
            3.2412154508827166e-8,
            1.0198482108839007e-7
          ],
          "B_x": 1.5962850145527839e-9,
          "B_y": 1.12150843771851e-8,
          "B_z": 7.927541051153356e-8,
          "rho_sx": 0.006598986611589326,
          "rho_sy": -0.11966109969031674,
          "rho_sz": 0.9925629062783262
        }]
  }
  ```
         
* `results/`: contains files produced by `visualise` function 
  * `parameters.json`: contains parameters used in simulation converted into interpretable units.
    * The keys in the json file are as follows
      * `B_0`: B_0`: the applied external static field (Teslas)
      * `lattice_points`: : the number of lattice points
      * `length`:  the length of the carbon nanotube (metres)
      * `eV_0`: Constant defining the energy of the alternating electric field (electron volts).
      * `E_sl`: the energy associated with the slanted magnetic field (electron volts).
      * `E_x`: the energy shift associated with the apparent oscillation of the slanted magnetic field (electron volts).
      * `E_delta`: the energy difference between the lowest two eigenstates (electron volts).
      * `perturbation`: the type of time-dependent perturbation used (sin or cos).
      * `potential_type`: the type of confining potential used (infinite-well or parabolic).
      * `effective_mass`: effective mass of electrons within quantum dot (atomic units).
  * `animation.mp4`: animation of the time evolution of the system.
  * `magnetic-fields.eps`: plot of magnetic fields vs time.
  * `spin-expectations.eps`: plot of spin matrix expectation values vs time.
      


## Variables


## Systen Class

### Class Functions

### Class Variables
