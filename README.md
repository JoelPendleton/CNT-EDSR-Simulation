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
      

## System Class

This class is used to define the system object upon which time-evolution simulations can be executed.

The `__init__` function
* <u>Arguments</u>:

    * `hamiltonian`: the full hamiltonian to describe the system.
    * `pertubation_type`: the type of time-dependent potential (either "sin" or "cos").
    * `magnetic_field_file`: path to simulated/real magnetic field file.
    * `number_of_lattices`: the number of lattice points to use in the simulation.
    * `potential_type`: the type of confinement potential to use (either 0, infinite square well or 1, parabolic potential). 


### Class Functions

Within the class there are a variety of functions. They are described here

* `cosine_v_ac(self, time, z, eV_0_au, pulse_frequency_au, total_length_au)` this function defines the alternating potential associated the alternating electric field as a cosine function. 
    * <u>Arguments</u>:

        * `time`: the current time step of the simulation
        * `z`: the positions along the nanotube
        * `eV_0_au`: Constant defining the energy of the alternating electric field (atomic units), i.e. the strength of the electric field.
        * `pulse_frequency_au`: the frequency of oscillation of the alternating potential.
        * `total_length_au`: the total length of the nanotube in atomic units.
    * <u>Returns:</u> the alternating potential at a given time across the nanotube.

* `sine_v_ac(self, time, z, eV_0_au, pulse_frequency_au, total_length_au)` this function defines the alternating potential associated the alternating electric field as a sine function. 
    * <u>Arguments</u>:

        * `time`: the current time step of the simulation
        * `z`: the positions along the nanotube
        * `eV_0_au`: Constant defining the energy of the alternating electric field (atomic units), i.e. the strength of the electric field.
        * `pulse_frequency_au`: the frequency of oscillation of the alternating potential.
        * `total_length_au`: the total length of the nanotube in atomic units.
    * <u>Returns:</u> the alternating potential at a given time across the nanotube.


* `tesla_to_au(self, tesla)` this function converts the magnetic field strength from teslas to atomic units.
    * <u>Arguments</u>:
        * `tesla`: the magnetic field strength in teslas to be converted.
    * <u>Returns:</u> the magnetic field strength in atomic units.


* `au_to_tesla(self, au)` this function converts the magnetic field strength from atomic units to teslas.
    * <u>Arguments</u>:
        * `au`: the magnetic field strength in atomic units to be converted.
    * <u>Returns:</u> the magnetic field strength in teslas.

* `second_to_au(self, time)` this function converts the time the simulation has been running from seconds to atomic units.
    * <u>Arguments</u>:
        * `time`: the time the simulation has evolved for in seconds.
    * <u>Returns:</u> the time the simulation has evolved for in atomic units.



* `au_to_second(self, time)` this function converts the time the simulation has been running from atomic units to seconds.
    * <u>Arguments</u>:
        * `time`: the time the simulation has evolved for in atomic units.
    * <u>Returns:</u> the time the simulation has evolved for in seconds.


* `hartree_to_ev(self, hartree)` this function converts energies from hartrees (atomic units) to electron volts.
    * <u>Arguments</u>:
        * `hartree`: the energy in hartrees.
    * <u>Returns:</u> the energy in electron volts.

* `ev_to_hartree_(self, ev)` this function converts energies from electron volts to hartrees (atomic units).
    * <u>Arguments</u>:
        * `ev`: the energy in electron volts.
    * <u>Returns:</u> the energy in hartrees.


* `au_to_m(self, au)` this function converts lengths from atomic units to metres.
    * <u>Arguments</u>:
        * `au`: the length in atomic units.
    * <u>Returns:</u> the length in metres.

* `m_to_au(self, m)` this function converts lengths from metres to atomic units.
    * <u>Arguments</u>:
        * `m`: the length in metres.
    * <u>Returns:</u> the length in atomic units.


* `hz_to_au(self, hz)` this function converts a frequency from hertz to atomic units.
    * <u>Arguments</u>:
        * `hz`: the frequency in hertz.
    * <u>Returns:</u> the frequency in atomic units.

* `import_mumax3_simulations(self)` this function imports mumax3 simulation files for use within the EDSR simulator.
    * <u>Returns:</u> 
    `B_x, B_y, B_z` the x,y and z components of the imported magnetic field.

 
* `potential(self, z, time)` this function combines a confinement potential, which is either a parabolic or hardwall potential, with the value for the alternating potential associated with the alternating electric field (at a given time).
    * <u>Arguments</u>:
        * `z`: the positions along the carbon nanotube.
        * `time`: the time the simulation has evolved for in atomic units.

    * <u>Returns:</u> `total_potential` the total potential experienced along the carbon nanotube.


* `kwant_shape(self, site)` function to define the shape of the scattering region of nanotube,
    * <u>Arguments</u>:
        * `site`: the current site.

    * <u>Returns:</u> a boolean which is true if the scattering site should be drawn.

* `make_system(self)` function to build the system in kwant and define the initial state of a variety of variables.
    * <u>Returns:</u> `self.syst` the built system object produced by kwant.

* `eigenstates(self)` function to compute the eigenstates, i.e. eigenvalues and eigenvectors for the defined system at some current time.
    * <u>Returns:</u> `eigenvalues, eigenvectors` the eigenvalues and eigenvectors associated with the system at some time.

* `eigenstates(self)` function to compute the eigenstates, i.e. eigenvalues and eigenvectors for the defined system at some current time.
    * <u>Returns:</u> `eigenvalues, eigenvectors` the eigenvalues and eigenvectors associated with the system at some time.


* `initial_pdfs(self)` Function to show the initial probability density functions of the spin-up and spin-down ground state (at time = 0). Saves plot to results folder.
    * <u>Returns:</u> `density_1, density_2` the initial probability density functions of the spin-up and spin-down ground state.

* `initial_energies(self)`Function to display initial energy levels of the system in eV.  Saves plot to results folder.
    * <u>Returns:</u> `y` a list of the eigenenergies of the system at time = 0.


* `evolve(self, time_steps=1000)` function to evolve/simulate the defined system through the pi rotation time.  Saves data to evolve-output folder.
    * <u>Arguments</u>:
        * `time_steps`: the number of time divisions/steps to use for the simulation.

    * <u>Returns:</u> a boolean which is true if the simulation is complete.


* `visualise(self, file_name)` function to visualise / produce plots of the evolved simulation.  Saves plots and animations to results folder.
    * <u>Arguments</u>:
        * `file_name`: the file in which the simulation/evolved data is stored.

    * <u>Returns:</u> a boolean which is true if the visualisation is complete.



### Class Variables

These class variables can be overwritten before evolving (simulating) to give your system the desired characteristics.

* `E_sl_eV`  the energy associated with the slanted magnetic field (electron volts).
* `omega_0_eV`  the energy spacing for the parabolic/harmonic potential (electron volts).
* `m_SI`  the effective mass of the electron (kg).
* `B_0_SI`  the external magnetic field, applied in the z-direction, used to induced zeeman splitting (T)

Note: the System objec could be modified so that the above class variables are instead define by passing arguments in the creation of the System class.


## Example of how to use System Class

```python
lattices = 100  # number of lattice points

potential = 0  # infinite square-well potential
# magnetic_field_file = "B_eff000000.npy" 
system = System("((A * k_z**2) + V(z, time)) * identity(2) + B(z) * sigma_z + C(z) * sigma_x + D(z) * sigma_y",
                pertubation_type="sin", number_of_lattices=lattices,
                potential_type=potential)  # call system object 
system.make_system()

# Plot the initial energies and probability density functions
# Run these both before you evolve.

system.initial_energies()
system.initial_pdfs()

# Evolve the simulation for 100 time steos
system.evolve(100)

# file_name = None
# Visualise the evolution by providing the file name of the generated output
# system.visualise(file_name)

```