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
docker run joelpendleton/cnt-edsr-simulation
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

### My Data

The results of my simulations used in my Master's thesis can be found on [Google Drive](https://drive.google.com/drive/folders/1aA9u9uZPkQBIRVceL2X89pe_fS5aaqBI?usp=sharing)


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
* <b>Arguments</b>:

    * `hamiltonian`: the full hamiltonian to describe the system.
    * `pertubation_type`: the type of time-dependent potential (either "sin" or "cos").
    * `magnetic_field_file`: path to simulated/real magnetic field file.
    * `number_of_lattices`: the number of lattice points to use in the simulation.
    * `potential_type`: the type of confinement potential to use (either 0, infinite square well or 1, parabolic potential). 


### Class Functions

Within the class there are a variety of functions. They are described here



* `cosine_v_ac(time, z, eV_0_au, pulse_frequency_au, total_length_au)`
</br>
</br> Function to define the alternating potential associated the alternating electric field as a cosine function. 
    * <b>Arguments</b>:
        * `time`: the current time step of the simulation
        * `z`: the positions along the nanotube
        * `eV_0_au`: Constant defining the energy of the alternating electric field (atomic units), i.e. the strength of the electric field.
        * `pulse_frequency_au`: the frequency of oscillation of the alternating potential.
        * `total_length_au`: the total length of the nanotube in atomic units.
    * <b>Returns:</b> the alternating potential at a given time across the nanotube.
</br></br>


* `sine_v_ac(time, z, eV_0_au, pulse_frequency_au, total_length_au)` </br> 
</br>Function to defines the alternating potential associated the alternating electric field as a sine function. 
    * <b>Arguments</b>:
        * `time`: the current time step of the simulation
        * `z`: the positions along the nanotube
        * `eV_0_au`: Constant defining the energy of the alternating electric field (atomic units), i.e. the strength of the electric field.
        * `pulse_frequency_au`: the frequency of oscillation of the alternating potential.
        * `total_length_au`: the total length of the nanotube in atomic units.
    * <b>Returns</b>: the alternating potential at a given time across the nanotube.
</br></br>

* `tesla_to_au(tesla)` </br> 
</br> Function to convert the magnetic field strength from teslas to atomic units. 
    * <b>Arguments</b>:
        * `tesla`: the magnetic field strength in teslas to be converted.
    * <b>Returns:</b> the magnetic field strength in atomic units.
</br></br>


* `au_to_tesla(au)` </br> 
</br>Function to convert the magnetic field strength from atomic units to teslas.
    * <b>Arguments</b>:
        * `au`: the magnetic field strength in atomic units to be converted.
    * <b>Returns:</b> the magnetic field strength in teslas.
</br></br>


* `second_to_au(time)` </br> 
</br> Function to convert the time the simulation has been running from seconds to atomic units.
    * <b>Arguments</b>:
        * `time`: the time the simulation has evolved for in seconds.
    * <b>Returns:</b> the time the simulation has evolved for in atomic units.
</br></br>



* `au_to_second(time)`  </br> 
</br>Function to convert the time the simulation has been running from atomic units to seconds.
    * <b>Arguments</b>:
        * `time`: the time the simulation has evolved for in atomic units.
    * <b>Returns:</b> the time the simulation has evolved for in seconds.
</br></br>



* `hartree_to_ev(hartree)` </br> 
</br> Function to convert energies from hartrees (atomic units) to electron volts.
    * <b>Arguments</b>:
        * `hartree`: the energy in hartrees.
    * <b>Returns:</b> the energy in electron volts.
</br></br>


* `ev_to_hartree_(ev)` </br> 
</br> Function to convert energies from electron volts to hartrees (atomic units).
    * <b>Arguments</b>:
        * `ev`: the energy in electron volts.
    * <b>Returns:</b> the energy in hartrees.
</br></br>


* `au_to_m(au)` </br> 
</br> Function to convert lengths from atomic units to metres.
    * <b>Arguments</b>:
        * `au`: the length in atomic units.
    * <b>Returns:</b> the length in metres.
</br></br>

* `m_to_au(m)` </br> 
</br> Function to convert lengths from metres to atomic units.
    * <b>Arguments</b>:
        * `m`: the length in metres.
    * <b>Returns:</b> the length in atomic units.
</br></br>


* `hz_to_au(hz)` </br> 
</br> Function to convert a frequency from hertz to atomic units.
    * <b>Arguments</b>:
        * `hz`: the frequency in hertz.
    * <b>Returns:</b> the frequency in atomic units.

* `import_mumax3_simulations()` </br></br>
Function to import mumax3 simulation files for use within the EDSR simulator.
    * <b>Returns:</b> 
    `B_x, B_y, B_z` the x,y and z components of the imported magnetic field.
</br></br>

 
* `potential(z, time)` </br></br>
Function to combine a confinement potential, which is either a parabolic or hardwall potential, with the value for the alternating potential associated with the alternating electric field (at a given time).
    * <b>Arguments</b>:
        * `z`: the positions along the carbon nanotube.
        * `time`: the time the simulation has evolved for in atomic units.

    * <b>Returns:</b>  `total_potential` the total potential experienced along the carbon nanotube.
</br></br>


* `kwant_shape(self, site)` function to define the shape of the scattering region of nanotube,
    * <b>Arguments</b>:
        * `site`: the current site.

    * <b>Returns:</b> a boolean which is true if the scattering site should be drawn.
</br></br>


* `make_system(self)` </br></br>
Function to build the system in kwant and define the initial state of a variety of variables.
    * <b>Returns:</b> `self.syst` the built system object produced by kwant.
</br></br>


* `eigenstates(self)` </br></br>
Function to compute the eigenstates, i.e. eigenvalues and eigenvectors for the defined system at some current time.
    * <b>Returns:</b> `eigenvalues, eigenvectors` the eigenvalues and eigenvectors associated with the system at some time.
</br></br>


* `eigenstates(self)` </br></br>
Function to compute the eigenstates, i.e. eigenvalues and eigenvectors for the defined system at some current time.
    * <b>Returns:</b> `eigenvalues, eigenvectors` the eigenvalues and eigenvectors associated with the system at some time.
</br></br>


* `initial_pdfs(self)` </br></br>
Function to show the initial probability density functions of the spin-up and spin-down ground state (at time = 0). Saves plot to results folder.
    * <b>Returns:</b> `density_1, density_2` the initial probability density functions of the spin-up and spin-down ground state.
</br></br>

* `initial_energies(self)` </br></br>
Function to display initial energy levels of the system in eV.  Saves plot to results folder.
    * <b>Returns:</b> `y` a list of the eigenenergies of the system at time = 0.
</br></br>


* `evolve(self, time_steps=1000)` </br></br>
Function to evolve/simulate the defined system through the pi rotation time.  Saves data to evolve-output folder.
    * <b>Arguments</b>:
        * `time_steps`: the number of time divisions/steps to use for the simulation.

    * <b>Returns:</b> a boolean which is true if the simulation is complete.
</br></br>


* `visualise(self, file_name)` </br></br>
Function to visualise / produce plots of the evolved simulation.  Saves plots and animations to results folder.
    * <b>Arguments</b>:
        * `file_name`: the file in which the simulation/evolved data is stored.

    * <b>Returns:</b> a boolean which is true if the visualisation is complete.
</br></br>



### Class Variables

These class variables can be overwritten before evolving (simulating) to give your system the desired characteristics.

* `E_sl_eV`: The energy associated with the slanted magnetic field (electron volts).
* `omega_0_eV`: The energy spacing for the parabolic/harmonic potential (electron volts).
* `m_SI`: The effective mass of the electron (kg).
* `B_0_SI`: The external magnetic field, applied in the z-direction, used to induced zeeman splitting (T)

<b>N.B.</b> The System object could be modified so that the above class variables are instead defined by passing arguments in the creation of the System class.


## Example of how to use System Class to Perform a Simulation

### Evolution

Here is an example you would evolve the system: 

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
```

### Visualisation

Suppose the above simulation produces an evolution file titled *[20220309-135220.json](./20220309-135220.json)*. The following command will allow you to visualise the output of this evolution.

```python
system.visualise("20220309-135220")
```

The visualisationproduces plots in `/results` folder along with an animation. An example of the animation is shown below:

https://user-images.githubusercontent.com/11929366/181641119-0a29b718-36b5-45cf-8f36-b9b28647e5b9.mp4
