# Electric Dipole Spin Resonance Simulator
A simulation of EDSR in a carbon nanotube Double Quantum Dot using kwant package.

Theory based on '_Coherent Single Electron Spin Control in a Slanting Zeeman Field_' by Tokura et. al.

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

* ```evolve-output/``` - the files produced by the the time evolution of the simulator.
  * Each file is named according the time of the simulation was initiated.
  * Each file is a json object with the following general structure
  ```
  {{  
            'B_0': #,
            'lattice_points': ,
            'length': s,
            'eV_0': self.hartree_to_ev(data['eV_0']),
            'E_sl': self.hartree_to_ev(data['E_sl']),
            'E_x': self.hartree_to_ev(data['E_x']),
            'E_z': self.hartree_to_ev(data['E_z']),
            'perturbation': data['perturbation'],
            'potential_type': data['potential_type'],
            'effective_mass': data['effective_mass'],
            'states': 
        }
  ```
* 


## Variables


## Systen Class

### Class Functions

### Class Variables
