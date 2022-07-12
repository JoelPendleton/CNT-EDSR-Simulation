# Electric Dipole Spin Resonance Simulator
A simulation of EDSR in a carbon nanotube Double Quantum Dot using kwant package.

Theory based on '_Coherent Single Electron Spin Control in a Slanting Zeeman Field_' by Tokura et. al.

* https://kwant-project.org
* https://tkwant.kwant-project.org

## My Development Environment
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

## Example Plot

Animation produced by simulation

![PDF](https://user-images.githubusercontent.com/11929366/159930136-39123816-2229-4735-9ed7-acb3d638ef2d.mp4)



