# CNT EDSR-Simulation
A simulation of EDSR in a carbon nanotube using kwant package.

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

## Example Plots

* Animation of the probability density function as time increases 

https://user-images.githubusercontent.com/11929366/149535321-635bcd71-e6c5-4bd5-aa7d-193fbba6de93.mp4



