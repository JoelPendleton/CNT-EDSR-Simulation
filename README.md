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

## Example Plot

Animation of the probability density function as time increases 

![PDF](https://github.com/JoelPendleton/CNT-EDSR-Simulation/blob/master/results/20220309-135220/animation.mp4)



