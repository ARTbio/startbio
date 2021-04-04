follow this [procedure](https://jupyterhub.readthedocs.io/en/stable/installation-guide-hard.html?highlight=add%20python%20module#install-jupyterhub-and-jupyterlab-from-the-ground-up)

Then if you wish to code with your own specific kernel, build it as a conda environment

----

Create a user kernel for jupyterhub:

1. Navigate in your conda environment
2. `conda create -n myconda python=3.7.5 ipykernel`
3. check your envs:
    
    `conda info -e`
    
4. upgrade your conda env as you like:
    
    `conda install python=3.8.5 # for instance`

5. inform jupyterhub of new potential kernel
    ```
    .conda/envs/myconda/bin/python -m ipykernel install --user --name 'chris-myconda' --display-name "Chris myconda Env"
    ```
6. check that the new kernel is there:
      
    `jupyter kernelspec list`
7. in case you wish to remove and dereference your kernel:
    
    `jupyter kernelspec uninstall chris-myconda`
