## Training Menu

Three methods of Galaxy server deployment are explained in this tutorial, which can be
used with personal computers, clusters of machines or virtual machines in cloud computing
environment.

All you need is **an ssh access** and **the root control** of the target machine.
These two conditions are far more easily fulfilled with virtual machines in clouds. This is
the reason why we are going to use
[virtual machines in the Google Cloud Engine](spin_off_VM.md) and/or [virtual machines in
the IFB core Cloud (Biosphere)](https://biosphere.france-bioinformatique.fr/cloud/). 

### Outline of the training session
-----

##### 1. Install a minimal standalone galaxy server in the `Google Cloud Platform`
  - Install a computational workflow
  - Install tools for proper execution of the workflow
  - Running the workflow.

##### 2. Install a Galaxy server with `Ansible` and the `GalaxyKickStart playbook`
##### 3. Use case of Galaxy administration
  We are going to import dataset in the server, uncompress them, manipulate collections,
  etc...
  This will pave the way to the subsequent analyses which you will have to perform.
  
##### 4. Deployment of a Galaxy server using `Docker`
  This is optional. We will do it if we have time (unlikely), or if you are eager to do
  it !
