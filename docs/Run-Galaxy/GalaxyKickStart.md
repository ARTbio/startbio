## Installation of a Galaxy server with Ansible and the GalaxyKickStart playbook

### What is Ansible ?

![Ansible](images/ansible.png)

Ansible is an automation engine that automates configuration management and application
deployment.

Ansible reads instructions (Tasks) from a playbook and performs the indicated tasks on
target machines (Hosts), through an ssh connection.

There is no magics: everything an "administrator" can do using command lines of a linux OS,
can be automated with ansible that "wraps" these command lines.
The power of Ansible (and similar orchestration software, ie Puppet, Chief, etc.) comes
from the abstraction of complex suite of commands in the Ansible syntax.
Moreover, automation allows to reproduce exactly the desired configuration.
Finally, Ansible is `idempotent`: whatever the initial configuration, it brings the target
to the exact same final state. This is useful to repair a broken configuration.

### Ansible playbook - GalaxyKickStart

The Ansible "language" (Striclty speaking, Ansible language is *not* a programming language)
is structured. Thus a playbook is not necessarily a single flat file. Multiple tasks can be gathered in a file, a "role" is the execution of a set of tasks, and a playbook can execute multiple roles.
 
GalaxyKickStart is an Ansible playbook that will

- install basic dependencies needed for Galaxy
- Create and manage all the linux users involved in the deployment of Galaxy
- Install and configure the services required for Galaxy:
    - postgresql (database engine)
    - nginx (web server)
    - docker (containers)
    - proftpd (ftp server)
    - slurm (job manager)
    - supervisor (service manager)
- Configure Galaxy for using these services
- Install tools and workflows using the [bioblend API](https://github.com/galaxyproject/bioblend).

The code of the GalaxyKickStart playbook is freely available at the ARTbio GitHub
Repository [https://github.com/ARTbio/GalaxyKickStart](https://github.com/ARTbio/GalaxyKickStart).

----
### Deployment

- start a GCE VM with the following characteristics

!!! info "Google Instance for Ansible deployment"
    - Name: `ansible-galaxy`
    - Region `europe-west6 (Zurich)` (or any region available with you Google coupon)
    - Zone: `europe-west6-a` (or `-b` or `-c`)
    - **Configuration de la machine**
        - `OPTIMISEE POUR LE CALCUL`
        - Série: `C2`
        - Type de machine: `c2-standard-16 (16 processeurs virtuels, 64 Go de mémoire)`
    - **Disque de démarrage (Modifier)**
        - `IMAGES PUBLIQUES`
        - Système d'exploitation: `Ubuntu`
        - Version*: `Ubuntu 20.04 LTS`
        - Type de disque de démarrage: `Disque persistant avec équilibrage`
        - Taille (Go): `200`
        - ==SELECTONNER==
    - **Pare-feu**
        - Check `Autoriser le trafic HTTP`

- connect to you VM using the Google ssh console
- start an interactive session as root using the command
```
sudo -i
```
- download the script `run_ansible_analyse_genomes_2019.sh` using the command
```
wget https://raw.githubusercontent.com/ARTbio/Run-Galaxy/master/deployment_scripts/run_ansible_analyse_genomes_2019.sh
```

- We are now ready to run this script. However this year there is bonus ! All trainees will participate to ==Pasteur 2019
Ansible Racing==.
In order to participate, you'll just have to put the `time` command just before the script invokation, as follows:

```
time sh run_ansible_analyse_genomes_2019.sh analyseGenomes_2019
```

!!! danger "The Ultimate Pasteur 2019 Ansible Racing"
    Please copy the time info returned by your console at the end of the deploymment.
    It shoud look like this:
    ```
    real	10m27.142s
    user	8m22.941s
    sys	1m16.409s
    ```
    Then Paste this time as a comment in this [GitHub issue](https://github.com/ARTbio/Run-Galaxy/issues/25)

- When the deployment is finished, connect to your ansible-deployed "GalaxyKickStart" instance:
    
    Just click on the url displayed in your Google Cloud Engine Console.
    
- Connect to your server as an admin:

    This time, ansible and the GalaxyKickStart playbook already programmatically registered
    an admin user. Just use the `admin@galaxy.org:admin` as credentials (user:password)
    
    When logged in, see that required tools as well as workflows are already installed !

!!! warning
    admin is not really a decent password,
    please ==c h a n g e  .  y o u r  .  p a s s w o r d==
    to avoid your Galaxy server getting hacked before the end of the course.

----
### Galaxy administration tasks
#### Transfert input data to you newly deployed Galaxy instance
that is :
    - a data set with reference sequences
    - a data set with small RNAseq files
    - a data set with RNAseq files

- Click the main menu `User` --> `Saved Histories`
- Press the top right button (above history list) `Import from file`
- copy this url :
```
https://galaxy.pasteur.fr/history/export_archive?id=4c5da5ad7355ff42
```
    
- repeat the same operation with: 
```
https://galaxy.pasteur.fr/history/export_archive?id=eb4c1d5564c9f78c
```
and
```
https://galaxy.pasteur.fr/history/export_archive?id=69a1b70d1c4a6bdb
```

??? bug "In case of emergency"
    
    ### NEXTCLOUD
    
    Small RNAseq
    ```
    https://usegalaxy.sorbonne-universite.fr/nextcloud/index.php/s/wC4DrxHN3gLtx4c/download
    ```
    Reférences
    ```
    https://usegalaxy.sorbonne-universite.fr/nextcloud/index.php/s/Sri6HKiwCn4RbSq/download
    ```
    RNAseq
    ```
    https://usegalaxy.sorbonne-universite.fr/nextcloud/index.php/s/4FWrWxZf72KDNji/download
    ```
    
    ### Amazon S3
    RNAseq
    ```
    https://mydeepseqbucket.s3.amazonaws.com/RNAseq
    ```
    small RNAseq
    ```
    https://mydeepseqbucket.s3.amazonaws.com/smallrnaseqs
    ```
    references
    ```
    https://mydeepseqbucket.s3.amazonaws.com/references
    ```
    
---
#### Content of the `run_ansible_analyse_genomes_2019.sh` script

``` bash
#!/usr/bin/env bash
set -e
apt update -y
apt install -y python-pip python-dev python-setuptools git htop
echo "Upgrading pip"
pip install -U pip
/usr/local/bin/pip --version
/usr/local/bin/pip install ansible==2.7.4
ansible --version
git clone https://github.com/ARTbio/GalaxyKickStart.git -b $1
cd GalaxyKickStart/
ansible-galaxy install -r requirements_roles.yml -p roles/ -f
cp scripts/8cpu_job_conf.xml roles/galaxyprojectdotorg.galaxy-extras/templates/job_conf.xml.j2
cp scripts/configure_slurm.py.j2 roles/galaxyprojectdotorg.galaxy-extras/templates/configure_slurm.py.j2
ansible-playbook -i inventory_files/analyseGenomes galaxy.yml
echo "end of deployment\n"
```

??? info "the `run_ansible_analyse_genomes_2019.sh` script explained"
    1. The shebang line (`#!`) says that it is a script code that has to be executed
    by the shell bash which can be found in the /usr/bin/env environment
    2. `set -e` says to the bash interpreter to exit the run at first error (to avoid catastrophes)
    3. update apt package database
    4. installs `python-pip`, `python-dev`, `python-setuptools` (these 3 packages are required to
    install pip), `git` (to clone and manage GitHub repositories) and `htop` (a monitoring tool)
    using the package installer `apt`
    5. Is just a command to inform the user about run state. This will prompt
    "Upgrading pip version" in the console
    6. does what is echoed before by the previous line : this is the command to upgrade the pip program that was
    installed with installation of `python-pip`.
    `pip` is a recursive acronym that can stand for either "Pip Installs Packages" or
    "Pip Installs Python".
    7. will prompt the version of pip in the console
    8. install `ansible`, version 2.7.4, using `pip`
    8. will prompt the version of ansible in the console
    10. clone the GalaxyKickStart Repository available at https://github.com/ARTbio/GalaxyKickStart.git,
    creating locally the `GalaxyKickStart` folder. The repository `branch` that is cloned is indicated
    as a parameter in the command line (the `$1`).
    11. Changes directory, i.e. goes to /root/GalaxyKickStart
    12. Says to ansible to install additional roles (collection of files to control ansible)
    which are not the the GalaxyKickStart repository but whose address is stated in the file
    `requirements_roles.yml`. These roles will be installed in the subdirectory
    `/root/GalaxyKickStart/roles/`. NB: `ansible-galaxy` has *nothing* to do with Galaxy,
    the name of this ansible command is serendipitous.
    13. triggers the play of the playbook `galaxy.yml` by ansible. The target host of the playbook
    is defined in the file `inventory_files/analyseGenomes`, as well as how ansible will interact with the target.
    Here, we play the playbook on the same computer (localhost).
    14. Prompts the end of the deployment
