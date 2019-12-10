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
- Install tools and workflows using the bioblend API.

The code of the GalaxyKickStart playbook is freely available at the ARTbio GitHub
Repository [https://github.com/ARTbio/GalaxyKickStart](https://github.com/ARTbio/GalaxyKickStart).

----
### Deployment

- start a GCE VM

!!! question "Google Instance"
    - **Série N2** --> `n2-standard-16 (16 processeurs virtuels, 64 Go de mémoire)`
    - **Disque de démarrage**  --> `Ubuntu 16.04 LTS`
    - **Taille (Go)** --> `200`
    - **Pare-feu** --> `Autoriser le trafic HTTP`

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
    
    ### Google Bucket
    RNAseq
    ```
    https://00e9e64bacea836fa47259198b2eab01124b904383636ff0aa-apidata.googleusercontent.com/download/storage/v1/b/artbio_genomic_analysis/o/export_archive%3Fid=4c5da5ad7355ff42?qk=AD5uMEsXTlqAydH0E_PKot_OEqkjwXrHZ8YtHLY3VQEYKtvaCNI65E8TwgKqD7r3IyGsuCsGYv0eY5i2tLrsA3UBUjLcj9anAslkSK7VTBKd0-_sDFsnv6cCdn2-ayIRtZLSVmMZHrLhaO0fUJLLx-FTb6VrdoUAXTcUSSZ5tw81ZV7XyWUvZatF68LRq_ysKLT8DC-H7JrT5FTsxUMLGojAmw1S0cMon9YgKh6-mOJIiSSPd1nFbgZHJmWG_4BZbARI_i_OgyDuwOvY9mAjJkP4z2vGST9kXhssCNHic7l1CzAtAtdBOWlnuCFhoSRJh-7UBtEJR-mzV0MxZJMApjwzgr6rZyUH7dyS5CQoxQ9bvAD7ht_nSPpaZyso9upXip3E9gLzlS58RntFFCWXZwC8GGZAZsBRPbrHVe5gpjTEPqOOG-E9x53D9Z81iXIg8VG19FXJoC0f6ckERtVLKcYcZqwX9M0Q7N9Pa915bHjU5XfWTIvX1_ZM7C2xtBsr-U92-PndloBe0RxSyrm-SEZPJNWDoQuCi50rSllyPxhx34m8NrCvYb3ybtwYAPDPqLajJHOifn-Noyd-j4D1pdBC71knirlSDZGPiye-xoiHjla41kclx_au_yHa5xylKGKVk7zS5ej8HrY4eVs6RrBD3UtJm_EltN7ryxeiGk4EeyS3ukneVjPZ60PvAN-7ETXYuEWNMir8ahMb-61v0MRkrLPGdz6ulze2d0sps9iz-iE0weDwsTYPaRpP-ZagBLBbhKMJHYHYAwAEEEgPqJQlCi6xAtLAohhqZsBvDrEOgAwqaUp3MV97oNuwDWq5FDogxSmg6uW32HC71U6x6BUQNqG0a6iQ_A
    ```
    small RNAseq
    ```
    https://00e9e64bac97b3729b9bfd29fb4f1ccf9f00fd36c006db080d-apidata.googleusercontent.com/download/storage/v1/b/artbio_genomic_analysis/o/export_archive%3Fid=eb4c1d5564c9f78c?qk=AD5uMEuLpiyeaRoOVUKn8Z9U7Bl-4EIywzZiqJikUmtRGzcccvax_2y4UFjj58AePA8Ti9bmlSC5zNxZgst-eo9eBiRD0FgtDOt04EZR9CBuVi2DGNH5uKsshZon4FZmj6_eSK-m-sVWHD4ADOOHOL_abOml2Y3BovzcIHM-rpAyXOmJl01xKDS9SgiwzU_y490EiWKN8Fx4Dld7JKGPtaUu_u0yqbYlIq3-pEKXlpNqaIDUHPdZsfWAzebA5Vin2D3WKkyq8a-Ikf6RBXr7wlh60LRna-lBK6o3wE-CQlFp9_t7MFgh1MAJ2_tFhNeozScOm60ICoFCWendWVBfo_BZBEg1gyodA1gBj-vmn9hV_DzRpqe0ZtGscBIgZszvx2PRuPNnmyyJ6SyR3IPwBME_K4wy7eG7OCCcP124eu48Kpx-xTsU0d9L3AfaS1WAm2h7G2UuHcmIWYq5npx_7Rh2ZZgy1z8JVyjSifzkPj9oeQUhEFrBtlSGw3yvKJytovj7c04pPDVIlWHslXUxKP4ndNbHUkt-IfGanAtXgGpDKGjLw5ultm1RR_n1NfvJ2CogpcIaH1vXlh_6z3Err6v83B2b9sc-goxjXIqF8kdT2-7I96f2T89R4lhz2mppkk7ag6vhf_UjMS-uhkCvnO1NHRFAL8T-PWy0ZCSC4OwdzoLUDrkXHrk3oba83g6MGLhWcgz075qL5KwDL9oY0ATILwU6VfOFV9lV65OcelhwFQ_TKOKmy95oTKqRWvX0cWlpAEUd_wHURsHGm403Sz9dvfMtYa2nsH97HEJ7u_k2BXr6Mb-UtluaGXF9KSh6nPUS7MB7Jm4XIfGAtB3VftpcLN3af9egxw
    ```
    references
    ```
    https://00e9e64bac20148f28e34cd40b4ae13615baaa028616851e2f-apidata.googleusercontent.com/download/storage/v1/b/artbio_genomic_analysis/o/export_archive%3Fid=69a1b70d1c4a6bdb?qk=AD5uMEubfUik_l0Qbfnp3EeO7ZW_kp7vUwVv_9grqwUl-c2AEU5BdftnS7n01EzwapUYDIjplZmOekyD16E5jYWIHYMczwZh-kytmVFmdV4pXzaEji2enWAFhhNQIX6t77Ixfu26ZtD8nfOT8JD4aQL02VyDHvpYb79pv7xLBFKWmPBigOP1Al0JDgEbhUgKvxeKBHtY62prehQAtzocv3Dctl_t5-hReT_ZySrNAU7rKirEvw1r4p7HZLxevno558zhms3hvGoldqfYNyQhPwtoPzvkYf4qHQFkgKpFE6PkOzrbph8UUEpTJg0O8IhVT2XRL-jPUiHiWae-JOgZdjWhduTEDncDINH8zOK7XXPN0nfh-qisau2zOdOlAK_evRP_X7ZrNob64th1hyoEJLPKY0MVcmKctjlEAU6SjT4KA3Xn3XriU38WyJUjQGLqJ9cU0pBug8n4Arfd4nmni9kerS1mcuzsmBS2m18QhG-KJ2NBcnhl79Dae9h1OrHw2C6n-_LkhqGfuj-Ms62FG1btMRp0FepmrhnG16mrxLWWApOzdb9phFS-MsaHU92ij_IZ3Q72AoLEGfMcnrUloPfutWSBI4TxroCXtjVjqZBviwSv3WnN53KnVWAXrpT4XH8re_A_NoeIr0OxYeakCx23UH5WXpd3vKHlRxfiL5LYW6YK329S4rvIS2bLU_tR5_RJbeSWG8-eWH5_h-hEWZT-0T3lAP_8CM3tNmcqcGa7AfrHOgxToskaX2Y7OCnZ7ohOPQDgsQHakuR8uq5M9RO3W-rCt77lYEXyqVxl7FuGxN8hqsnbaPxSnf4r8_3r--pujo3VI93Yu44ebGPhEBS35jyQg3-p1w
    ```
    
    Mississippi
    
    ```
    https://mississippi.snv.jussieu.fr/history/export_archive?id=a407203305eb0400
    ```
    ```
    https://mississippi.snv.jussieu.fr/history/export_archive?id=8916e356920c765e
    ```
    ```
    https://mississippi.snv.jussieu.fr/history/export_archive?id=e6cf496eef19f5b8
    ```
    
    ARTbio
    
    ```
    https://artbio.snv.jussieu.fr/history/export_archive?id=b1d0a8de1661cd49
    ```
    ```
    https://artbio.snv.jussieu.fr/history/export_archive?id=d1ad3dbae5fee6ca
    ```
    ```
    https://artbio.snv.jussieu.fr/history/export_archive?id=b4076b0ca537f521
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
