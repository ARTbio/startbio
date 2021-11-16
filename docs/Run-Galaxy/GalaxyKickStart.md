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
- download the script `run_ansible_analyse_genomes_2021.sh` using the command:

??? warning "for Girls :woman:"
    ## Only for :woman:
    ```
    wget https://raw.githubusercontent.com/ARTbio/Run-Galaxy/master/deployment_scripts/run_ansible_analyse_genomes_2021-F.sh
    ```
??? warning "for Boys :man:"
    ## Only for :man:
    ```
    wget https://raw.githubusercontent.com/ARTbio/Run-Galaxy/master/deployment_scripts/run_ansible_analyse_genomes_2021-M.sh
    ```

- We are now ready to run this script. In addition, all trainees may participate to the ==Pasteur 2021
Ansible Racing==.
In order to participate, you'll just have to put the `time` command just before the script invokation, as follows:

```
time sh run_ansible_analyse_genomes_2021*.sh
```

!!! danger "The Ultimate Pasteur 2021 Ansible Racing"
    Please copy the time info returned by your console at the end of the deploymment.
    It shoud look like this:
    ```
    real    37m23.924s
    user    17m26.569s
    sys     2m33.091s
    ```
    Then Paste this time as a comment in this [GitHub issue](https://github.com/ARTbio/Run-Galaxy/issues/30)

- When the deployment is finished, connect to your ansible-deployed "GalaxyKickStart" instance:
    
    Just click on the url displayed in your Google Cloud Engine Console.
    
- Connect to your server as an admin:

    This time, ansible and the GalaxyKickStart playbook already programmatically registered
    an admin user. Just use the `admin@galaxy.org:artbio2020` as credentials (user:password)
    
    When logged in, see that required tools as well as workflows are already installed !

!!! warning
    artbio2020 is not really a decent password,
    please ==c h a n g e  .  y o u r  .  p a s s w o r d==
    to avoid your Galaxy server getting hacked before the end of the course 😉
