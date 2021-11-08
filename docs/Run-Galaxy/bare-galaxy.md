!!! info "Differences between Google Cloud Engine and IFB Cloud"
    The procedure to start a VM is not the same whether you are using Google Cloud Engine
    or the IFB Cloud.
    If you are working with GCE, follow 1-option A and 2-option A.
    Whereas if you are working with the IFB Cloud, follow 1-option B and 2-optionB


### 1A. Spin off a virtual Machine `bare-galaxy` with ![](images/google-padok.png){: style="width:30px"} Google Cloud Engine

??? tip "Tip"
    Start and stop of Google Virtual Machines is also described in [Annex 2](spin_off_VM.md)

- Connect to your Google Compute Instances
  [dashboard](https://console.cloud.google.com/compute/instances)

- Create a Virtual Machine Instance
 

!!! question "with the following settings"
    - Name: `bare-galaxy`
    - Region `europe-west6 (Zurich)` (or any region available with you Google coupon)
    - Zone: `europe-west6-a` (or `-b` or `-c`)
    - **Configuration de la machine**
        - `OPTIMISEE POUR LE CALCUL`
        - Série: `C2`
        - Type de machine: `c2-standard-8 (8 processeurs virtuels, 32 Go de mémoire)`
    - **Disque de démarrage (Modifier)**
        - `IMAGES PUBLIQUES`
        - Système d'exploitation: `Ubuntu`
        - Version*: `Ubuntu 20.04 LTS`
        - Type de disque de démarrage: `Disque persistant avec équilibrage`
        - Taille (Go): `100`
        - ==SELECTONNER==
    - **Pare-feu**
        - Check `Autoriser le trafic HTTP`

This settings should look like:
    
![](images/GCE_spin.png){: style="width:450px"}
![](images/GCE_OS.png){: style="width:450px"}
![](images/GCE_firewall.png){: style="width:450px"}

### 1B. Spin off a virtual Machine `bare-galaxy` with the ![](images/biosphere.png){: style="width:70px"} core-IFB cloud

??? tip "Tip"
    Start and stop of IFB-Cloud Machines is also described in [Annex 3](spin_off_IFB-VM.md)

- Connect to your the [biosphere](https://biosphere.france-bioinformatique.fr/), and click
  on **RAINBio** menu.

- Choose the ![](images/ubuntu_rainbio.png){: style="width:200px"} virtual image.
- Choose "Déploiement Avancé" in the menu "Lancer"

![](images/deploiement_avance.png){: style="width:200px"}

- Give a name to your VM, choose `IFB-core` as a cloud region, `ifb.m4.2xlarge
  (8 vcpu, 32Go GB RAM, 470Go GB local disk)` for the machine, and press `Lancer`.
  
  ![](images/appliance_IFB.png)
 
- Follow the deployment of your VM in the `myVM` menu. In contrast to the Google Cloud platform,
  this may take more that 10 min.

### 2A. Connect to the VM using the ssh web console

!!! question "ssh connection"
    Roll down the `ssh` menu in the control pannel and select the first option
    `Ouvrir dans une fenêtre du navigateur`

    ![Select ssh session in browser](images/select_ssh.png)
    Here is your web ssh console to control your VM
    ![](images/web_ssh_console.png)

### 2B. SSH connect to the IFB VM using your terminal

- Be sure that your `private` key (`mykey`) is in your ~/.ssh/folder.
- The corresponding `public`key (`mykey.pub`) _should have been deposited/uploaded to biosphere_,
  otherwise, it cannot work.
- Type the following command
```
ssh -A ubuntu@134.158.247.168 # replace by the IP of your deployed VM
```
- If this command does not work (it happens...), type instead:
```
ssh -i ~/.ssh/mykey ubuntu@134.158.247.168 # ! mykey, NOT mykey.pub. And replace by your IP
```
- You should see a shell in your connected VM, which looks like:
??? info "Terminal"
    ```
    imac-chris:~ aligre$ ssh -A ubuntu@134.158.247.168
    The authenticity of host '134.158.247.168 (134.158.247.168)' can't be established.
    ECDSA key fingerprint is SHA256:WdN9NuYfDgj0DM0r78fH5rUkSwuQK3IIH+H4FmkGpOM.
    Are you sure you want to continue connecting (yes/no/[fingerprint])? yes
    Warning: Permanently added '134.158.247.168' (ECDSA) to the list of known hosts.
    Welcome to Ubuntu 20.04.3 LTS (GNU/Linux 5.4.0-88-generic x86_64)
    
     * Documentation:  https://help.ubuntu.com
     * Management:     https://landscape.canonical.com
     * Support:        https://ubuntu.com/advantage
    
      System information as of Mon Nov  8 18:29:33 UTC 2021
    
      System load:  0.07               Processes:                153
      Usage of /:   17.0% of 19.21GB   Users logged in:          0
      Memory usage: 3%                 IPv4 address for docker0: 172.17.0.1
      Swap usage:   0%                 IPv4 address for ens3:    10.158.16.9
    
    0 updates can be applied immediately.
    
    
    *** System restart required ***
    
    The programs included with the Ubuntu system are free software;
    the exact distribution terms for each program are described in the
    individual files in /usr/share/doc/*/copyright.
    
    Ubuntu comes with ABSOLUTELY NO WARRANTY, to the extent permitted by
    applicable law.
    ```



**==From there, the procedure is the same, whether with the Google VM or the IFB VM==**

### 3. Installation of the Galaxy server

In this first approach "==bare-galaxy==", everything is made super simple:

- We are gonna become `root` unix user. This is easier because installation
of new programs as well as manipulations of network interfaces is generally permitted only
to users with administration rights.

- We are gonna check that all software needed to deploy galaxy are there (they are with
Ubuntu 20.04 !)

- Finally, we will run the automated deployment of Galaxy

So let's do this, step by step:

  1.
    
  ```
  sudo -i
  ```
  This command open a new "shell" where you are root. You can check this by typing `pwd` that
  should return `/root/`, meaning that you are now working in the directory of the `root` user.
  
  2.
  ```
  python3 --version && git --version && nano --version
  ```
  This command checks that the only 2 programs required for the deployment are already there
  
  3.
  ```
  git clone https://github.com/galaxyproject/galaxy.git -b release_21.05
  ```
  This command says to use `git` to `clone` the code repository located at
  `https://github.com/galaxyproject/galaxy.git`.
  
  In addition the `-b release_21.05` option specifies that only the version `release_21.05`
  will be cloned locally in your virtual machine. You may try to visualize the URL
  [https://github.com/galaxyproject/galaxy.git](https://github.com/galaxyproject/galaxy.git)
  in your web browser. You will, literally, see the code of Galaxy. It is Open Source, as
  you can notice.
  
  4.
  ```
  cd galaxy
  ```
  This command shift you in the `galaxy` directory that was created by git and the
  `git clone` command in 3.
  
  5.
  ```
  cp config/galaxy.yml.sample config/galaxy.yml
  ```
  This command makes a copie of the `galaxy.yml.sample` file to `galaxy.yml` - in the
  directory `config` that is in the `galaxy` directory.
  
  6.
  ```
  nano config/galaxy.yml
  ```
  Using this command, we are going to edit some important settings that are required to
  run our Galaxy fresh instance.
!!! question ""
    - Find the line
    ```
    http: 127.0.0.1:8080
    ```
    (you can use the editor command ++ctrl+w++, paste the previous line and press enter)
    
    and edit it to
    ```
    http: 0.0.0.0:80
    ```
    By doing this, we ensure that we will be able to reach the galaxy web server on our
    virtual machine using the usual web port `80`.

    - Find the line
    ```
    #admin_users: ''
    ```
    delete the `#` character and type your email address between the two single quotes.
    
    Any email address is ok (admin@galaxy.org for instance). It is just used here as
    an admin identifier.
    
    - save your changes by pressing the key combination ++ctrl+o++
    - quit nano by pressing the key combination ++ctrl+x++

??? tip "OPTIONAL but SAVING US 20 min of deployment !"
    Before starting the deployment of Galaxy, we are going to use a little trick to
    bypass the step of compilation of html and javascript codes which are used to
    render the Galaxy graphic interface.
    
    Indeed, modern web application use a lot of cached codes that speed up
    the user experience. However, this implies that this code cache is built during the
    deployment of the application.
    
    For Galaxy, building/caching the client codes for the web server takes about 20 min
    and this is increasing with newer galaxy versions.
    
    To save us these 20 min, we are going to remove the web client folders and replace them
    by already built client folders, prepared by your trainer...
    
    1. Remove the web client folders
    ```
    rm -rf ~/galaxy/client ~/galaxy/static
    ```
    2. Download the cached web client folders
    ```
    cd ~/galaxy && wget https://mydeepseqbucket.s3.amazonaws.com/client.tar.gz https://mydeepseqbucket.s3.amazonaws.com/static.tar.gz
    ```
    3. Uncompress the cached client folders
    ```
    cd ~/galaxy && tar -xvf static.tar.gz && tar -xvf client.tar.gz
    ```
    
    Last note: this tip is TOTALLY optional, if you run the next command without doint it,
    everything will go the same, but it will take more time...
  7.
  Ready for deploying Galaxy ?
  
??? bug "One more thing if you are using an IFB VM"
    In IFB cloud, VM instances have a very small system volume, and we have installed
    the galaxy git repository on this volume.
    You can check this using the command
    ```
    df -h
    ```
    Since the deployment of Galaxy will increase the size of the galaxy folder, it is
    safer to move this folder on a larger volume that is mounted at /mnt
    To do so, just type:
    ```
    mv /root/galaxy /mnt/mydatalocal/ && cd /mnt/mydatalocal/galaxy
    ```
  
  Then type `sh run.sh` and press the `enter` key !
  
  You should see an abundant log scrolling down. Don't worry !

  - All Galaxy dependencies required for the Galaxy server instance are being downloaded and installed
  - The Galaxy computing environment is automatically set up
  - the Galaxy web server is installed and static pages are built (this step specifically takes more and more time)
  - The Galaxy database (sqlight) is automatically upgraded to its latest structure/model
  - The package manager Conda, which is heavily used by Galaxy to install its tools is installed.
  
  After 5-10 minutes, you should see the log stopping with:

```
Starting server in PID 3813.
serving on http://127.0.0.1:80
```
### 4. Connect to your living Galaxy instance

You should now be able to access to you Galaxy instance in a your web browser window.

  - Go back to your Google Cloud Engine control panel.
  - Find the `External IP address` / `Adresse IP externe` in the 7th column of the dashboard
  (to the left of the ssh menu that you used before).
  - Click on the hyperlink.

!!! tip "In the new browser window, connect as an admin to your Galaxy server instance"
    - Register to your instance using the email address you put in the galaxy.yml at step 7
    (menu "Authentification et enregistrement --> Enregistrement)
    - Now that you are registered, you can log in using the same login and password you
    have chosen.
    - After login, you will see the admin tab in the top menu of the Galaxy interface.
    
    ==You are connected to Galaxy as an admin !==
