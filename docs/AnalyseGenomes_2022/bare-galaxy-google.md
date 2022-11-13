### 1. Spin off a virtual Machine `bare-galaxy` with ![](images/google-padok.png){: style="width:30px"} Google Cloud Engine

- Connect to your Google Compute Instances
  [dashboard](https://console.cloud.google.com/compute/instances)

- Create a Virtual Machine Instance
 

!!! info "with the following settings"
    - Name: `bare-galaxy`
    - Region `europe-west6 (Zurich)` (or any region available with you Google coupon). ==As
    it is very unlikely that a single Google zone will be able to provide enough resources
    to support 18 virtual machines at the same time, we will have to coordinate to
    distribute our instances to different zones in Europe and USA==.
    - Zone: `europe-west6-a` (or `-b` or `-c`)
    - **Configuration de la machine**
        - `OPTIMISEE POUR LE CALCUL` (or `COMPUTE-OPTIMISED`) ==:link:[in case of trouble](#trouble-shouting)==
        - Série: `C2`
        - Type de machine: `c2-standard-8 (8 processeurs virtuels, 32 Go de mémoire)`
    - **Disque de démarrage (Modifier)**
        - `IMAGES PUBLIQUES`
        - Système d'exploitation: `Ubuntu`
        - Version*: `Ubuntu 20.04 LTS`
        - Type de disque de démarrage: `Disque persistant avec équilibrage`
        - Taille (Go): ==`200`==
        - ==SELECTIONNER==
    - **Pare-feu**
        - Check `Autoriser le trafic HTTP`

This settings should look like:
    
![](images/GCE_spin.png){: style="width:450px"}
![](images/GCE_OS.png){: style="width:450px"}
![](images/GCE_firewall.png){: style="width:450px"}

??? bug "Trouble shouting"
    **In some occasions, launching of your VM may fail** as illustrated bellow:
    ![](images/instance_failing.png){: style="width:600px"}
    
    1. Maybe you are not, indeed, using the billing account associated to your
    Google coupon, but instead using a billing account associated to a "Free Trial".
        
        If this is the case, try either of the following solutions:
        
        - If it is not already done, activate your coupon by following the received
        instructions, and be sure that you activate a project associated with the billing
        account of the coupon.
        - Instead a selecting `OPTIMISEE POUR LE CALCUL` (or `COMPUTE-OPTIMISED`), select
        `USAGE GENERAL` (or `GENERAL-PURPOSE`) and scroll-down the **Machine-type** menu
        to select `e2-standard-8 (8 vCPU, 32 GB memory)`
    2. The Region and Zone which you have chosen (in the example, `europe-west6-a`) is
    overloaded.
        
        In this case, try another `Zone` (-b or -c), and/or another `Region`, in Europe or
        America.

### 2. Connect to the VM using the ssh web console

!!! info "ssh connection"
    Roll down the `ssh` menu in the control pannel and select the first option
    `Ouvrir dans une fenêtre du navigateur`

    ![Select ssh session in browser](images/select_ssh.png)
    
    **This opens a web ssh shell session to control your VM:**
    
    ![](images/web_ssh_console.png)


### 3. Installation of the Galaxy server

We have automated the installation of Galaxy on your Google Virtual Machine.
All you need is to (i) taking the control of the machine as root and (ii) downloading a 
bash script and running it.
So let's do this, step by step, using the ssh Terminal:

    
  ```Console
  sudo -i
  ```
!!! info ""
    This command open a new `shell` where you are root. You can check this by typing `whoami`
    that should return `root`, meaning that you are now working as `root` user.
    
    This is required because installation of new programs as well as manipulations of network
    interfaces is permitted only to users with administration rights.

```
wget https://raw.githubusercontent.com/ARTbio/AnalyseGenome/main/GalaxyServer/deploy_galaxy.sh
```
!!! info ""
    This command is downloading an installation script located in the GitHub repository
    @artbio/AnalyseGenome

```
sh deploy_galaxy.sh
```
!!! info ""
    Finally, this command runs the sh script
    [deploy_galaxy.sh](https://raw.githubusercontent.com/ARTbio/AnalyseGenome/main/GalaxyServer/deploy_galaxy.sh)
    
    Note that this occurs inside a `screen` session, so that you can detach from this session
    without killing the serving process.
    
    See bellow a few tips to navigate between `screen` sessions

As prompted by the screen, you can now follow the deployment of Galaxy by entering the 
command 
```
screen -r galaxyscreen
```

  You should see an abundant log scrolling down. Don't worry !

  - All Galaxy dependencies required for the Galaxy server instance are being downloaded and installed
  - The Galaxy computing environment is automatically set up
  - the Galaxy web server is installed and static pages are built
  - The Galaxy database (sqlight) is automatically upgraded to its latest structure/model
  - The package manager Conda, which is heavily used by Galaxy to install its tools is installed.
  
  After 5-10 minutes, you should see in galaxyscreen screen session the log:


```{.bash title="Terminal"}
Galaxy server instance 'main.1' is running
serving on http://0.0.0.0:80
```

??? info "The `screen` program"
    As we are going to launch Galaxy as a "daemon" server, we need a special session that we
    can leave or reattach without interrupting the running Galaxy server.
    
    Therefore, we are going to use the linux program `screen` which create, attach, detach or
    reattach a "virtual" sessions that end only when you decide to kill them.
    
    A handful of screen commands you should know:
    
    - `screen -ls` lists all screen sessions available, attached (currently active) or
      detached in the background
    - `screen -r <session>` reattaches a detached screen session.
    - _within an active screen session_ ++ctrl++++a++ then ++d++ detaches the active session
    - _within an active screen session_ type `exit` to terminate the active session
    - `screen -S <session>` creates a new session
    - `screen -d -r <session>` detaches the current session and reattach another one


### 4. Connect to your living Galaxy instance

You should now be able to access to you Galaxy instance in a your web browser window.

- Go back to your Google Cloud Engine control panel.
- Find the `External IP address` / `Adresse IP externe` in the 7th column of the dashboard
  (to the left of the ssh menu that you used before).
  
  ![externIP](images/externIP.png)
  
- Click on the hyperlink.
- In the new browser window, follow the menu `Authentification et enregistrement`
  --> `Enregistrement` --> `Register here`
  
  ![register](images/register.png){ width="300" }

  and  **register** to your instance using the email address
  `admin@galaxy.org` and the password of your choice (:warning: don't forget it)
  
- After login, you should see the admin tab in the top menu of the Galaxy interface.
  
  ![](images/admin_menu.png){ width="600" }
  
  ==You are connected to Galaxy as an admin !==
