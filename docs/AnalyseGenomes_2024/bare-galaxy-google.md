### 0. Prerequisite

- [x] You have obtained and activated your Google Coupon for this training as described in
[Appendix 1](../Google_cloud_Account)
- [x] You have accessed to the Google dashboard and tested Starting and Stopping a virtual
machine (VM) instance as described in [Appendix 1](../Google_cloud_Account)

### 1. Spin off a virtual Machine `bare-galaxy` with ![](images/google-padok.png){width="30px" align="bottom"} Google Cloud Engine

Before starting, we recommend you to pay extra attention any time you see the :warning:
signal.

- [x] Connect to your Google Compute Instances
  [dashboard](https://console.cloud.google.com/compute/instances)

- [x] Create a Virtual Machine Instance

!!! info "with the following settings"
    - Name: `ansible-galaxy`
    - Region `europe-west6 (Zurich)` :point_left: Check your Region in the popup table bellow 
    - Zone: `europe-west6-a` (or `-b` or `-c`) :point_left: Check your Zone in the popup table bellow
    - **Configuration de la machine**
        - `USAGE général`
        - Série: `E2`
        - Type de machine: `PRÉDEFINI` :arrow_forward: `Standard` :arrow_forward: `e2-standard-8`
    - **Disque de démarrage (Modifier)**
        - `IMAGES PUBLIQUES`
        - Système d'exploitation: `Ubuntu`
        - Version: `Ubuntu 20.04 LTS` :point_left: :warning: Watch to the version number (20.04).
          When several processor types are available (eg x86, amd, ...) you can choose anyone.
        - Type de disque de démarrage: `Disque persistant avec équilibrage`
        - Taille (Go): ==`200`==
        - ==SELECTIONNER==
    - **Pare-feu**
        - Check `Autoriser le trafic HTTP`

??? warning "Region and Zone assignments to students :warning:"
    ==As it is possible that a single Google zone is not able to provide enough resources
    to support 18 virtual machines at the same time, we will distribute our instances to
    different zones in Europe==.
    
    The following table assigns the instances by name to different Regions and Zones.

    Please respect this attribution for your final instance, the one you will use during
    your practical work.
    
    |    Email prefix    |         Region           |      Zone         |
    |--------------------|--------------------------|-------------------|
    |alleon.gaelle       |europe-west1 (Belgique)   |europe-west1-b     |
    |enzo.becherel       |europe-west1 (Belgique)   |europe-west1-d     |
    |emma.benbakir       |europe-west1 (Belgique)   |europe-west1-c     |
    |samuel.bensoussan   |europe-west2 (Londres)    |europe-west2-c     |
    |tberthom            |europe-west2 (Londres)    |europe-west2-b     |
    |gregblavier76       |europe-west2 (Londres)    |europe-west2-a     |
    |faroukbouraima      |europe-west3 (Francfort)  |europe-west3-c     |
    |lunadebarbarin      |europe-west3 (Francfort)  |europe-west3-a     |
    |baptiste.demaret    |europe-west3 (Francfort)  |europe-west3-b     |
    |nicolas.doucet      |europe-west6 (Zurich)     |europe-west6-a     |
    |maeva.drai          |europe-west6 (Zurich)     |europe-west6-b     |
    |yoann.gonneau       |europe-west6 (Zurich)     |europe-west6-c     |
    |sarah.graine        |europe-west9 (Paris)      |europe-west9-a     |
    |margot.hully        |europe-west9 (Paris)      |europe-west9-b     |
    |nathan.lacombe      |europe-west9 (Paris)      |europe-west9-c     |
    |jules.richez.22     |europe-west10 (Berlin)    |europe-west10-a    |
    |loann.paterour      |europe-west10 (Berlin)    |europe-west10-b    |
    |mathilde.quibeuf    |europe-west10 (Berlin)    |europe-west10-c    |
    |michiel.tawdarous   |europe-southwest1 (Madrid)|europe-southwest1-a|
    |oceane.wauthier     |europe-southwest1 (Madrid)|europe-southwest1-b|

These settings should be similar to this:
    
![](images/GCE_spin.png){width="600px"}
![](images/GCE_OS.png){width="450px"}
![](images/GCE_firewall.png){width="450px"}

**When**

- [x] you have double-checked all indicated settings
- [x] you are sure that your instance will start in the zone assigned to you

**Then** you can start you instance by clicking the button

![](images/creer_instance.png){width="350px"}


??? bug "Trouble shouting"
    **In some occasions, launching of your VM may fail** as illustrated bellow:
    ![](images/instance_failing.png){: style="width:600px"}
    
    1. Maybe you are not, indeed, using the billing account associated to your
    Google coupon, but instead using a billing account associated to a "Free Trial".
        
        - [x] If it is not already done, activate your coupon by following the received
        instructions, and be sure that you activate a project associated with the billing
        account of the coupon.
    2. The Region and Zone which you have chosen (in the example, `europe-west6-a`) is
    overloaded.
        
        - [x] In this case, try another `Zone` (-b or -c), and/or another `Region`, in Europe or
        America.

### 2. Connect to the VM using the ssh web console

Roll down the `ssh` menu in the control pannel and select the first option

`Ouvrir dans une fenêtre du navigateur`

![Select ssh session in browser](images/select_ssh.png)
    
**This opens a web ssh shell session to control your VM:**

![](images/web_ssh_console.png)

---