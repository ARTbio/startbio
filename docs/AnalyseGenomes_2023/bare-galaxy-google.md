### 0. Prerequisite

- [x] You have obtained and activated your Google Coupon for this training as described in
[Appendix 1](../Google_cloud_Account)
- [x] You have accessed to the Google dashboard and tested Starting and Stopping a virtual
machine (VM) instance as described in [Appendix 1](../Google_cloud_Account)

### 1. Spin off a virtual Machine `bare-galaxy` with ![](images/google-padok.png){width="30px" align="bottom"} Google Cloud Engine

- [x] Connect to your Google Compute Instances
  [dashboard](https://console.cloud.google.com/compute/instances)

- [x] Create a Virtual Machine Instance
 

!!! info "with the following settings"
    - Name: `bare-galaxy`
    - Region `europe-west6 (Zurich)` (or any region available with you Google coupon). ==As
    it is unlikely that a single Google zone will be able to provide enough resources
    to support 18 virtual machines at the same time, we will have to
    distribute our instances to different zones in Europe and USA==.
    - Zone: `europe-west6-a` (or `-b` or `-c`)
    - **Configuration de la machine**
        - `USAGE général`
        - Série: `E2`
        - Type de machine: `PRÉDEFINI` :arrow_forward: `Standard` :arrow_forward: `e2-standard-8`
    - **Disque de démarrage (Modifier)**
        - `IMAGES PUBLIQUES`
        - Système d'exploitation: `Ubuntu`
        - Version*: `Ubuntu 20.04 LTS`
        - Type de disque de démarrage: `Disque persistant avec équilibrage`
        - Taille (Go): ==`200`==
        - ==SELECTIONNER==
    - **Pare-feu**
        - Check `Autoriser le trafic HTTP`

These settings should look like:
    
![](images/GCE_spin.png){: style="width:600px"}
![](images/GCE_OS.png){: style="width:450px"}
![](images/GCE_firewall.png){: style="width:450px"}

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