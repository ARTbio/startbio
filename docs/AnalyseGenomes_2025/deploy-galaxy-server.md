### 1. Connect your VM as Root

If you have previously suspended your VM, reactivate it and open a WEB SSH window.
If you just completed the previous section, type :

```
sudo -i
```

??? info "What does `sudo -i` command ?"
    This command open a new `shell` where you are root. You can check this by typing `whoami`
    that should return `root`, meaning that you are now working as `root` user.
    
    This is required because installation of new programs as well as manipulations of network
    interfaces is permitted only to users with administration rights.

### 2. Installation of the Galaxy server
??? warning "Recommendations before starting"
   
    The creation of your Galaxy server includes the setup of the Galaxy Services and the
    installations of ~28 bioinformatics tools to analyse sequencing datasets.
    
    Although it is completely scripted and requires minimal intervention from your part,
    this process **takes ==~1 hour in total==, once,** and the deployed server will serve you
    for ==_the rest of the training week_==.
    
    Therefore, we ask you **extra focus** on this section as well as ==preparing your
    Galaxy server in advance of the Galaxy training week== **starting on Monday
    2nd, December - 2024**.
    
    **Practical recommendation about internet connection:**
    
    The deployment of the Galaxy server and the installation of Galaxy tools
    in the server involves **remote execution** of scripts in your Virtual Machine.
    Therefore, it is much better that the internet connection between your local
    terminal (where you are physically working) and the remote VM ==**STAYS UP**==
    and that you stay physically around your screen during these two phases.
    
    Some local machines are configured to sleep after a certain amount of time of inactivity.
    This sleeping process may stop your connection with the VM and consequently
    give you the inconfortable feeling that you "lost the thread".
    However, should you loose your SSH connection,
    for any reason, do not [panic and restart everything from scratch.]
    
    Indeed, once the installation software is triggered, it will continue in the
    background of ubuntu until its completion and independently from the SSH
    connection. 
    
    --> Just reconnect to your VM instance and follow the running installation
    whose log is keept in /root/install_log.txt, using the command
    `tail -f /root/install_log.txt`

We have automated the installation of Galaxy on your Google Virtual Machine.
All you need is to clone a `galaxyXpand` folder in your VM and run a bash script,
using a single command.

--> Copy the full content of the box below and paste it in your ssh terminal.

```
git clone https://github.com/artbio/galaxyXpand -b ubuntu24.04 && \
screen -d -m bash ~/galaxyXpand/scripts/deploy_ag2025.sh && \
sleep 5 && tail -f ~/install_log.txt
```
??? info "What is `git` command doing ?"
    This command is cloning the GitHub repository @artbio/galaxyXpand into a
    local folder named `galaxyXpand`.
    
    galaxyXpand is a software developped to quickly and easily install a Galaxy
    server. It is based upon the ansible framework for software deployment.
??? info "What is `screen -d -m` doing ? (:metal: Linux geek corner)"
    `screen -d -m <command>` is starting the `<command>` in a separate child shell
    and a "detached" mode. In this particular case, `deploy_ag2025.sh` is run in a
    "orphean child shell". This way, interruption of your ssh connection with the
    parent shell will not interrupt the detached shell process.
    
    You can see it as a small "daemon" program :smiling_imp:.

??? info "What is `sh ~/galaxyXpand/scripts/deploy_ag2025.sh` doing ?"
    This command runs the script
    [deploy_ag2025.sh](https://github.com/ARTbio/galaxyXpand/blob/ubuntu24.04/scripts/deploy_ag2025.sh)


### 3. Monitoring the deployment of the Galaxy server

The tasks executed by the `deploy_ag2024.sh` are displayed in your terminal
(thanks to the `tail -f ~/install_log.txt` command) as well as logged in the
file `install_log.txt`.

:warning: Although the installation log in your terminal may seem to stop for
several minutes (because of long internal steps), it is only when the following
lines show up that the Galaxy Installation is finished. Keep cool :sunglasses:
and remember that the full deployment lasts for ~50 mins.

??? info "Last lines of install_log.txt"
    
    ```bash
    changed: [localhost] => (item={'name': 'sambamba', 'owner': 'artbio', 'tool_panel_section_id': 'samtools', 'tool_panel_section_label': 'Samtools', 'tool_shed_url': 'https://toolshed.g2.bx.psu.edu/'})
    changed: [localhost] => (item={'name': 'bedtools', 'owner': 'iuc', 'tool_panel_section_id': 'bedtools', 'tool_panel_section_label': 'Bedtools', 'tool_shed_url': 'https://toolshed.g2.bx.psu.edu/'})
    
    TASK [install.galaxy-tools : include_tasks] ************************************
    kipping: [localhost]
    
    TASK [install.galaxy-tools : include_tasks] ************************************
    skipping: [localhost]
    
    PLAY RECAP *********************************************************************
    localhost                  : ok=11   changed=5    unreachable=0    failed=0    skipped=3    rescued=0    ignored=0   
    
    Wed Nov 13 18:05:40 UTC 2024
    -- Installation is complete --
    ```

??? info "The main steps of the Galaxy server deployment"
    - The Ubuntu system is updated at its latest version
    - Python dependencies required for the Galaxy server instance are downloaded and installed
    - The ansible framework v3.0.0 is installed for running the ansible playbooks
    - The Galaxy computing environment (virtualenv) is automatically set up
    - The Galaxy web server is installed (nginx reverse proxying gunicorn) and static pages are built
    - The Galaxy database Postgresql is installed and upgraded to its latest structure/model
    - The package manager Conda, which is heavily used by Galaxy to install its tools, is installed.
    - Plus many other tasks : a high-performance server relies on complex software.
    - [x] The final step in Galaxy deployment is the automated installation of
      28 tools which you will need for your analyses.
    
    **Naturally, this deployment will happen once. The next time you connect to your
    Galaxy server, you'll be ready to use it !**

### 4. :point_right: Galaxy Server Check

- [x] Type ++ctrl++++c++ to get the hand back over your web terminal
- [x] Copy and paste the last lines of the installation log in
  a new post in this
  [GitHub Discussion](https://github.com/ARTbio/AnalyseGenome/discussions/40){:target="_blank"}.
  These lines should be similar to the textbox above ("Last lines of install_log.txt")
- [x] Enter the following command line :
  ```
  galaxyctl status
  ```
- [x] copy the returned output (:warning: *copy* is not *screenshot*) and paste it in
  the same post in the [GitHub Discussion](https://github.com/ARTbio/AnalyseGenome/discussions/40){:target="_blank"}
<center>
![](images/checkpoint.png){width="80"}
</center> 

### 5. Connect to your living Galaxy instance

You should now be able to access to you Galaxy instance in a web browser window.

- Go back to your Google Cloud Engine control panel.
- Find the `External IP address` / `Adresse IP externe` in the 7th column of the dashboard
  (to the left of the ssh menu that you used before).
  
  ![externIP](images/externIP.png)
  
- Click on the hyperlink.
- In the new browser window, follow the menu `Authentification et enregistrement`
  --> `Enregistrement` --> `Don't have an account? Register here.`
  
  ![register](images/register.png){ width="300" }

  and  **register** to your instance using the email address
  ```
  admin@galaxy.org
  ```
  and the password of your choice (:warning: don't forget it). Use an easy to remember one as you will need to type it frequently.
  
- After login, you should see the admin tab in the top menu of the Galaxy interface.
  
  ![](images/admin_menu.png){ width="600" }
  
  ==You are connected to Galaxy as an admin !==

??? bug "If you do not see the `admin` menu"
    This is most likely due to an error in the login which must be
    exactly `admin@galaxy.org`
    
    --> Just register another user with the proper email address
---
