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


  You should see an abundant log scrolling down. Don't worry !

  - All Galaxy dependencies required for the Galaxy server instance are being downloaded and installed
  - The Galaxy computing environment is automatically set up
  - the Galaxy web server is installed and static pages are built
  - The Galaxy database (sqlight) is automatically upgraded to its latest structure/model
  - The package manager Conda, which is heavily used by Galaxy to install its tools is installed.
  
  After ~10 minutes, you should see in the ssh console the following log:

```{.bash title="Terminal"}
Galaxy is now running as a daemon in the background
```
The log of the Galaxy server can be displayed using the command

```
Galaxy server instance 'main.1' is running
serving on http://0.0.0.0:80

```

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
