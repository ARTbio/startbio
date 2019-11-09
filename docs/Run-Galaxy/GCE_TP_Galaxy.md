### Spin off a virtual Machine
##### 1. Go to the Google Cloud Dashboard and select "Compute Engine" on the left hand menu bar
##### 2. Select the submenu "Instances de VM"

![Instances](images/IntancesVM.png)

##### 3. Click on the top bar menu the "CREER UNE INSTANCE" panel

![create instance](images/create_instance.png)
    
##### 4. Put name `my-galaxy-server`, Zone `europe-west1-b (or c)`, Type de machine `8 vCPU` + `30 Go` de mémoire.
![zone, cpu and ram](images/instance_naming.png)
##### 5. Disque de Démarrage: Click on `Modifier`
![click on modify startup image](images/modify_startup_image.png)
##### 6. Select the top menu `images personnalisées` (`custom images`)
![click the images personnalisées menu](images/startup_disk_panel.png)

##### 7. Click on the rolling menu `Afficher les images de` and select the `My Project - main-sunset-133416`
![select My Project](images/select_my_project.png)
    
What is important here is the identifier `main-sunset-133416`

##### 8. Check the button to select `galaxy-image-pasteur`
![Check Galaxy Image](images/check_galaxy_image.png)

##### 9. At the bottom of the same form, choose `100 Go` for the Disk Size (Taille). Note that this size should be already selected.
![choose size disk and select out](images/choose_galaxy_disk_size.png)
    
Click the `Sélectionner` button to leave the selection `Disque persistant standard` / `Standard persistant drive`
    
##### 10. Back to the main form, Click `Authorize HTTP traffic` / `Autoriser le traffic HTTP`

![Select HTTP traffic](images/http.png)

##### 11. Click `Créer` / `Create`
##### 12. After ~1 minute or so, the VM turns on "green" and an `ssh` menu becomes selectable

![Instance fire up](images/instance_turns_green.png)
    
    
##### 13. Click on the http link provided in the `Adress IP externe` column
You should now be able to access to your own Galaxy server instance, but not that this
phase can take an additional minute or so, this is the time to start all the galaxy services
in the new server instance.
    
##### 14. Immediately Log in to your server as the administrator

![first login](images/Galaxy_first_login.png)
   	
And log in with `admin@galaxy.org` : `admin`
   	
![authenticate](images/Admin_authentication.png)
   		
## YOU ARE READY TO USE GALAXY !
