You can control you galaxy server by connecting to a ssh shell session in your VM.
Just remember that you can easily open this type of session using the Google
Compute Engin dashboard, as you did before.

The preferred program to exert this control is `galaxyctl` which has to be
invoked with the root rights.

!!! info "The main options available with `galaxyctl`"
    :warning: The following commands won't work if you are not logged as root.
    (`sudo -i` at your initial connection)
    
    - [x] Check the state of the Galaxy server
    ```
    galaxyctl status
    ```
    - [x] Restart **gently** the Galaxy server
    ```
    galaxyctl graceful 
    ```
    - [x] Restart "firmly" the Galaxy server (if for whatever reason it is frozen)
    ```
    galaxyctl restart
    ```
    - [x] Follow the activity log of the Galaxy WEB server:
    ```
    galaxyctl follow
    ```
---