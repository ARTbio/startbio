!!! info "Galaxy server commands, quick reminder"
    Check that the server is running and see last logs
    
    :warning: The following commands won't working if you are not logged as root.
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