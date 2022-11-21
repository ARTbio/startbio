??? info "The `screen` program"
    `screen` is a useful linux program that creates, attach, detach or
    reattach "virtual" shell sessions.
    `screen` allows to run simultaneous linux processes in isolated environments which can
    be put in the background while working with the console at other tasks
    
    A handful of screen commands you should know:
    
    - `screen -ls` lists all screen sessions available, attached (currently active) or
      detached in the background
    - `screen -r <session>` reattaches a detached screen session.
    - _within an active screen session_ ++ctrl++++a++ then ++d++ detaches the active session
    - _within an active screen session_ type `exit` to terminate the active session
    - `screen -S <session>` creates a new session
    - `screen -d -r <session>` detaches the current session and reattach another one


!!! bug "See the Galaxy server logs"
    When some tools are failing, you may grab useful information.
    ```
    tail -f ./galaxy/database/gravity/log/gunicorn.log
    ```

!!! bug "If tools fail with libssl / openssh issue in the bug report"
    ```
    cd /lib/x86_64-linux-gnu
    ln -s libssl.so.1.1 libssl.so.1.0.0
    ln -s libcrypto.so.1.1 libcrypto.so.1.0.0
    ```
    We are not fan of this, it is a rather dirty turnaround

!!! bug "For Conda Geek only"
    In case of problems with some conda packages, you may try command like:
    ```
    conda install -c bioconda samtools=1.9 --force-reinstall
    conda install -c bioconda --force-reinstall ucsc-genepredtobed ucsc-gtftogenepred
    ```
    after activating the proper conda environment 
    
!!! bug "To shrink you `_conda` dependencies folder"
    This folder is located at `/home/galaxy/tool_dependencies/_conda` and you must first
    activate the Galaxy conda base environment using the command (from anywhere):
    ```
    source /root/galaxy/database/dependencies/_conda/bin/activate
    ```
    Then
    ```
    conda clean --dry-run --all
    ```

:warning: the following content was not updated in 2022, do not use it, it is deprecated
!!! bug "In case of many issues with many tools ==which are not experienced by most of other students=="
    There is a problem with your conda tool dependencies.
    
    You can use the :dizzy: ==**magic patch**== :dizzy:
    
    Copy and paste the code below in the (root) terminal. This will swap your _conda environment
    with a clean, fresh one.
    
    ```
    cd ~ && \
    wget https://raw.githubusercontent.com/ARTbio/Run-Galaxy/master/deployment_scripts/patch_conda_env.sh && \
    sh patch_conda_env.sh
    ```
