All commands below should be run as the admin user (`sudo -i`)

!!! bug "see your slurm jobs triggered by Galaxy"
    ```
    watch "squeue --format '%.18i %.9P %.40j %.8u %.2t %.10M %.4C %.6D %R' && sinfo"
    ```
!!! bug "If your slurm cluster is stuck (the datasets stay in the "grey state" forever)"

    - the previous command should show it
    - then try this:
    
    ```
    scontrol update NodeName=<name of your instance> State=Resume
    ```

!!! bug "See the Galaxy server logs"
    When some tools are failing, you may grab useful information.
    ```
    tail -f /home/galaxy/galaxy/uwsgi.log
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
!!! bug "To shrink you `_conda` dependencies folder"
    This folder is located at `/home/galaxy/tool_dependencies/_conda` and you must first
    activate the Galaxy conda base environment using the command (from anywhere):
    ```
    source /home/galaxy/tool_dependencies/_conda/bin/activate
    ```
    Then
    ```
    conda clean --dry-run --all
    ```

!!! bug "In case of many issues with many tools"
    There is a problem with your conda tool dependencies.
    
    You can use the :dizzy: ==**magic patch**== :dizzy:
    
    Copy and paste the code below in the (root) terminal. This will swap your _conda environment
    with a clean, fresh one.
    
    ```
    cd ~ && \
    wget https://raw.githubusercontent.com/ARTbio/Run-Galaxy/master/deployment_scripts/patch_conda_env.sh && \
    sh patch_conda_env.sh
    ```
