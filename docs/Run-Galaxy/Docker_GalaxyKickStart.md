## Installation of a Galaxy server with Docker

### What is Docker ?

![Docker](images/docker.png)

#### Virtual machines

Virtual machines (VMs) are an abstraction of physical hardware turning one server into
many servers. The hypervisor allows multiple VMs to run on a single machine.
Each VM includes a full copy of an operating system, one or more apps, necessary binaries
and libraries - taking up tens of GBs. VMs can also be slow to boot.

#### Containers

Containers are an abstraction at the app layer that packages code and dependencies
together. Multiple containers can run on the same machine and share the OS kernel with
other containers, each running as isolated processes in user space. Containers take up
less space than VMs (container images are typically tens of MBs in size), and start
almost instantly.

###  GalaxyKickStart Docker **Container**

Instead of using the GalaxyKickStart playbook in a VM, the playbook can be used to build
a Docker container image that will be an almost exact mirror of the GalaxyKickStart VM
you have just built.

You are not going to do that today (although you should be able to do it by reading the instructions).

Instead, you are going to

- Install the `docker` system
- pull the GalaxyKickStart docker container that is deposited in the [`Docker Hub`](https://hub.docker.com/u/artbio/)
- run this docker container and connect to the deployed GalaxyKickStart server instance

### Deployment

- There is actually no need for a new VM, the ansible already installed the docker service
in the VM used to deploy GalaxyKickStart.
- If not already, be connected to you VM as root user using the Google ssh console (`sudo -i`)
- download the script `run_docker_analyse_genomes_2019.sh` using the command
`wget https://raw.githubusercontent.com/ARTbio/Run-Galaxy/master/deployment_scripts/run_docker_analyse_genomes_2019.sh`
- run the script using the command `sh run_docker_analyse_genomes_2019.sh`
- Connect to your docker-deployed "GalaxyKickStart" instance:
    
    Just click on the url displayed in your Google Cloud Engine Console
    and connect using the login:password `admin@galaxy.org:admin`
    
    
#### The docker_galaxy.sh script explained

NB: in the following code, numbers in line heads should be removed to run the script.

```
 1 #!/usr/bin/env sh
 2 set -e
 3 echo "Now pulling the galaxykickstart docker image from DockerHub\n"
 4 supervisorctl stop all
 5 docker pull $1
 6 echo "Running galaxykickstart docker container\n"
 7 export DOCKER_INSTANCE=`docker run -d -p 80:80 -p 21:21 -p 8800:8800 \
 8           --privileged=true \
 9           -e GALAXY_CONFIG_ALLOW_USER_DATASET_PURGE=True \
10           -e GALAXY_CONFIG_ALLOW_LIBRARY_PATH_PASTE=True \
11           -e GALAXY_CONFIG_ENABLE_USER_DELETION=True \
12           -e GALAXY_CONFIG_ENABLE_BETA_WORKFLOW_MODULES=True \
13           -v /tmp/:/tmp/ \
14           -v /export/:/export \
15           $1`
16 echo "The $1 docker container is deploying...\n"
17 sleep 90
18 echo "The $1 docker container is up and running\n"
19 docker logs  $DOCKER_INSTANCE
20 docker exec $DOCKER_INSTANCE sudo su galaxy -c '/home/galaxy/galaxy/.venv/bin/pip install cryptography==2.2.2'
21 docker exec $DOCKER_INSTANCE sudo su galaxy -c 'cd ~/galaxy/config && wget https://raw.githubusercontent.com/ARTbio/Run-Galaxy/master/deployment_scripts/sanitize_whitelist.txt'
22 echo "Galaxy in container will restart to take into account new settings\n"
23 sleep 120
24 docker exec $DOCKER_INSTANCE sudo supervisorctl restart galaxy:
```

1. The shebang line. Says that it is a script code and that the interpreter to execute the
code is sh and can be found in the /usr/bin/env environment
2. set -e says to the sh interpreter to exit the run at first error (to avoid catastrophes)
3. prompts "Now pulling the galaxykickstart docker image from DockerHub"
4. stop the galaxy services (galaxy, postgresql, slurm, nginx,...) that were deployed before with ansible (to liberate ports)
5. Pulls (Downloads) the Docker Image specified as parameter to the script (either`artbio/biologiegenome`
or `artbio/rna-biologie-genome`) from the [DockerHub repository](https://hub.docker.com/r/artbio/biologiegenome/)
6. reports this action to the terminal
7. Lines 7 to 15 are actually a single command to run an instance of the galaxykickstart
docker image. Note the `\` at ends of lines 7 to 14: this character `\` specify that the
code line is continued without line break for the bash interpreter.

    The line 7 starts with an `export DOCKER_INSTANCE=` instruction. This means that the result
    of the command between \` after the sign `=` will be put in the environmental variable
    `DOCKER_INSTANCE`, available system-wide.

    Now, the docker command (between \`) itself:

    Still in line 7, we have `docker run -d -p 80:80 -p 21:21 -p 8800:8800`.

    This means that a container will be run as a deamon (`-d ` option) and that the internal
    TCP/IP ports 80 (web interface) and 21 (ftp interface) of the docker instance will be mapped
    to the ports 80 and 21 of your machin (The VM in this case). Note that in the syntax `-p 80:80`,
    the host port is specified to the left of the `:` and the docker port is specified to the right
    of the `:`.
    
8. docker command continued: here we specify that the docker container acquires the root privileges
9. docker command continued: `-e GALAXY_CONFIG_ALLOW_USER_DATASET_PURGE=True`.

    The -e option specifies an environmental variable `GALAXY_CONFIG_ALLOW_USER_DATASET_PURGE`
    passed (`e`xported) to the docker container with the value `True`
    `galaxy_manage_trackster: true` with the string `galaxy_manage_trackster: false`
    in the ansible configuration file `groups/all`.

10. The environmental variable `GALAXY_CONFIG_ALLOW_LIBRARY_PATH_PASTE` is exported to
the docker container with the value `True`
11. The environmental variable `GALAXY_CONFIG_ENABLE_USER_DELETION` is exported to
the docker container with the value `True`
12. The environmental variable `GALAXY_CONFIG_ENABLE_BETA_WORKFLOW_MODULES` is exported to
the docker container with the value `True`

    Note that all these exports in the docker command correspond to advanced boiling/tuning of the docker container.
    You are not obliged to understand the details to get the container properly running.
13. Now the -v is important, better to understand it !

    -v stands for "volume". the `-v` option says to export the /tmp directory of the docker container
    to the /tmp directory of the host.
14. we also export the /export directory of the container (any docker container has or
should have by default an /export directory) to an /export directory of the host (your VM here).

    Note that if the /export directory does not exists at docker run runtime, it will be created.

    So it is important to understand the -v magics: every directory specified by the -v option will be shared
    between the docker container filesystem and the host filesystem. It is a mapping operation, so that
    the same directory is accessible either from inside the docker container or from inside the host.

    Now, if you stop and remove the docker container, all exported directory will persist in the host.
    If you don't do that, all operations performed with a container are lost when you stop this container !

15. This is the end of the docker run command. The docker image to be instantiated is specified by $1 variable,
the parameter passed to the script at runtime.
16. reports to the terminal user
17. wait 90 sec during the docker container deployment
18. reports to the terminal user
19. Now that the docker container is launched, you can access its logs with the command
`docker logs ` followed by the identification number of the docker container.
We have put this ID in the variable `DOCKER_INSTANCE`
and we access to the content of this variable by prefixing the variable with a `$ `:
`docker logs $DOCKER_INSTANCE`
20. This is just an update of the python package `cryptography` required for the galaxy server. You see
here an example of how a user can interact with the docker container to adjust the services it
is providing.
21. another update of a configuration file for galaxy (job_conf.xml)
24. Restart Galaxy inside the container to take the setting changes into consideration


