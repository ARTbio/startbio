## Install a minimal galaxy server with git

### Spin off a virtual Machine `bare-galaxy`
You may have already done this in the [previous section](spin_off_VM.md). If not, refer to this section
We are going to use a GCE VM 
- under `Ubuntu 16.04 LTS`
- `2` processors
- `7.5` Gb RAM
- a `50 Gb` Volume (more than enough)

### Connect to the VM as explained in [this section](spin_off_VM.md) using the ssh web console

### Installation of dependencies

To install Galaxy, we only need a few dependencies (i.e. pre-installed programs which Galaxy needs) and a limited number of command line.
We are going to execute these instruction as the `root` unix user. This is easier because installation
of new programs as well as manipulations of network interfaces is generally permitted only
to users with administration rights.

So let's do this step by step:

##### 1. `sudo -i`

This command open a new "shell" where you are root. You can check this by typing `pwd` that
should return `/root/`, meaning that you are now working in the directory of the `root` user.

#### 2. `apt update -y`
This command will just trigger an update of the program database in the Ubuntu OS.

##### 3. `apt install -y python-dev python-pip nano git` 

This command install some python programs (python-dev python-pip) that are intensively used by Galaxy,
`nano`, a simple text editor we need, and the `git` program which is the software to fetch,
and update Galaxy (i.e. a sort of "installer" program). The `-y` option specifies to the `apt-get`
package installer that no confirmation is needed for this command.

##### 4. `git clone https://github.com/galaxyproject/galaxy.git -b release_18.09`

This command says to use `git` to `clone` the code repository located at `https://github.com/galaxyproject/galaxy.git`.
In addition the `-b release_18.09` option specifies that only the version `release_18.09` will be cloned locally in your virtual machine.
You may try to visualize the URL [https://github.com/galaxyproject/galaxy.git](https://github.com/galaxyproject/galaxy.git)
in your web browser. You will literally see the code of Galaxy !

##### 5. `cd galaxy`

This command shift you in the `galaxy` directory that was created by git

##### 6. `cp config/galaxy.yml.sample config/galaxy.yml`

This command makes a copie of the `galaxy.yml.sample` file into `galaxy.yml` - in the
directory `config` that is in the `galaxy` directory.

##### 7. `nano config/galaxy.yml`

Using this command, we are going to edit some important settings that are required to run our Galaxy fresh instance.

- Find the line `  http: 127.0.0.1:8080` and edit it to `  http: 0.0.0.0:80`

By doing this, we ensure that we will be able to reach the galaxy web server on our virtual machine using the usual web port `80`.

- Find the line `#admin_users: ''`, delete the `#` character and type your email address between the two single quotes.

Any email address is ok (me@myname.fr is you want). It is just used here as an admin identifier.

- save your edits by pressing the key combination `Ctrl`+`o`
- quit nano by pressing the key combination `Ctrl`+`x`

##### 8. Ready for the Big Bang ?

Then type `sh run.sh` and press the `enter` key !

You should see an abundant log scrolling down. Don't worry !
- All Galaxy dependencies required for the Galaxy server instance are being downloaded and installed
- The Galaxy computing environment is automatically set up
- the Galaxy web server is installed
- The Galaxy database is automatically upgraded to the latest structure
- Various tools are upgraded.

After 8 minute or so, you should see the log freezing with

```
Starting server in PID 3813.
serving on http://127.0.0.1:80
```
##### 8. Connect to your living Galaxy instance

If so, this is all good, and you can now access to you Galaxy instance in a you web browser window:

Go back to your Google Cloud Engine control panel. Find the `External IP address` / `Adresse IP externe`
in the 7th column of the Dashbord (to the left of the ssh menu that you used before. And just click on the hyperlink.

##### 9. Connect as an admin of your Galaxy server instance

- Register to your instance using the email address you put in the galaxy.yml at step 7 (menu "Authentification et enregistrement --> Enregistrement)
- Now that you are registered, you can log in using the same login and password you have chosen.
- After login, you will see the admin tab in the top menu of the Galaxy interface.

###### You are connected to Galaxy as an admin !








