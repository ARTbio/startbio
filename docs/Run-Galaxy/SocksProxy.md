## How to access to a machine deployed in the IFB core cloud ?

### The problem

All virtual machines deployed in the IFB core are located in a subnetwork whose access is
limited to

- the port 22, for ssh connections
- the port 443, for https (web) connections.

Note that the port 80 _is not open_, precluding connection through the "insecure" http port of your
web Galaxy server.

Accessing a web server running on a virtual instance through https (443) requires
that each machine has declared its own SSL certificate and most preferably owns a unique
domain name, in the form of `mymachine.ifb.fr`.
Although there are turnarounds for generating self-signed SSL certificate for cloud instances,
this implies manipulations which are far beyond the scope on this training for beginners.

There is a least 2 ways for circumventing the https limitation for this training.

### 1. Running a SOCKS proxy on your local machine

  A. using your remote virtual machine
    
  Type in a terminal session  (and leave it alive):
    
  ```bash
  ssh -A -D 9900 ubuntu@134.158.247.85
  ```
  **OR**
  
  B. using another virtual machine of the IFB cloud
    
  Type in a terminal session (and leave it alive):
    
  ```bash
  ssh -i .ssh/ifbsocks -D 9900 ubuntu@<communicated.ifb.ip.address> # with the ifbsocks private key which you will be given
  ```

**THEN**

- open your network settings
- go to your proxy settings
- Check the box for SOCKS Proxy (v4 or v5)
- in the field for the Server Proxy SOCKS address, enter `localhost`
- in the field for the Server Proxy SOCKS port, enter `9900`

**From this point**

You should be able to access directly to your cloud Galaxy server by typing 

`http://<IFB.IP.your.server>`

This IP address is available at the [IFB biosphere interface](https://biosphere.france-bioinformatique.fr/cloud/deployment/)

![Docker](images/docker.png)


### 2. Tunneling of the unaccessible port 80 through an accessible ssh (22) port


