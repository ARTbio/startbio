## Go back to your Google ssh terminal

- Run the following script using the bash interpreter:
```
source /root/.bashrc && \
bash /root/AnalyseGenome/GalaxyServer/install_galaxy_tools.sh
```
You should now be asked for an API key:
```
please enter your admin API key: 
```
The next section explains how to generate and copy this API key to interact programmatically
with the Galaxy server and be able to install tools

## Go back to your Galaxy web window

In the main menu `User` --> `Preferences`

![user preferences](images/user_preferences.png){width="200"}

Select `Manage API Key`, click `Create a new Key`, and copy the current API key

## Paste the copied API key in the Google ssh terminal

and press the ++return++ key

The tool installation should start immediately and last for about 15 minutes (you should
see tools installing one at a time)

---