## Go back to your Google ssh terminal

- Detach from the current galaxyscreen session by using the key combination

++ctrl+a++ then ++d++

- Run the following script:
```
sh /root/AnalyseGenome/GalaxyServer/install_galaxy_tools.sh
```
You should now be asked for an API key, paste your API key and press ++enter++

## Go back to your Galaxy web window

In the main menu `User` --> `Preferences`

![user preferences](images/user_preferences.png){width="200"}

Select `Manage API Key`, click `Create a new Key`, and copy the current API key

## Paste the copied API key in the Google ssh terminal

and press the ++return++ key

The script will install a set of 21 Galaxy tools which are needed for the rest of your
training (in the `ephemeris` screen session).
This may take more than 15 min, be patient...

after pressing the ++return++ key, you can follow the tool installation by :

```
tail -f tools.log
```

Use ++ctrl+c++ to interrupt the log scrolling (it does not affect the tool installation)