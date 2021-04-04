## Stockage et transfert de fichiers volumineux


# Faire une introduction,  Parler de sécurité - de chriffrement des données...
Le chiffrement avec des outils comme ZED! 

### Les clouds
  - [France Grilles](http://www.france-grilles.fr/catalogue-de-services/fg-cloud/)
  - [Eudat](https://b2drop.eudat.eu) Cloud européen qui met à disposition 20 Go à tout le monde.
  - Drive (Google)
  - Dropbox
  
  


### Les transferts


  - [FileSender (Renater)](https://www.renater.fr/fr/filesender)
  - [MyCore (CNRS)](https://mycore.core-cloud.net/)
  - (Inserm)
  

  Liste d'outils gratuits : https://www.codeur.com/blog/envoi-gros-fichier/
  
  - [Smash](https://fromsmash.com/) Pas de limite en Go
  - [TransferNow](https://www.transfernow.net/) 5Go 
  - [WeTransfer](https://wetransfer.com/) 2Go
  - [pCloud](https://transfer.pcloud.com/)
  - [FileMail](https://fr.filemail.com/) 
  - [GrosFichiers](https://www.grosfichiers.com/) 4Go
  - [Sendbox](https://www.sendbox.fr/) 3Go
  - [Firefox Send](https://send.firefox.com/) 
  - [Free](http://dl.free.fr)
  
  
  
  
  Aussi creuser un peu la dedans :
 On
peu faire ça aussi en python, des serveurs web python existent
(https://bogotobogo.com/python/python_network_programming_server_client_file_transfer.php
,
https://www.ostechnix.com/how-to-setup-a-file-server-in-minutes-using-python/).
En bas niveau, sans serveur web, on peut aussi faire du transfert avec
netcat directement :
https://tutorials.technology/tutorials/How-to-transfer-files-over-the-network-using-Netcat.html
(faire des tar.gz/checksums avant). Sinon, je suis tombé là-dessus
récemment : https://github.com/warner/magic-wormhole , je ne sais pas ce
que ça vaut.
> >>
C'est un outil très efficace et extrêmement bien sécurisé. Parfait sous
Linux et MacOS, mais assez difficile à installer sous Windows (pour moi,
je ne suis pas un expert).

C'est pour un envoi quasi-synchrone : vous lancez le transfert, vous
appelez le destinataire pour lui donner le mot de passe, et il doit
lancer l'outil de son côté pour récupérer le fichier.

>>>>

>>>

dans la même veine que magic-wormhole, il y a croc:
- https://github.com/schollz/croc

(dont l'installation et/ou la compilation croisée est bien plus aisée que magic-wormhole)
>>>



Il existe d'autres solutions, mais il faudrait que les admins de labo
s'y collent, s'ils n'en ont pas déjà (serveurs lufi (voir une liste ici
https://pads.tedomum.net/p/T%C3%A9l%C3%A9travail_et_les_outils_libres),
serveurs (s)ftp(s), serveur http(s) de partage, etc...).
