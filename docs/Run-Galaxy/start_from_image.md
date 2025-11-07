
[comment]: <> (# You are fucked up, have a :beer: :stuck_out_tongue_winking_eye:)
[//]: <> (this is another hidden comment)
In Case of bad issue or impossibility of deploying your instance with the instructions
in this manual, we have provisioned an `Image` of a VM which is equivalent to the one you
should be able to build during the training.

At any time you can deployed this Galaxy server, and the instructions to do it are bellow.

However, to be able to access to this image, you must provide the Gmail address which is
necessarily associated to your project/billing account (even if your academic email address
was required to obtain you Google coupon).

Remember: to see the VM `Image` which will be used here to deploy a Galaxy server, ==send
me, or provide me by any mean (Slack, GitHub, etc.), your Gmail adress identifier==.

# Starting a new VM instance from the image `analyse-genomes-v6`

## 1. open the form for creating an instance
This will be mostly done as in [this section](GalaxyKickStart.md#deployment) which we
are reproducing bellow, except that for the `Disque de démarrage` (`Boot disk`), we are
not going to choose `IMAGES PUBLIQUES`/`PUBLIC IMAGES` :

!!! info "Google Instance from the **Image** `analyse-genomes-v6`"
    - Name: `ansible-galaxy`
    - Region `europe-west6 (Zurich)` (or any region available with you Google coupon)
    - Zone: `europe-west6-a` (or `-b` or `-c`)
    - **Configuration de la machine**
        - `OPTIMISEE POUR LE CALCUL` (or `USAGE GENERAL`)
        - Série: `C2`
        - Type de machine: `c2-standard-16 (16 processeurs virtuels, 64 Go de mémoire)`
    - **Disque de démarrage** ==Here, this is changing==
        - This time, select the tab `CUSTOM IMAGES`/`IMAGES PERSONNALISÉES`
        - Click the button `SELECT A PROJECT`:`SELECTIONNER UN PROJET`
        - Click on the tab `ALL`: A this stage _you should see the project of your instructor_
          `Analyse Genome Coupon 1`, if, and only if you provided him with your Gmail
          address identifier of your project, and if he did not forget to record it :innocent:
        - Select `Analyse Genome coupon 1`
        - Now if you drop-down the menu `Image *`, you should immediately be able to select
          `analyse-genome-v3 (Created on nov. 28, 2021, 12:08:27)` (or another version if I
          upgraded the image)
        - Automatically, the boot disk type should be `Balanced persistent disk` and the
          Size (GB) should be `200`. Leave it that way.
        - ==SELECTONNER/SELECT== (no need to modify the advanced configuration)
    - **Pare-feu** (Do not forget, back to the main instance panel)
        - Check `Autoriser le trafic HTTP`/`Allow HTTP traffic`
    - Click on the `CREATE`/`CREER` button !


## 2. Deploy !

Your VM should take just a little more time to deploy, otherwise it is like the start of
a regular VM instance.

Using this Image, you will benefit from already provisioned (and uncompressed) references,
small RNA datasets and RNA datasets, as well as already buit dataset collections !

# <center>:nail_care:</center>
:warning: Do not forget to attach your static IP address to this VM (as described
[here](Preparing_reference.md#a-reserve-a-static-ip-address)) if you wish to use it for
your analyses.


