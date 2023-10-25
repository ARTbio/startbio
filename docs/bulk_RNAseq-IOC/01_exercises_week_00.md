---
hide:
  - navigation
---

### Slack

During the incoming week, please communicate exclusively with us using the IOC Slack workspace.

We will guide you on the best practices to adopt with Slack, and show you
several tips to exchange information between all participants or between only two participants.

If you have any problem to perform the exercices below, there is only one address,
the [IOC Bulk RNAseq Slack workspace](https://iocbulkrnaseq.slack.com)

### Exercises with GitHub (web interface only)

- Go to https://github.com/ARTbio/ARTbio_064_IOC_Bulk-RNAseq
- Check that you are invited to contribute to the repo (write in it) through a Pull request
- Create a readme.md file in this path /participants/<yourname>/readme.md
  
  Note that the subdirectory <yourname> will be automatically created with you state the
  path for the new file you wish to edit
- Upload a galaxy workflow file (the one that you have generated during the first online
  meeting, with an extension .ga) in your newly created folder.


### Data upload in PSILO, then in Galaxy from Psilo

During this first week, the essential work is to gather your data a some places,
such as you will be confortable to manipulate and analyse them.

- Start by uploading your raw/starting data on you Psilo account
  Put everything you need in your PSILO account: fastq and fasta files, genome fasta
  reference(s), tables of metadata describing your data (a data without metadata is most
  often useless).
  
  Don't clutter your mind with unnecessary mental maps: if your data is all in one place,
  you will have no problem finding it later!
- Then, practice downloading your PSILO data into a Galaxy story that you will name
  “Input dataset”

### Pretreatments: data and metadata organisation in Galaxy

####First of all, follow the
[training to collection operations](https://training.galaxyproject.org/training-material/topics/galaxy-interface/tutorials/collections/tutorial.html)
prepared by the Galaxy Training Network.

**Do it in your Galaxy account, of course !**

In this training (plan at least one hour, this is a minimum), play a particular attention
to

- The section "**Relabel identifiers**", in the submenu `Tools that manipulate elements
within a collection`
- A tool that *is not mentioned in the training* but that we find very useful:
  **Extract element identifiers** (of a list collection). Test this tool, first on the
  list dataset 16: `Call variants on collection 11`, then on the paired-list dataset:9
  `M117-collection`.

#### To do
Post one of the two dataset/file generated above in you GitHub folder ! 

#### Last But Not least

1. Transfer all the files you will need for you analysis from PSILO, to a Galaxy history
which your are going to name `Input datasets`
- Group your datasets in collections as you have learned in the
  [Galaxy Training](https://training.galaxyproject.org/training-material/topics/galaxy-interface/tutorials/collections/tutorial.html)
  
  You may have to create either list of datasets or list of dataset pairs, depending on
  data structures.
  
### Good Work ! :construction_worker: :construction_worker: :construction_worker:
----