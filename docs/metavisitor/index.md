![metavisitor_logo](images/metavisitor_logo.png)


[Metavisitor](http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0168397) is
a user-friendly and adaptable software to provide biologists, clinical researchers and
possibly diagnostic clinicians with the ability to robustly detect and reconstruct viral
genomes from complex deep sequence datasets. A set of modular bioinformatic tools and
workflows was implemented as the Metavisitor package in the Galaxy framework. Using the
graphical Galaxy workflow editor, users with minimal computational skills can use existing
Metavisitor workflows or adapt them to suit specific needs by adding or modifying analysis
modules.

# Reference viral database

Metavisitor's workflows use a home-made reference viral database **vir2**. This database was
made using *Galaxy-Workflow-Metavisitor__Workflow_for_nucleic_vir2_generation* and
*Galaxy-Workflow-Metavisitor__Workflow_for_proteic_vir2_generation*, that can both be found
in [Metavisitor's Github](https://github.com/ARTbio/metavisitor/tree/master/extra-files/metavisitor).

#### How was *nucleic vir2* generated?

- Downloading every viral sequence from NCBI's *nuccore* database with the following
  queries (2018/03/21):

  >*txid10239[Organism] NOT txid131567[Organism] NOT phage[All Fields] NOT patent[All Fields] NOT chimeric[Title] NOT vector[Title] NOT method[Title] NOT X174[All Fields] AND 301:10000[Sequence length]*

  >*txid10239[Organism] NOT txid131567[Organism] NOT phage[All Fields] NOT patent[All Fields] NOT chimeric[Title] NOT vector[Title] NOT method[Title] NOT X174[All Fields] AND 10001:1300000[Sequence length]*

- Clustering sequences with 95% identity and shorter than 10 001 bp using vclust.

**vir2** is  available for download in [Figshare](https://figshare.com/articles/vir2_NCBI_21-03-2018/6106892)

# Quick Start

Users who want to use Metavisitor on the [Galaxy Mississippi Server](https://mississippi.snv.jussieu.fr),
or got already the Metavisitor suite of tools installed in their own Galaxy server, can
jump to the next chapter [Prepare input data histories](use_cases_input_data).


## Availability of Metavisitor tools and workflows

Metavisitor has been developed at the [ARTbio platform](http://artbio.fr). Its tools are
primarily available in [GitHub](https://github.com/ARTbio/tools-artbio). Its workflows are
primarily available in the
[metavisitor repository](https://github.com/ARTbio/metavisitor/tree/master/extra-files/metavisitor)

Metavisitor tools and workflows are also available in the [toolshed](https://toolshed.g2.bx.psu.edu/)

#### Metavisitor tools developed by ARTbio in the [ARTbio GitHub](https://github.com/ARTbio/tools-artbio)

- [`blast_to_scaffold`](https://github.com/ARTbio/tools-artbio/tree/master/tools/blast_to_scaffold)
- [`blastx_to_scaffold`](https://github.com/ARTbio/tools-artbio/tree/master/tools/blastx_to_scaffold)
- [`blastparser_and_hits`](https://github.com/ARTbio/tools-artbio/tree/master/tools/blastparser_and_hits)
- [`blast_unmatched`](https://github.com/ARTbio/tools-artbio/tree/master/tools/blast_unmatched)
- [`cap3`](https://github.com/ARTbio/tools-artbio/tree/master/tools/cap3)
- [`cherry_pick_fasta`](https://github.com/ARTbio/tools-artbio/tree/master/tools/cherry_pick_fasta)
- [`cat_multi_datasets`](https://github.com/ARTbio/tools-artbio/tree/master/tools/concat_multi_datasets)
- [`fetch_fasta_from_ncbi`](https://github.com/ARTbio/tools-artbio/tree/master/tools/fetch_fasta_from_ncbi)
- [`oases`](https://github.com/ARTbio/tools-artbio/tree/master/tools/oases)
- [`sequence_format_converter`](https://github.com/ARTbio/tools-artbio/tree/master/tools/sequence_format_converter)
- [`small_rna_maps`](https://github.com/ARTbio/tools-artbio/tree/master/tools/small_rna_maps)
- [`sr_bowtie`](https://github.com/ARTbio/tools-artbio/tree/master/tools/sr_bowtie)
- [`yac_clipper`](https://github.com/ARTbio/tools-artbio/tree/master/tools/yac_clipper)

Tools from other developers are used in the suite metavisitor-2. These tools are
available from the [main Galaxy toolshed](https://toolshed.g2.bx.psu.edu/):

     name="bowtie2" owner="devteam"
     name="data_manager_bowtie2_index_builder" owner="devteam"
     name="data_manager_fetch_genome_dbkeys_all_fasta" owner="devteam"
     name="fasta_compute_length" owner="devteam"
     name="fasta_filter_by_length" owner="devteam"
     name="fastx_trimmer" owner="devteam"
     name="ncbi_blast_plus" owner="devteam"
     name="data_manager_bowtie_index_builder" owner="iuc"
     name="khmer_normalize_by_median" owner="iuc"
     name="sra_tools" owner="iuc"
     name="trinity" owner="iuc"
     name="vsearch" owner="iuc"
     name="regex_find_replace" owner="galaxyp"
     name="spades" owner="nml"

#### Availability of Metavisitor tools and workflows for **Galaxy instance administrators**

All metavisitor tools are available from the
[suite_metavisitor_2](https://toolshed.g2.bx.psu.edu/view/artbio/suite_metavisitor_2/1570b18266be)
Galaxy Admin can just install this suite of tools by using the `Install new tools` menu in
their Admin panel, searching for "metavisitor", and installing the `suite_metavisitor_2`
tool suite.

Galaxy Admins can install the workflows from the `metavisitor_workflows` repository in the
[main Galaxy toolshed](https://toolshed.g2.bx.psu.edu/), which will
install in addition all tools needed for Metavisitor.

#### Availability of Metavisitors workflows for any Galaxy instance user.
We have deposited the Metavisitors workflows in the
[myexperiment server](http://www.myexperiment.org/workflows #TODO: this will need to be updated too),
where they are searchable with "metavisitor" keyword and can be downloaded and reuploaded
to the Galaxy instance.

## Starting a Metavisitor Galaxy server from scratch

In the [last section](https://artbio.github.io/metavisitor/install_metavisitor/) of this
documentation, we provide instructions to set up Galaxy server instances *from scratch*
with *pre-installed* Metavisitor tools and workflows:

- Based on Ansible: see [Metavisitor with GalaxyKickstart (Ansible)](metavisitor_ansible.md)
- Based on Docker: see [Metavitor with Docker](metavisitor_docker.md)
