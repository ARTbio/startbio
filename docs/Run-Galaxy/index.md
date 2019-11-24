This training documentation is coded in this GitHub
[repository](https://github.com/ARTbio/Run-Galaxy)

# Why Running Galaxy as an administrator ?

You may be wondering: "Why doing all this geeky IT stuff when I have access to Galaxy
servers administrated by professional ?"

It is true that there is a lot of powerful Galaxy instances, and at first, 
the [main Galaxy instance](https://usegalaxy.org/). The expanding list of
public galaxy servers is available [here](https://galaxyproject.org/public-galaxy-servers/).

However, a number of issues can be successfully addressed
if you are able to administrate your own Galaxy server, including:

1. **Storage/Disk Space.**
    
    Most of Public Galaxy Servers provide their users with a quota that rarely exceed 200-300
    Giga-bytes. Although this may seem a lot, it is not unfrequent that analyses that deal with
    numerous samples require 1 Tera-bytes or more.
    
    When you administrate your Galaxy server, you control your storage space. Of course,
    since nothing in free in this world, keep in mind that you will have to support for the
    cost of this storage.
    
2. **Isolation.**
    
    If you control a Galaxy server for a given analysis project and only for this project,
    you can argue that you benefit from an analysis environment that is isolated.
    
3. **Accessibility and Reproducibility**
    
    Whenever you need to give access to collaborators or reviewers to your work, giving access
    to your Galaxy server is enough to provide high-quality transparency and reproducibility.
    This is far better than just sharing public histories, since when you are not administrator,
    you do not have access to all computational details that are logged for Galaxy admins.
    Moreover, if you deploy you Galaxy server in a virtual environment (VM or docker containers)
    you can preserve the whole environment in an archive and redeploy this environment latter
    on in another infrastructure.
    
4. **Computational Resources.**
    
    Galaxy public servers are generally hosted in high performance computing
    infrastructures whose resources are shared between users.
    
    For instance, the main Galaxy server is hosted by
    [a network of US supercomputers](https://galaxyproject.org/main/). Nevertheless, the
    computational walltime for a user to execute standards analyses (BWA, bowtie, Tophat,
    Trinity, etc.) may exceed 5 or 6 hours.
    
    Likewise, some metagenomic or *de novo* assembly approaches may require a substantial
    amount of memory that is not necessarily provided by public Galaxy server.
    
- **Full control on installed tools**
    
    You may need a particular combination of tools for your analysis, and this combination
    may not be available in any public server. Although Galaxy admin are generally happy to
    install new tools for their users, other considerations that have to be taken into account
    in a public resources may limit installation of new tools: "not considered as harmless for
    the server", "to much resource-demanding for the infrastructure", "unable to provide support
    to the users of this tool", "not in the policy of the thematic Galaxy server", etc.
    
    When you administrate your Galaxy server, you can install any tool you need.
    You can even modify tools, or code your own tools and test these tools in live in your
    Galaxy instance.
    
    Last, but not least, when you are administrator, you have access to information on tool
    & workflow runs you cannot access to when you are regular users (some metadata, including
    running times, command lines, etc.)
    
- **Full Control on computational workflows.**
    
    Galaxy workflows can be exchanged between researchers and between Galaxy instances.
    However, to be effective, this interoperability requires that the tools called by an
    imported workflow are installed in the new Galaxy instance.
    You can only do that if your are administrator of this Galaxy instance.

    
- **Help your community.**
    
    Galaxy server administration is a very useful expertise: you can greatly
    help your colleagues if you are able to run a Galaxy server for them !
