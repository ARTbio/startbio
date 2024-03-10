## Galaxy Workflows

At this point, you should be more familiar with

- importing and manipulating datasets in Galaxy
- using tools in single, consecutive steps
- visualising the metadata associated to inputs, computational steps, and outputs.

You could arguably point out that all of these actions can be performed (maybe faster) either
in a Linux terminal (for Linux tools), the R environment (for R packages), or in a python
environment for python scripts.

It can be noted, however, that using several completely separate environments would make
the analysis difficult to understand, compared to reading an analysis in a single Galaxy
history.

Much worse, if you opt to use multiple environments with command lines you will
not maintain the links that connect the inputs, the computational tool and the outputs and
you will have to guess them based on their plausibility. On the contrary, in a Galaxy
hisstory, all these links are kept in a database (postgresql) and they can be retrieved
(even years later) by clicking on the galaxy datasets information icon.

Having said that, the accumulation of computational steps in histories is not the
culmination of an argument in favor of Galaxy.

You've likely noticed that analysis histories can become quite complex. Numerous
trial-and-error iterations and datasets accumulate, making it difficult to recall their
origins and purposes after a surprisingly short period.

Scripting these analyses into
computational workflows offers the most effective solution for preserving them.

**These workflows are the foundation of Galaxy, that streamlines their creation, execution,
and management**.

### Building and Sharing Analyses with Galaxy Workflows

Galaxy workflows offer a powerful solution for managing complex analyses. You can either:

- [x] Extract a workflow from an existing history:
    
    This captures the steps you've taken in your analysis, making it easy to replicate.
    
- [x] Build a workflow from scratch using the Galaxy workflow editor:
    
    This allows you to design custom workflows for specific analyses.

The beauty of workflows lies in their reusability. You can:

- [x] Replay a workflow at any time:
    
    Simply run the workflow again to regenerate your analysis, saving time and effort.
    
- [x] Export workflows as shareable .ga files:
    
    This allows you to export your workflows and import them into other Galaxy servers. As
    long as the new server has the required data and tools, the analysis will run identically.

### Key Takeaway
Advanced Galaxy users leverage workflows to capture their analyses, ensuring transparency,
reproducibility, and reusability of their computational protocols.

### Looking Ahead:
The next section will explore... (insert what the next section covers).

, you will test 2 workflows that are available in your
Galaxy server and recapitulate most of the analyses you have performed today.