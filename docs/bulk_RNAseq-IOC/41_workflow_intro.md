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
    
- [x] Use a combination of both approaches !
    
    Beginners tend to start with the first approach since it allows to automatically build
    a workflow without interacting too much with the workflow `editor`. However, in use this
    proves difficult, because the stories are often cluttered with several trials and
    errors or trials and successes, with different parameter values for the same tool.
    
    Thus, a workflow built from a story can be difficult to untangle.
    
    On the other hand, experts in using the workflow editor favor creating workflows from
    scratch. This mode requires you to have an analysis plan in mind, whereby workflow
    editing is literally akin to the graphic writing of a computer script. Testing this
    workflow can be done as it is written, by running it in Galaxy and verifying that the
    outputs are valid and conform to what is expected.
    
    In real life, it is often a combination of the two approaches that is implemented: you
    can start a workflow from a not too complicated story and correct / develop it later
    by first using the editor before testing it
    
    Along the same lines, Galaxy masters will also rely on already existing workflows to
    avoid reinventing what has already been done and save time. It is also possible to use 
    a workflow as a tool in another workflow, and thus to build very complex and elaborate
    workflows by structuring them as `workflows of workflows`.

The beauty of workflows lies in their reusability. You can:

- [x] Replay a workflow at any time:
    
    Simply run the workflow again to regenerate your analysis, saving time and effort.
    
- [x] Export workflows as shareable .ga files:
    
    This allows you to export your workflows and import them into other Galaxy servers. As
    long as the new server has the required data and tools, the analysis will run identically.

### Workflow reports
Another essential aspect of Galaxy workflows is that their invocations are logged and
accessible in the menu `User` --> `Workflow invocations`

In addition, a report is automatically generated for each workflow invocation. A minimal
default report is generated for each workflow invocation and give access to inputs, outputs
and the workflow ==in its runtime version==. You can customize and enrich this automated
report using the Galaxy workflow editor.

:warning: Reports cannot still be considered as a Material and Methods section for your
scientific manuscripts with computational analyses but they clearly make this section more
accurate and easier to write ! Moreover, the goal of reports is clearly to generate this
section in a fully automated manner, and Galaxy development is happening at a rapid pace !

### Key Takeaway
Advanced Galaxy users leverage workflows to capture their analyses, ensuring transparency,
reproducibility, and reusability of their computational protocols.
