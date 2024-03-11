# A workflow of your use-case

The exercise of this week is difficult:

You are going to prepare a complete workflow of your analysis.

Depending on your model organisms, you may not have been able to perform all of the
analyses covered in this training. This is not a problem: you are expected to create a
workflow from what you have actually been able to do.

In order to make a sustainable, reproducible and transparent workflow, you should meet the
following requirements:

## Workflow inputs

Best inputs are

- [x] Completely unprocessed data (i.e. fastq files)
- [x] Preferably accessible through a sustainable URL. If it is not possible, they should
  be at least easily accessible (i.e. gathered in a single folder, whose location is
  precisely described)
- [x] reference data (GTF, bed, etc...) should be precisely annotated, date, organisation,
  version, etc... Importantly, a **direct** URL to the original reference should be included
- [x] :warning: Unless impossible to do, do not use processed data as inputs of your
  workflow. If you think this is impossible to do, **let's discuss it** !
- A lot of good workflows stand on a metadata table, which describes input data, their
  names, labels if required, replicate status, etc. This metadata table may be considered
  as a genuine dataset which can be used by the workflow to perform some operations.
  
## Computational steps

- [x] Whenever a computational step applies to multiple sample, think "**Collections**"
- [x] A good clue that you should switch to collections is when your workflow contains
  twice or more the same step with the same parameters (or almost the same)
- [x] Take the time, for each step, to carefully fill the tool form at the right hand-side
  of the workflow editor.
- [x] There are several fields in this tool form that *must* be used to clarify the step:
  The `Label` field at the top of the tool form, the `Step Annotation` field, and the
  `Configure Output: xxx` fields and their sub-fields `Label`, `Rename dataset` and `Change
  datatype`
  
  Experiment theses fields with your workflow !
  
- [x] Workflow **can use parameters** at their runtime. If you are interested by this functionality,
  let's discuss it !
  
## Workflow outputs

- [x] You can hide some output datasets for better readability of the workflow by
  unchecking this outputs in the tool items of the workflow.
      
      :warning: By default all outputs are visible although unchecked. This is only when you
      check a first output that unchecked outputs become hidden.
      
      :warning: Hidden does not mean deleted: all workflow outputs are still there and you can
      reveal them in the Galaxy history.
  
- [x] Whenever possible, rename your datasets in the workflow using the `Configure Output: xxx`
  fields in the tool forms
  
## Your objective:

Is that you generate the complete analysis in a **single** workflow run, with the minimal
number of inputs.

This way, you can even loose/trash your Galaxy history :
Just having the inputs plus the workflow should be enough to regenerate the analysis.

Consider that it is also a **huge** gain in term of data storage.
