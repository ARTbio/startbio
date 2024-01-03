![](images/galaxylogo.png)

# `DESeq2`

----
![](images/tool_small.png)

  1. Let's create a clean fresh history (`wheel` --> `Create New`) and name it DESeq2 ![](images/wheel.png)
  2. Copy the `.Counts`datasets from your `STAR`/ `HISAT2` history to this new history
  (`wheel` --> `Copy datasets`)
  3. Select the `DESeq2` tool with the following parameters:
      1. `how`: Select group tags corresponding to levels
      2. In `Factor`:
          1. In `1: Factor`
              - `Specify a factor name`: Treatment
              - In `Factor level`:
                  - In `1: Factor level`:
                      - `Specify a factor level`: treated
                      - `Counts file(s)`: the 3 gene count files with `treat` in their name
                  - In `2: Factor level`:
                      - `Specify a factor level`: untreated
                      - `Counts file(s)`: the 4 gene count files with `untreat` in their name
          2. Click on `Insert Factor` (not on `Insert Factor level`)
          3. In `2: Factor`
              - `Specify a factor name` to Sequencing
              - In `Factor level`:
                  - In `1: Factor level`:
                      - `Specify a factor level`: Paired
                      - `Counts file(s)`: the 4 gene count files with `paired` in their name
                  - In `2: Factor level`:
                      - `Specify a factor level`: Single
                      - `Counts file(s)`: the 3 gene count files with `single` in their name
      3. `Files have header?`: Yes
      4. `Output normalized counts table`: Yes
      5. `Execute`
  
