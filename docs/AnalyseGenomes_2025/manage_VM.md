## Management of your Google Virtual Machine

By now, you have likely launched and suspended (or stopped) your Galaxy Virtual Machine at least once.

Here are a few guidelines to follow for the rest of the training.

### **1.==Stopping== your VM versus ==Suspending== it.**
  
  **Stopping** your VM is like stopping your PC or you laptop.
  
  You will stop everything and will have to literally reboot everything, including
  the Galaxy server. It is not that difficult actually, but it takes more time.
  
  Instead, **Suspending** your VM is like putting your PC in sleeping mode, or closing
  the lid of your laptop.
  
  In both Suspend and Stop cases, you stop paying for Google computing resources
  (this payment is occuring through your 50$ coupon)

  In contrast, you will continue to pay for the storage resources as long as you do not
  *destroy* your VM. ==But, please, don't do that before the end of the course !==
  
  Thus, remember, at the end of the day or whenever you are not going to use your VM for a long
  time, use:
  
  ![suspend](images/suspend.png)


### **2. ==Protect your instance from unwanted destruction==**
    
**:warning: :warning: :warning: Accidents happen quickly :warning: :warning: :warning:**
  
- Go to the Compute Engine management web page.
- Click on the **name** of your VM.
- Click on the top menu :pencil2:`Modifier`
- Check the `Activer la protection contre la suppression` box (At the end of the section `Informations générales`).
    
![](images/self_destruction.png){width=600}

- [x] and do not forget to save this new setting (button `Enregistrer` at the bottom of the page).
  
  From this point, you will need to uncheck the box to destroy the instance and your are
  protected against unwanted manifestations of bad karma. :imp:

### **3. Stalled jobs in your Galaxy histories**

When you trigger the `Run` button of a Galaxy tool, this generate one or several new datasets
in your current history, which will pass through different states that can be recognized from
their color.

**The normal life cycle of a dataset.**

|   Color                      |  State                 |    Events                                                                                        |
|------------------------------|------------------------|--------------------------------------------------------------------------------------------------|
|![](images/Slurm_drained.png) |Job waiting for starting|Decompression in the backgroud (if required). Attribution of the job to the slurm manager         |
|![](images/running_job.png)   |Job is running          |Slurm attributed the job to the compute node. Job currently executed                              |
|![](images/successful_job.png)|Successful Job          |Slurm sent back the job to Galaxy. Job is successful                                              |
|![](images/failed_job.png)    |Failed Job              |Tool is broken, Parameters were wrong (80% probability), memory overflow, etc.|                               

- [x] Failed jobs may seem freaky at first glance. However, do not panic.
      - Broken tools are very unlikely with your Galaxy server;
      - if the case happens, all students should experience the same problem. Same reasoning
      applies to the "memory overflow" hypothesis.
      - From my long experience in Galaxy server management, I can say that in most cases the problem arose
      from the use of wrong parameters in the tool form. Stay cool, review this parameters carefully (you can compare
      them with those of other students), change something in the form and run again.
- [x] The most tricky case is when your dataset seems stalled in "grey" state forever.
      - :warning: "long" is not "forever". Indeed, with some tools, a compressed (gz or zip) dataset must be
        uncompressed before to be sent to the slurm job manager. With large datasets (eg 1 GB), this may take 30" or 1 min.
      - Beyond 1 min of "grey" state, it is likely that something is going wrong with the communication between
        Galaxy and Slurm.
        Clearly, this situation may happen in the course of this Galaxy training, in particular if you just restart your
        Galaxy server.
- [x] **:warning: What you should do in this latest case:**
      - Let your "grey" dataset as it is. Do not delete it.
      - Open an ssh session in your VM, using the menu `ssh` --> `Ouvrir dans une fenêtre du navigateur`
      - Type `sudo -i` to have the root admin rights
      - Type `sinfo`. If the `STATE` is NOT `idle` or `mix` (eg `down*` or `Drain`), you found the problem:
        le slum node is indeed disconnected from its slurmctl controller and cannot
        communicate with Galaxy anymore.
      - Note or copy the name in the `NODELIST` (should be something like `batman-galaxy`)
      - Then type the following command, adapting the `<your-instance-name>` to your real case:
        ```
        scontrol update NodeName=<your-instance-name> State=Resume
        ```
      - type `sinfo` again. If the previous command solved the problem, the `STATE` should be now `Idle` or `mix`
      - ... and your "grey" dataset should soon turn to orange then green.


---