## Management of your Google Virtual Machine

By now, you have likely launched and suspended (or stopped) your Galaxy Virtual Machine at least once.

Here are a few guidelines to follow for the rest of the training:

- [x] ==**Stopping** your VM versus ==**Suspending**== it.
  
  Stopping your VM is like stopping your PC or you laptop.
  
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


- [x] **==Protect your instance from unwanted destruction==**
    
    :warning: :warning: :warning: :warning: :warning: :warning: :warning: :warning: :warning: :warning: :warning: :warning: Accidents happen quickly...
  
    - Go to the Compute Engine management web page.
    - Click on the **name** of your VM.
    - Click on the top menu :pencil2:`Modifier`
    - Check the `Activer la protection contre la suppression` box (At the end of the section `Informations générales`).
        
    ![](images/self_destruction.png){width=600}
    
    - [x] and do not forget to save this new setting (button `Enregistrer` at the bottom of the page).
  
  From this point, you will need to uncheck the box to destroy the instance and your are
  protected against unwanted manifestations of bad karma. :imp:

---