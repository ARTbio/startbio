## Management of your Google Virtual Machine

If you read this, you have probably launched at least one time a Google Virtual Machine.

A few rules to get your life easier during the rest of this training:

- ==**Avoid** stopping your VM, instead **suspend** it==
  
  Stopping your VM is like stopping your PC or you laptop.
  
  You will stop everything and will have to literally reboot everything, including
  the Galaxy server. It is not that difficult actually, but it takes a bit more time.
  
  Instead, **Suspending** your VM is like putting your PC in sleeping mode, or closing
  the lid of your laptop.
  
  Thus, remember, at the end of the day or whenever you are not going to use your VM for a long
  time, use:
  
  ![suspend](images/suspend.png)
  
- **==Protect your instance from unwanted destruction==**
    
    An accident happens so quickly...
  
    - Go to the Google Cloud Platform management web page.
    - Click on the **name** of your VM.
    - Click on the top menu :pencil2:`Modifier`
    - Edit the `Protection contre la suppression` option as follows:
    
    ![](images/self_destruction.png){width=600}
    
    (just at the end of the section **Informations générales**) and do not forget to save
    this new setting.
  
  From this point, you will need to uncheck the box to destroy the instance and your are
  protected against unwanted manifestations of bad karma :imp:!

---