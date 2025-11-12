## Management of your Google Virtual Machine

By now, you have likely launched your Galaxy Virtual Machine at least once.
The next section will cover how to access the web interface of your Galaxy server.

Before we continue, here are a few important guidelines to follow for the rest of the training:

- [x] ==**Do not** stop your VM, instead **suspend** it==
  
  Stopping your VM is like stopping your PC or you laptop.
  
  You will stop everything and will have to literally reboot everything, including
  the Galaxy server. It is not that difficult actually, but it takes a bit more time.
  
  Instead, **Suspending** your VM is like putting your PC in sleeping mode, or closing
  the lid of your laptop.
  
  Thus, remember, at the end of the day or whenever you are not going to use your VM for a long
  time, use:
  
  ![suspend](images/suspend.png)
  
- [x] **==Protect your instance from unwanted destruction==**
    
    Accidents happen quickly...
  
    - Go to the Compute Engine management web page.
    - Click on the **name** of your VM.
    - Click on the top menu :pencil2:`Modifier`
    - Check the `Activer la protection contre la suppression` box (At the end of the section `Informations générales`).
        
    ![](images/self_destruction.png){width=600}
    
    - [x] and do not forget to save this new setting (button `Enregistrer` at the bottom of the page).
  
  From this point, you will need to uncheck the box to destroy the instance and your are
  protected against unwanted manifestations of bad karma. :imp:

---