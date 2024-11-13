!!! question "My question to [Bard](https://bard.google.com/chat)"
    HI,
    
    I am teaching a bioinformatics class using Google Cloud Education Credits for Students.
    Today a student trash its VM by error and I wanted to share with him a VM image prepared for this king of emergency situation..
    
    However, the student benefits from google cloud engine through its academic email addresse. I was not able to designate her as a recipient of image sharing rights since in the IAM board it is stated that a gmail address must be associated to the google cloud recipient.
    
    Do you know a turn-around for this issue ?

!!! info "The Bard response"
    Yes, there are a couple of workarounds to share the VM image with the student using
    their academic email address:
    
    1. **Create a shared project**: You can create a Google Cloud project specifically for this
       class and add both yourself and the student as owners. This will allow the student
       to access the VM image in the shared project. To do this, follow these steps:
       
         1. Go to the Google Cloud Console and create a new project.
         - Add yourself as an owner of the project.
         - Add the student's academic email address as an owner of the project.
         - Share the VM image with the shared project.
    2. **Use a service account**: You can create a service account specifically for this class
       and grant it access to the VM image. Then, you can share the service account's
       credentials with the student. To do this, follow these steps:
       
         1. Go to the Google Cloud Console and create a new service account.
         - Grant the service account access to the VM image.
         - Create a service account key for the service account.
         - Share the service account key's JSON file with the student.
