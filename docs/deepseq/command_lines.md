![](images/toolbox.jpg){width="150"}

## 1. ssh Connection to the Analysis Server

- Open/Start your ssh terminal (PuTTY, MacOS terminal, Chrome “Secure Shell” extension, Firefox “sshGate” extension)

- Type the command line in your terminal
```
ssh <prenom en minuscule>@35.205.151.231
```
and press the ++enter++ key to run the command.

For instance, Maëva is going to type: `ssh maëva@35.205.151.231` in her terminal, followed
by ++enter++.

:warning: Be careful, sometime your first name is written with accents and sometime not...
To be sure you should open this [shared Google doc](https://docs.google.com/spreadsheets/d/
1BYd5O-p94l6ToZPOwCAi7XcPtEdnq1FR/edit?usp=sharing&ouid=115530273158026068814&rtpof=true&sd=true)
where you are going to find ==your== personalised command for ssh connection to the server
with IP address `35.205.151.231`.

If your first command is successful, you should be prompted with the following message:
```
The authenticity of host '35.205.151.231 (35.205.151.231)' can't be established.
ECDSA key fingerprint is SHA256:T12UFAjo+DBgo5LT8CXoEMKI9MQCw57iPf+1SH+zPYI.
Are you sure you want to continue connecting (yes/no/[fingerprint])?
```
You can trust this server and answer `yes`

Then you will be asked for your password. Copy and paste this password from the Google sheets
document and press the ++return++ key.

:warning: your password typing is blind, this is normal.

From there, you are connected to you `home` directory, which is also symbolised by `~`.
Your `prompt` the text just before the command typing area, is constructed as follows:
```
carolina@instance-1:~$
```
where, for instance here, `carolina`is the login, `instance-1` is the name of the machine,
and the path is `~`, ie your home directory.

:warning: For the rest of the training, the commands will be always indicated in command
fields (the grey box with a copy/paste icon in the upper right corner), without the `prompt`.

## 2. Basic "navigation"

- Type (or copy and paste)
```
pwd
```
in your terminal and press ++enter++.

??? question "What do you read ?"
    you see the full path of your current position in your file system. `pwd` stands for
    **print working directory**
    
    You should be in `/home/<user name>` !

- Type (or copy and paste)
```
ls -la
```
++enter++
??? question "What do you read ?"
    the `ls` command print the elements in your current (working) directory.
    
    you expect something very close to:
    ```
    total 24
    drwxr-xr-x  3 carolina students 4096 Dec  7 23:58 .
    drwxr-xr-x 45 root     root     4096 Dec  8 00:38 ..
    -rw-r--r--  1 carolina students  220 Feb 25  2020 .bash_logout
    -rw-r--r--  1 carolina students 3771 Feb 25  2020 .bashrc
    drwx------  2 carolina students 4096 Dec  7 23:58 .cache
    -rw-r--r--  1 carolina students  807 Feb 25  2020 .profile
    ```
    Here, the options -la are formatting the output of `ls` in a certain way.
    
    If you want to know the options available for a command, type:
    ```
    man ls
    ```
    You can also try:
    ```
    ls
    ```
    alone and see how it affects the output of the command.

- Type (or Copy/Paste):
```
cp ../GKG-13.fastq.gz .
```
++enter++
??? question "What have you done with this command ?"
    You have copied the file `GKG-13.fastq.gz` which is located in the parent directory of
    your current directory (ie `/home/`) to your current directory.
    
    Note that in POSIX commands, `..` stand for the directory "above" your current directory,
    whereas `.` (a single dot) stands for the current directory.

- Type (or Copy/Paste):
```
ls -la
```
??? question "What do you read and what can you deduce from it?"
    You now see a _copy_ of the GKG-13.fastq.gz file in your home directory !

- Type (or Copy/Paste):
```
gunzip GKG-13.fastq.gz
```
??? question "What is it doing ?"
    the GKG-13.fastq.gz file is a compressed file (format gzip), as indicated by the `.gz`
    extension.
    The gunzip command has uncompressed the file to `GKG-13.fastq`
    
    You can verify it by typing
    ```
    ll
    ```
    which should show:
    ```
    drwxr-xr-x  3 carolina students       4096 Dec  8 01:26 ./
    drwxr-xr-x 45 root     root           4096 Dec  8 00:38 ../
    -rw-r--r--  1 carolina students        220 Feb 25  2020 .bash_logout
    -rw-r--r--  1 carolina students       3771 Feb 25  2020 .bashrc
    drwx------  2 carolina students       4096 Dec  7 23:58 .cache/
    -rw-r--r--  1 carolina students        807 Feb 25  2020 .profile
    -rwxr-xr-x  1 carolina students 1043785660 Dec  8 01:25 GKG-13.fastq*
    ```
    :warning: as you can deduce, the `ll` command is an _alias_ to `ls -laF`

## 3. What is this fastq file containing ?

Type (or Copy/Paste):
```
more GKG-13.fastq
```
??? question "What is doing the `more` command ?"
    It reads the file by chunks of your screen size. Each time you type the ++space++ bar,
    you print the next chunk.
    To exit from this read mode, just press the ++q++ key.

??? question "How is made a fastq file ?"
    ```
    @HWIEAS210R_0028:2:1:3019:1114#AGAAGA/1 Header
    TNGGAACTTCATACCGTGCTCTCTGTAGGCACCATCAA  Sequence
    +HWIEAS210R_0028:2:1:3019:1114#AGAAGA/1 Header
    bBb`bfffffhhhhhhhhhhhhhhhhhhhfhhhhhhgh  Sequence Quality (ASCII encoded)
    @HWIEAS210R_0028:2:1:3925:1114#AGAAGA/1
    TNCTTGGACTACATATGGTTGAGGGTTGTACTGTAGGC
    +HWIEAS210R_0028:2:1:3925:1114#AGAAGA/1
    ]B]VWaaaaaagggfggggggcggggegdgfgeggbab
    @HWIEAS210R_0028:2:1:6220:1114#AGAAGA/1
    TNGGAACTTCATACCGTGCTCTCTGTAGGCACCATCAA
    +HWIEAS210R_0028:2:1:6220:1114#AGAAGA/1
    aB^^afffffhhhhhhhhhhhhhhhhhhhhhhhchhhh
    @HWIEAS210R_0028:2:1:6252:1115#AGAAGA/1
    TNCTTGGACTACATATGGTTGAGGGTTGTACTGTAGGC
    +HWIEAS210R_0028:2:1:6252:1115#AGAAGA/1
    aBa^\ddeeehhhhhhhhhhhhhhhhghhhhhhhefff
    @HWIEAS210R_0028:2:1:6534:1114#AGAAGA/1
    TNAATGCACTATCTGGTACGACTGTAGGCACCATCAAT
    +HWIEAS210R_0028:2:1:6534:1114#AGAAGA/1
    aB\^^eeeeegcggfffffffcfffgcgcfffffR^^]
    @HWIEAS210R_0028:2:1:8869:1114#AGAAGA/1
    GNGGACTGAAGTGGAGCTGTAGGCACCATCAATAGATC
    +HWIEAS210R_0028:2:1:8869:1114#AGAAGA/1
    aBaaaeeeeehhhhhhhhhhhhfgfhhgfhhhhgga^^
    ```
    One sequence read is encoded by a block of ==**4**== lines, 1 for the header (the name
    of the read), 1 for the nucleotide sequence, 1 starting with a `+` which is a copy of
    the header, and the last line for the quality of the base calling at each position,
    encoded by an ASCII character.

## 4. How many sequence reads in my file ?
Type (or Copy/Paste):
```
wc -l GKG-13.fastq
```
!!! info "the `wc` command"
    prints the number of newlines, words, and bytes for each file in argument (here,
    GKG-13.fastq)
    
    :warning: As you used the `-l` option, you only print the number of newlines in the
    file.

??? question "And then... How many _**sequences**_ in GKG-13.fastq ?"
    
     **6425957**, NOT 25703828, because each sequence read is encoded by 4 lines !

## 5. Are my sequence reads containing the adapter ?

This fastq file corresponds to the sequencing of a small RNA library, whose 3' adapter
contains the sequence 5’-**CTGTAGG**CACCATCAAT-3’
Type (or Copy/Paste):
```
cat GKG-13.fastq | grep CTGTAGG | wc -l
```
This should return:
```
6355061
carolina@instance-1:~$
```

??? Question "A lot of things to comment in the previous command !"
    - `cat` print the _total_ content of the file in argument
    - the sign ++pipe++ is important, we call it the `pipe`. the `|` takes the ==output== of
    the _upstream_ command (here `cat GKG-13.fastq`) and gives it as ==input== to the
    downstream command (here `grep CTGTAGG`).
    - `grep` prints the lines of the input (or the argument if used without pipe) _only_ if
    these lines contains the string `CTGTAGG`.
    - the second pipe `|` sign sends the output of the grep command and counts the number
    of lines in _this_ output.
    
    Brillant isn't it ?

### Why doing it simple when you can do it complicated ?

Yeah... As a matter of fact, you can obtain exactly the same information from `GKG-13.fastq`
by typing (or copying and pasting) the command:
```
grep -c "CTGTAGG" GKG-13.fastq
```

:warning: Here, the option -c is passed to `grep` to ask for only counting and not _printing_
the lines that contain the string pattern `CTGTAGG`

Check it out !


