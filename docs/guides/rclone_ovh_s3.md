# Transfer data to an OVH S3 bucket with rclone

This guide shows how to copy data from a Linux server, a Mac or a Windows PC to an
OVHcloud Object Storage bucket (S3) with [rclone](https://rclone.org/), and how to check
that everything arrived intact. You always start with an exercise: you transfer dummy
data, generated for the occasion, with exactly the same commands as for your real data.
On a Mac or a PC, you can also try a graphical client, Cyberduck, on the same data.

Moving the data to cold archive afterwards is optional and covered separately at
the end of the page ([Part 2](#part-2-cold-archive-optional)).

## Your values

<div id="rclone-form" class="rclone-form" markdown>

First choose where you will type the commands, then fill in the fields: every command
on this page updates as you type, and you only have to copy it. Fields highlighted in
yellow still hold an example value. Your values are kept in this browser for your next
visit; nothing is sent anywhere.

The exercise data and the logs always go to your home folder, in the folders
`rclone-exercise` and `rclone-logs`: they need no field.

<div class="rclone-os" role="group" aria-label="Where you run the commands">
<button type="button" data-os-pick="linux">Linux server<small>in a terminal, often through SSH</small></button>
<button type="button" data-os-pick="mac">Mac<small>Terminal app</small></button>
<button type="button" data-os-pick="win">Windows<small>PowerShell, Windows 10 or 11</small></button>
</div>

<div class="rclone-fields"></div>

<button type="button" class="md-button rclone-reset">Reset this system's fields to example values</button>

!!! warning "Never type your access keys here"
    Keys are typed only in the terminal, at step 2 (and in Cyberduck's own login window
    if you try it).

Example: your values in the command that lists the content of the bucket. It changes
nothing, and works once steps 1 and 2 are done.

``` { .text .rclone-tpl }
rclone lsf {{remote}}:{{bucket}}
```

</div>

<p class="rclone-choose">Choose Linux server, Mac or Windows above: the steps below then show
the commands for your system only.</p>

## Before you start

- Three pieces of information from whoever manages your OVHcloud Public Cloud project:
  the bucket name (<span data-rclone="bucket">my-bucket</span>), an access key and a
  secret key.
- <span data-only="linux">An account on the Linux server that holds the data, with
  Internet access, and at least 2 GB of free space in your home directory for the
  exercise.</span><span data-only="mac">The Terminal app (Applications, then Utilities)
  and at least 2 GB of free space in your home folder for the exercise.</span><span
  data-only="win">Windows PowerShell (Start menu, type PowerShell) and at least 2 GB of
  free space in your home folder for the exercise.</span><span data-only="mac win">
  Keep the computer plugged in and awake during transfers.</span>
- Administrator rights on your machine to install rclone (`sudo` on Linux and Mac).
  <span data-only="linux">On a shared server where you do not have them, use the
  variant without administrator rights given at step 1.</span>

## Part 1 — Transfer the data

### Step 1 — Install rclone

<div data-only="linux mac" markdown>

With rclone's official installation script, which asks for your password (`sudo`) and
installs the latest version:

``` { .text .rclone-tpl }
sudo -v ; curl https://rclone.org/install.sh | sudo bash
rclone version
```

The last command must print a version number.

</div>

<div data-only="linux" markdown>

??? note "Without administrator rights (shared server)"

    This installs rclone in your own `~/bin` folder; it works on a standard x86-64
    server (`uname -m` prints `x86_64`).

    ``` { .text .rclone-tpl }
    mkdir -p ~/bin && cd /tmp
    curl -O https://downloads.rclone.org/rclone-current-linux-amd64.zip
    unzip -o rclone-current-linux-amd64.zip
    cp rclone-*-linux-amd64/rclone ~/bin/
    cd -
    echo 'export PATH="$HOME/bin:$PATH"' >> ~/.bashrc && source ~/.bashrc
    hash -r && which rclone && rclone version
    ```

    The output must show a path in your home directory ending in `/bin/rclone`, then a
    version number.

</div>

<div data-only="win" markdown>

With winget, which is part of recent Windows versions:

``` { .text .rclone-tpl }
winget install Rclone.Rclone
```

Then close PowerShell, open a new PowerShell window, and check:

``` { .text .rclone-tpl }
rclone version
```

</div>

### Step 2 — Configure access to the bucket

Paste this block. It asks for the access key, then for the secret key, which is not
displayed while you type. If a configuration with the same name already exists, the
block stops without changing anything: choose another name in "Your values" and copy
the block again. Other rclone configurations are left untouched.

<div data-only="linux mac" markdown>

``` { .text .rclone-tpl }
if rclone listremotes | grep -qx "{{remote}}:"; then
  echo "A configuration named {{remote}} already exists: choose another name in Your values."
  rclone listremotes
else
  printf "Paste the access key, then Enter: "; read -r AK
  printf "Paste the secret key (nothing is displayed, this is normal), then Enter: "
  stty -echo; read -r SK; stty echo; echo
  rclone config create {{remote}} s3 provider=Other env_auth=false \
    access_key_id="$AK" secret_access_key="$SK" \
    region={{region}} endpoint=https://{{endpointhost}} \
    acl=private no_check_bucket=true > /dev/null
  unset AK SK
  chmod 600 "$(rclone config file | tail -1)"
  rclone listremotes
fi
```

</div>

<div data-only="win" markdown>

``` { .text .rclone-tpl }
& {
  if ((rclone listremotes) -contains "{{remote}}:") {
    Write-Host "A configuration named {{remote}} already exists: choose another name in Your values."
  } else {
    $AK = Read-Host "Paste the access key, then Enter"
    $SS = Read-Host "Paste the secret key (nothing is displayed, this is normal), then Enter" -AsSecureString
    $SK = [Runtime.InteropServices.Marshal]::PtrToStringBSTR([Runtime.InteropServices.Marshal]::SecureStringToBSTR($SS))
    rclone config create {{remote}} s3 provider=Other env_auth=false access_key_id=$AK secret_access_key=$SK region={{region}} endpoint=https://{{endpointhost}} acl=private no_check_bucket=true | Out-Null
  }
  rclone listremotes
}
```

</div>

The last command must print <code><span data-rclone="remote">ovh</span>:</code> (and
your other rclone configurations, if any). The secret key is now stored in rclone's
configuration file, readable by you only: do not copy it anywhere else.

### Step 3 — Check access

``` { .text .rclone-tpl }
echo "test" | rclone rcat {{remote}}:{{bucket}}/access_test.txt
rclone lsf {{remote}}:{{bucket}}
rclone deletefile {{remote}}:{{bucket}}/access_test.txt
```

The second command must print `access_test.txt`. If a command prints
`AccessDenied`, check the keys (step 2), then contact whoever gave them to you.

### Step 4 — Practise with dummy data (mandatory)

This step is not optional. It checks that your machine, the network and the bucket
work together, without any risk for your real data.

<div data-only="linux mac" markdown>

**a) Create the generation script.** Paste this block: it creates the file
`generate_test_data.sh` in the current folder and makes it executable.

??? example "Show the block that creates the script (60 lines)"

    ``` bash
    cat > generate_test_data.sh <<'EOF'
    --8<-- "docs/guides/rclone/generate_test_data.sh"
    EOF
    chmod +x generate_test_data.sh
    ```

You can also [download the script](rclone/generate_test_data.sh) and read it first.

**b) Generate the data.** They go to <span data-rclone="exo">$HOME/rclone-exercise</span>.
The script refuses to start if the folder is not empty or if there is not enough space.

``` { .text .rclone-tpl }
mkdir -p "{{logs}}"
./generate_test_data.sh "{{exo}}" 1
```

</div>

<div data-only="win" markdown>

**a)** Nothing to create on Windows: generation is a single block, pasted at step b).

**b) Generate the data.** They go to <span data-rclone="exo">rclone-exercise</span>.
The block refuses to start if the folder is not empty or if there is not enough space.

??? example "Show the generation block (33 lines)"

    ``` { .text .rclone-tpl }
    & {
      $Dest = "{{exo}}"; $Vol = 1
      New-Item -ItemType Directory -Force -Path "{{logs}}" | Out-Null
      if ((Test-Path -LiteralPath $Dest) -and (Get-ChildItem -LiteralPath $Dest -Force | Select-Object -First 1)) {
        Write-Host "Error: $Dest exists and is not empty."; return
      }
      New-Item -ItemType Directory -Force -Path $Dest | Out-Null
      $Base = if ($Vol -lt 10) { 2 } else { $Vol }
      $Need = [math]::Floor($Base * 11 / 10)
      $Free = [math]::Floor((Get-Item -LiteralPath $Dest).PSDrive.Free / 1GB)
      if ($Free -lt $Need) {
        Write-Host "Error: $Need GB of free space needed, only $Free GB available on this disk."
        Remove-Item -LiteralPath $Dest; return
      }
      if ($Vol -lt 10) { $NM = 3; $SM = 100MB; $NL = 1; $SL = 600MB }
      else { $NM = [math]::Floor($Vol * 70 / 100); $SM = 1GB; $NL = [math]::Floor($Vol * 30 / 1000); $SL = 10GB }
      $rng = [Security.Cryptography.RandomNumberGenerator]::Create()
      function Write-Random([string]$Path, [long]$Bytes) {
        $buf = New-Object byte[] ([math]::Min($Bytes, 1MB))
        $fs = [IO.File]::Create($Path)
        try { $left = $Bytes; while ($left -gt 0) { $n = [int][math]::Min($left, $buf.Length); $rng.GetBytes($buf); $fs.Write($buf, 0, $n); $left -= $n } }
        finally { $fs.Dispose() }
      }
      foreach ($d in "run_01/small", "run_01/medium", "run_02/large", "folder with spaces/données accentuées") {
        New-Item -ItemType Directory -Force -Path (Join-Path $Dest $d) | Out-Null
      }
      for ($i = 1; $i -le 2000; $i++) { Write-Random (Join-Path $Dest ("run_01/small/read_{0:D4}.bin" -f $i)) 100KB }
      for ($i = 1; $i -le $NM; $i++) { Write-Random (Join-Path $Dest ("run_01/medium/sample_{0:D3}.fastq.gz" -f $i)) $SM }
      for ($i = 1; $i -le $NL; $i++) { Write-Random (Join-Path $Dest ("run_02/large/alignment_{0:D3}.bam" -f $i)) $SL }
      Set-Content -LiteralPath (Join-Path $Dest "folder with spaces/données accentuées/échantillon (copie) n°1.txt") -Value "witness file"
      $f = Get-ChildItem -LiteralPath $Dest -Recurse -File
      "Done: {0:N2} GB in {1} files" -f (($f | Measure-Object Length -Sum).Sum / 1GB), $f.Count
    }
    ```

</div>

<p class="rclone-hint-line">1 GB: about 2000 small files, three medium ones and one
large one, plus a file with spaces and accents in its name.</p>

<div data-only="linux" markdown>

**c) Upload the data.** First open a `screen` session, which keeps running even if you
log out:

``` { .text .rclone-tpl }
screen -S transfer
```

Then, inside this session:

``` { .text .rclone-tpl }
BW="{{bw}}"
rclone copy "{{exo}}" {{remote}}:{{bucket}}/EXERCISE \
  --fast-list --transfers 4 \
  --bwlimit "$BW" -P \
  --log-file "{{logs}}/exercise_copy.log" --log-level INFO || tail -n 5 "{{logs}}/exercise_copy.log"
```

To leave `screen` without stopping the transfer: ++ctrl+a++ then ++d++. To come back:
`screen -r transfer`.

</div>

<div data-only="mac" markdown>

**c) Upload the data.** `caffeinate -i` keeps the Mac awake while rclone runs; keep the
lid open.

``` { .text .rclone-tpl }
BW="{{bw}}"
caffeinate -i rclone copy "{{exo}}" {{remote}}:{{bucket}}/EXERCISE \
  --exclude .DS_Store --fast-list --transfers 4 \
  --bwlimit "$BW" -P \
  --log-file "{{logs}}/exercise_copy.log" --log-level INFO || tail -n 5 "{{logs}}/exercise_copy.log"
```

</div>

<div data-only="win" markdown>

**c) Upload the data.** Keep the PowerShell window open until the end.

``` { .text .rclone-tpl }
$BW = "{{bw}}"
rclone copy "{{exo}}" "{{remote}}:{{bucket}}/EXERCISE" --fast-list --transfers 4 --bwlimit $BW -P --log-file "{{logs}}\exercise_copy.log" --log-level INFO
if ($LASTEXITCODE -ne 0) { Get-Content "{{logs}}\exercise_copy.log" -Tail 5 }
```

</div>

- Progress is displayed live (option `-P`); the detailed log is
  <code><span data-rclone="logcopy">$HOME/rclone-logs/exercise_copy.log</span></code>.
- If rclone fails, the end of the command prints the last lines of this log, which say
  why (see Troubleshooting).
- The `BW` line limits the bandwidth to <span data-rclone="day">50M</span>iB/s on
  weekdays from 8 am to 8 pm, and to <span data-rclone="night">80M</span>iB/s at night
  and at weekends, so as not to saturate the site network. Times are those of the
  machine.

**d) Check the copy**, once the transfer is over:

<div data-only="linux" markdown>

``` { .text .rclone-tpl }
rclone check "{{exo}}" {{remote}}:{{bucket}}/EXERCISE --one-way --fast-list \
  --combined "{{logs}}/exercise_check.txt"
rclone size "{{exo}}"
rclone size {{remote}}:{{bucket}}/EXERCISE
```

</div>

<div data-only="mac" markdown>

``` { .text .rclone-tpl }
rclone check "{{exo}}" {{remote}}:{{bucket}}/EXERCISE --one-way --fast-list \
  --exclude .DS_Store --combined "{{logs}}/exercise_check.txt"
rclone size "{{exo}}" --exclude .DS_Store
rclone size {{remote}}:{{bucket}}/EXERCISE --exclude .DS_Store
```

</div>

<div data-only="win" markdown>

``` { .text .rclone-tpl }
rclone check "{{exo}}" "{{remote}}:{{bucket}}/EXERCISE" --one-way --fast-list --combined "{{logs}}\exercise_check.txt"
rclone size "{{exo}}"
rclone size "{{remote}}:{{bucket}}/EXERCISE"
```

</div>

The exercise succeeded if the check ends with `0 differences found` and if the last
two commands print the same number of objects and the same total size.
<span data-only="mac">`--exclude .DS_Store` leaves out the small files that Finder
creates in every folder you open: they are useless in the bucket.</span>

??? tip "Going further: measure the throughput"
    With 1 GB, the exercise checks that everything works, but says little about how long
    a large transfer will take. To find out, generate a bigger set: the last number of the
    generation command is the volume in GB<span data-only="win"> (`$Vol` at the start of
    the generation block)</span>. As a reference, 100 GB took 6 minutes to generate and 35
    minutes to upload from a lab server.<span data-only="linux"> On a server, put such a
    volume on a data disk, never in your home directory, which is often on the system
    disk: filling it up can freeze the server.</span> Then adapt the path in the copy,
    check and clean-up commands. The speed is displayed live by `-P`; for the average,
    divide the volume transferred by the elapsed time, both given in the last lines of
    the log.

<div data-only="mac win" markdown>

**e) Optional: the same upload with Cyberduck.** Cyberduck is a free graphical client
for Mac and Windows. This tests it on the same data, then checks with rclone that what
it sent is intact.

1. Install Cyberduck from [cyberduck.io/download](https://cyberduck.io/download/).
2. Click "Open Connection" and choose the protocol **S3 (HTTPS)** in the list at the top
   (if it is not there, enable it in the preferences, under Profiles). Server:
   <code><span data-rclone="endpointhost">s3.eu-west-par.io.cloud.ovh.net</span></code>,
   port 443, then the Access Key ID and the Secret Access Key, typed in this window
   only.
3. Open the bucket, create a folder named `CYBERDUCK-EXERCISE` (Action menu, then New
   Folder) and open it. Select everything inside the folder
   `rclone-exercise` of your home folder, in your file manager, and drag it into
   the Cyberduck window. Note how long the transfer takes.
4. Check with rclone. The `--download` option makes rclone compare the content byte by
   byte: for large files, which Cyberduck sends in several parts, the bucket keeps no
   checksum that rclone can compare.

<div data-only="mac" markdown>

``` { .text .rclone-tpl }
rclone check "{{exo}}" {{remote}}:{{bucket}}/CYBERDUCK-EXERCISE --one-way --download \
  --exclude .DS_Store --combined "{{logs}}/cyberduck_check.txt"
rclone size "{{exo}}" --exclude .DS_Store
rclone size {{remote}}:{{bucket}}/CYBERDUCK-EXERCISE --exclude .DS_Store
```

</div>

<div data-only="win" markdown>

``` { .text .rclone-tpl }
rclone check "{{exo}}" "{{remote}}:{{bucket}}/CYBERDUCK-EXERCISE" --one-way --download --combined "{{logs}}\cyberduck_check.txt"
rclone size "{{exo}}"
rclone size "{{remote}}:{{bucket}}/CYBERDUCK-EXERCISE"
```

</div>

The check must end with `0 differences found`, and the last two commands must print the
same number of objects and the same total size. Keep the folder `CYBERDUCK-EXERCISE` for
now: step 5 uses it.

</div>

### Step 5 — Practise managing files in the bucket (recommended)

Day-to-day work is rarely one big transfer: a new file arrives, a batch of files must be
stored, a file was put in the wrong place, another one must go. This step practises
these operations on the exercise data, which are still in the bucket under `EXERCISE`.

In a bucket, a folder is only a part of the files' paths. You never need to create a
folder before copying a file into it, and a folder disappears by itself when its last
file is moved or deleted.

**a) Prepare the new files**: a 5 MB file, `new_sample.bam`, and a folder of ten 1 MB
files, `new_batch`, both in your home folder.

<div data-only="linux mac" markdown>

``` { .text .rclone-tpl }
head -c 5242880 /dev/urandom > "{{newfile}}"
mkdir -p "{{newbatch}}"
for i in $(seq -w 1 10); do head -c 1048576 /dev/urandom > "{{newbatch}}/file_$i.bin"; done
ls -l "{{newfile}}" "{{newbatch}}"
```

</div>

<div data-only="win" markdown>

``` { .text .rclone-tpl }
& {
  $rng = [Security.Cryptography.RandomNumberGenerator]::Create()
  function Write-Random([string]$Path, [int]$Bytes) { $b = New-Object byte[] $Bytes; $rng.GetBytes($b); [IO.File]::WriteAllBytes($Path, $b) }
  Write-Random "{{newfile}}" 5MB
  New-Item -ItemType Directory -Force -Path "{{newbatch}}" | Out-Null
  for ($i = 1; $i -le 10; $i++) { Write-Random (Join-Path "{{newbatch}}" ("file_{0:D2}.bin" -f $i)) 1MB }
  Get-ChildItem -LiteralPath "{{newfile}}", "{{newbatch}}" | Select-Object Name, Length
}
```

</div>

**b) Upload a new file**, first into an existing folder:

``` { .text .rclone-tpl }
rclone copy "{{newfile}}" {{remote}}:{{bucket}}/EXERCISE/run_02/large
rclone lsf {{remote}}:{{bucket}}/EXERCISE/run_02/large
```

`rclone copy` takes a destination folder. The `lsf` must list `new_sample.bam` next to
the `alignment_…bam` files already there.

Then into a folder that does not exist yet:

``` { .text .rclone-tpl }
rclone copyto "{{newfile}}" {{remote}}:{{bucket}}/EXERCISE/run_03/new_sample.bam
rclone lsf {{remote}}:{{bucket}}/EXERCISE/run_03
```

`rclone copyto` takes the full name of the file at destination, which also lets you
rename it on the way. The `lsf` must print `new_sample.bam`: the folder `run_03` did not
exist, the copy created it.

**c) Upload a new folder**, then check it:

``` { .text .rclone-tpl }
rclone copy "{{newbatch}}" {{remote}}:{{bucket}}/EXERCISE/new_batch
rclone check "{{newbatch}}" {{remote}}:{{bucket}}/EXERCISE/new_batch --one-way
```

!!! warning "Name the folder in the destination"
    `rclone copy` copies the *content* of a folder. Without `/new_batch` at the end of
    the destination, the ten files would land loose in `EXERCISE`, mixed with the rest.

**d) Move and rename.** The files are moved inside the bucket, without being uploaded
again:

``` { .text .rclone-tpl }
rclone moveto {{remote}}:{{bucket}}/EXERCISE/run_03/new_sample.bam {{remote}}:{{bucket}}/EXERCISE/run_01/sample_renamed.bam
rclone move {{remote}}:{{bucket}}/EXERCISE/new_batch {{remote}}:{{bucket}}/EXERCISE/run_02/new_batch
rclone lsf {{remote}}:{{bucket}}/EXERCISE
rclone lsf {{remote}}:{{bucket}}/EXERCISE/run_02
```

`moveto` moves one file and can rename it; `move` moves the content of a folder, with
the same rule as `copy`: name the folder in the destination. The first `lsf` no longer
lists `run_03/` nor `new_batch/`: once empty, they are gone. The second lists
`new_batch/` in its new place.

**e) Delete.** Always run the command with `--dry-run` first: it deletes nothing and
only shows what would be deleted.

``` { .text .rclone-tpl }
rclone deletefile {{remote}}:{{bucket}}/EXERCISE/run_01/sample_renamed.bam --dry-run
```

The output must name `sample_renamed.bam` only, followed by
`Skipped delete as --dry-run is set`. Then, and only then, delete it for real:

``` { .text .rclone-tpl }
rclone deletefile {{remote}}:{{bucket}}/EXERCISE/run_01/sample_renamed.bam
```

`deletefile` deletes exactly one file and refuses a folder: it is the safest way to
delete. To delete a whole folder, use `rclone purge`:

!!! danger "`rclone purge` permanently deletes everything under the given path, without confirmation and without a trash"
    A mistake in the path can wipe out the whole bucket, with all the data it holds.
    Read the path three times before pressing Enter, and run the command with
    `--dry-run` first.

``` { .text .rclone-tpl }
rclone purge {{remote}}:{{bucket}}/EXERCISE/run_02/new_batch --dry-run
```

The output must end with `EXERCISE/run_02/new_batch: Skipped purge directory as
--dry-run is set`: check that this path is the folder to delete, and nothing above it.
Then delete it for real, and check that only `large/` remains in `run_02`:

``` { .text .rclone-tpl }
rclone purge {{remote}}:{{bucket}}/EXERCISE/run_02/new_batch
rclone lsf {{remote}}:{{bucket}}/EXERCISE/run_02
```

<div data-only="mac win" markdown>

**f) Optional: the same operations with Cyberduck**, in the folder `CYBERDUCK-EXERCISE`:

1. Drag `new_sample.bam` from your home folder into `run_01`.
2. Rename it: select it, press Return, type `renamed.bam`, press Return again.
3. Move it: Edit menu, Cut; open `run_02`; Edit menu, Paste ("Édition", "Couper",
   "Coller" if Cyberduck is in French).
4. Drag the folder `new_batch` from your home folder into `CYBERDUCK-EXERCISE`.
5. Look at the result with rclone:

    ``` { .text .rclone-tpl }
    rclone lsf {{remote}}:{{bucket}}/CYBERDUCK-EXERCISE/run_01
    rclone lsf {{remote}}:{{bucket}}/CYBERDUCK-EXERCISE/run_02
    rclone lsf {{remote}}:{{bucket}}/CYBERDUCK-EXERCISE/new_batch
    ```

    `run_02` must list `renamed.bam` next to `large/`, and `new_batch` its ten files.
    With an OVHcloud bucket, Cyberduck can leave the original file in place after a
    rename or a move: if `run_01` still lists `new_sample.bam` or `renamed.bam`, delete
    them in Cyberduck. rclone does not have this problem (step d).

6. Delete `renamed.bam`: select it, then Delete ("Effacer" in French), and confirm.

</div>

**g) Clean up** the bucket<span data-only="mac win"> (both exercise folders)</span>
and your computer. Same care as above: the path must end with `/EXERCISE`.

``` { .text .rclone-tpl }
rclone purge {{remote}}:{{bucket}}/EXERCISE
```

<div data-only="mac win" markdown>

In Cyberduck, select the folder `CYBERDUCK-EXERCISE` at the root of the bucket, choose
Delete ("Effacer" in French) and confirm: Cyberduck removes the folder and everything
in it. Do it in Cyberduck rather than with `rclone purge`, which would leave the empty
folder created by Cyberduck behind.

</div>

<div data-only="linux mac" markdown>

``` { .text .rclone-tpl }
rm -rf "{{exo}}" "{{newfile}}" "{{newbatch}}"
```

</div>

<div data-only="win" markdown>

``` { .text .rclone-tpl }
Remove-Item -LiteralPath "{{exo}}", "{{newfile}}", "{{newbatch}}" -Recurse -Force
```

</div>

Then send the files of the
<code><span data-rclone="logs">$HOME/rclone-logs</span></code> directory to whoever
supports you.

### Step 6 — Transfer the real data

The same commands, with the directory of the real data
(<span data-rclone="src">/path/to/data</span>) and the destination directory in the
bucket (<span data-rclone="dst">my-dataset</span>). First look for symbolic links,
which rclone ignores by default: if the following command prints anything, ask for
advice before going further.

<div data-only="linux mac" markdown>

``` { .text .rclone-tpl }
find "{{src}}" -type l | head
```

</div>

<div data-only="win" markdown>

``` { .text .rclone-tpl }
Get-ChildItem -LiteralPath "{{src}}" -Recurse -Attributes ReparsePoint | Select-Object -First 10 -ExpandProperty FullName
```

</div>

<div data-only="linux" markdown>

Open a `screen` session, then start the transfer inside it:

``` { .text .rclone-tpl }
screen -S transfer
```

``` { .text .rclone-tpl }
BW="{{bw}}"
rclone copy "{{src}}" "{{remote}}:{{bucket}}/{{dst}}" \
  --fast-list --transfers 4 \
  --bwlimit "$BW" -P \
  --log-file "{{logs}}/transfer_copy.log" --log-level INFO || tail -n 5 "{{logs}}/transfer_copy.log"
```

</div>

<div data-only="mac" markdown>

``` { .text .rclone-tpl }
BW="{{bw}}"
caffeinate -i rclone copy "{{src}}" "{{remote}}:{{bucket}}/{{dst}}" \
  --exclude .DS_Store --fast-list --transfers 4 \
  --bwlimit "$BW" -P \
  --log-file "{{logs}}/transfer_copy.log" --log-level INFO || tail -n 5 "{{logs}}/transfer_copy.log"
```

</div>

<div data-only="win" markdown>

``` { .text .rclone-tpl }
$BW = "{{bw}}"
rclone copy "{{src}}" "{{remote}}:{{bucket}}/{{dst}}" --fast-list --transfers 4 --bwlimit $BW -P --log-file "{{logs}}\transfer_copy.log" --log-level INFO
if ($LASTEXITCODE -ne 0) { Get-Content "{{logs}}\transfer_copy.log" -Tail 5 }
```

</div>

Once the transfer is over, check:

<div data-only="linux" markdown>

``` { .text .rclone-tpl }
rclone check "{{src}}" "{{remote}}:{{bucket}}/{{dst}}" --one-way --fast-list \
  --combined "{{logs}}/transfer_check.txt"
rclone size "{{src}}"
rclone size "{{remote}}:{{bucket}}/{{dst}}"
```

</div>

<div data-only="mac" markdown>

``` { .text .rclone-tpl }
rclone check "{{src}}" "{{remote}}:{{bucket}}/{{dst}}" --one-way --fast-list \
  --exclude .DS_Store --combined "{{logs}}/transfer_check.txt"
rclone size "{{src}}" --exclude .DS_Store
rclone size "{{remote}}:{{bucket}}/{{dst}}" --exclude .DS_Store
```

</div>

<div data-only="win" markdown>

``` { .text .rclone-tpl }
rclone check "{{src}}" "{{remote}}:{{bucket}}/{{dst}}" --one-way --fast-list --combined "{{logs}}\transfer_check.txt"
rclone size "{{src}}"
rclone size "{{remote}}:{{bucket}}/{{dst}}"
```

</div>

If the transfer is interrupted, run exactly the same command again: rclone only sends
what is missing.

!!! warning "Always `rclone copy`, never `rclone sync`"
    `sync` deletes from the bucket whatever no longer exists on your machine.

### Troubleshooting

- **The copy command returns at once, or ends with `Errors: 1` or more**: read the last
  lines of the log, which the command prints in that case. `didn't find section in
  config file` means that the rclone configuration name in "Your values" is not the one
  created at step 2 (the command `rclone listremotes` lists the existing ones).
  `directory not found` means that the source directory does not exist on this
  machine: check its path in "Your values".
- **`AccessDenied`**: the keys were mistyped at step 2. Delete the configuration with
  <code>rclone config delete <span data-rclone="remote">ovh</span></code>, then paste
  the step 2 block again.
- **`rclone` is not found** after installation: open a new terminal window and try
  again.
- **The site network or the machine slows down**: stop rclone (++ctrl+c++), lower both
  bandwidths in "Your values", copy the command again and restart it. Once all
  transfers are over, this command deletes interrupted uploads older than 24 hours,
  which would otherwise still be billed in the bucket:

    ``` { .text .rclone-tpl }
    rclone backend cleanup {{remote}}:{{bucket}}
    ```

- **The check finds differences**: run the `rclone copy` command again with
  `--checksum` added, then check again.

## Part 2 — Cold archive (optional)

This part only concerns data that have already been transferred and checked, and that
you decide to keep for a long time without reading them. The decision belongs to the
people responsible for the data.

OVHcloud offers its cheapest storage class, Cold Archive, only in the Paris 3-AZ region
(`eu-west-par`): the bucket must have been created there. Three constraints to know
before deciding:

- **180 days minimum**: a file archived, then deleted or replaced before 180 days is
  billed until the end of that period.
- **Delayed reading**: an archived file must be restored before it can be read; the
  restoration is billed and can take up to 48 hours.
- **File-by-file restoration**: it is better to bundle small files into tar archives
  (for instance one per sequencing run) before the transfer.

Two ways to do it:

- **Lifecycle rule (recommended)**: in the OVHcloud Control Panel, *Lifecycle* tab of
  the bucket, create a rule that transitions objects to Cold Archive after N days.
  Files move automatically; the transition itself is free.
- **Direct upload to the archive**: add `--s3-storage-class DEEP_ARCHIVE` to the
  `rclone copy` command. Simpler, but leaves no room for mistakes.

To restore an archived directory and keep it readable for 7 days: the restoration does
not touch the archived version, which stays in Cold Archive. It is not a swap. It
creates a temporary readable copy alongside, billed for its lifetime, which disappears
by itself after the 7 days.

``` { .text .rclone-tpl }
rclone backend restore "{{remote}}:{{bucket}}/{{dst}}" -o priority=Standard -o lifetime=7
```

See OVHcloud's [Cold Archive FAQ](https://docs.ovhcloud.com/en/guides/storage-and-backup/object-storage/cold-archive-faq)
for current conditions and prices.

## Appendix — The generation script (Linux and Mac)

The same script as in step 4 a), without the command that creates the file.

??? example "Show the script"

    ``` bash
    --8<-- "docs/guides/rclone/generate_test_data.sh"
    ```
