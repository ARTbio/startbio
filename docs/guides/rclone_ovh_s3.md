# Transfer data to an OVH S3 bucket with rclone

This guide shows how to copy data from a Linux server to an OVHcloud Object Storage
bucket (S3) with [rclone](https://rclone.org/), and how to check that everything
arrived intact. You always start with an exercise: you transfer dummy data, generated
for the occasion, with exactly the same commands as for your real data.

Moving the data to cold archive afterwards is optional and covered separately at
the end of the page ([Part 2](#part-2-cold-archive-optional)).

## Your values

<div id="rclone-form" class="rclone-form" markdown>

Fill in these fields: every command on this page updates as you type, and you only
have to copy it. Fields highlighted in yellow still hold an example value. Your values
are kept in this browser for your next visit; nothing is sent anywhere.

<div class="rclone-fields"></div>

<button type="button" class="md-button rclone-reset">Reset to example values</button>

!!! warning "Never type your access keys here"
    Keys are typed only in the server's terminal, at step 2.

Preview: the real transfer command, with your values (step 5).

``` { .text .rclone-tpl }
BW="{{bw}}"
rclone copy "{{src}}" "{{remote}}:{{bucket}}/{{dst}}" \
  --fast-list --transfers 4 \
  --bwlimit "$BW" -P \
  --log-file "{{logs}}/transfer_copy.log" --log-level INFO
```

</div>

## Before you start

- An account on the Linux server that holds the data, with Internet access.
- Three pieces of information from whoever manages your OVHcloud Public Cloud project:
  the bucket name (<span data-rclone="bucket">my-bucket</span>), an access key and a
  secret key.
- At least <span data-rclone="need">110</span> GB of free space on a data disk of the
  server for the exercise, not on the system disk.
- About one hour for the 100 GB exercise. Measured on a lab server: 6 minutes to
  generate the data and 35 minutes to upload them, plus the time to check, clean up
  and read the commands.

## Part 1 — Transfer the data

### Step 1 — Install rclone

In a terminal on the server (no administrator rights needed). These commands are for
a standard x86-64 server (`uname -m` prints `x86_64`).

``` { .text .rclone-tpl }
mkdir -p ~/bin && cd /tmp
curl -O https://downloads.rclone.org/rclone-current-linux-amd64.zip
unzip -o rclone-current-linux-amd64.zip
cp rclone-*-linux-amd64/rclone ~/bin/
cd -
echo 'export PATH="$HOME/bin:$PATH"' >> ~/.bashrc && source ~/.bashrc
hash -r && which rclone && rclone version
```

The output must show a path in your home directory ending in `/bin/rclone` (not
`/usr/bin/rclone`), then a version number.

### Step 2 — Configure access to the bucket

Paste this block into the terminal. It asks for the access key, then for the secret
key (which is not displayed while you type). Other rclone configurations already
present on the server are left untouched.

``` { .text .rclone-tpl }
read -rp "Paste the access key, then Enter: " AK
read -rsp "Paste the secret key (nothing is displayed, this is normal), then Enter: " SK; echo
rclone config create {{remote}} s3 provider=Other env_auth=false \
  access_key_id="$AK" secret_access_key="$SK" \
  region={{region}} endpoint=https://s3.{{region}}.io.cloud.ovh.net \
  acl=private no_check_bucket=true > /dev/null
unset AK SK
chmod 600 "$(rclone config file | tail -1)"
rclone listremotes
```

The last command must print <code><span data-rclone="remote">ovh</span>:</code> (and
your other rclone configurations, if any). The secret key is stored in a file that
only you can read: do not copy it anywhere else.

### Step 3 — Check access

``` { .text .rclone-tpl }
echo "test" | rclone rcat {{remote}}:{{bucket}}/access_test.txt
rclone lsf {{remote}}:{{bucket}}
rclone deletefile {{remote}}:{{bucket}}/access_test.txt
```

The second command must print `access_test.txt`. If a command prints
`AccessDenied`, check the keys (step 2), then contact whoever gave them to you.

### Step 4 — Practise with dummy data (mandatory)

This step is not optional. It checks that the server, the network and the bucket work
together, and it tells you how long the real transfer will take, without any risk for
your real data.

**a) Create the generation script.** Paste this block into the terminal: it creates
the file `generate_test_data.sh` and makes it executable.

??? example "Show the block that creates the script (60 lines)"

    ``` bash
    cat > generate_test_data.sh <<'EOF'
    --8<-- "docs/guides/rclone/generate_test_data.sh"
    EOF
    chmod +x generate_test_data.sh
    ```

You can also [download the script](rclone/generate_test_data.sh) and read it first.

**b) Generate the data.** They go to <span data-rclone="exo">/data_disk/rclone-exercise</span>,
on the disk set in "Your values". Never use your home directory, which is often on the
system disk: filling it up can freeze the server. The script refuses to start if
there is not enough space. If you get "Permission denied", ask the server's
administrator for a directory in your name on that disk.

``` { .text .rclone-tpl }
mkdir -p "{{logs}}"
./generate_test_data.sh "{{exo}}" {{vol}}
```

Generation takes about 6 minutes for 100 GB on a typical lab server; it depends on the
server's disks. For a more thorough exercise, set "Exercise volume" to 500.

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
  --log-file "{{logs}}/exercise_copy.log" --log-level INFO
```

- To leave `screen` without stopping the transfer: ++ctrl+a++ then ++d++. To come
  back: `screen -r transfer`.
- Progress is displayed live (option `-P`); the detailed log is in
  <code><span data-rclone="logs">$HOME/rclone-logs</span>/exercise_copy.log</code>.
- The `BW` line limits the bandwidth to <span data-rclone="day">50M</span>iB/s on
  weekdays from 8 am to 8 pm, and to <span data-rclone="night">80M</span>iB/s at night
  and at weekends, so as not to saturate the site network or the server. Times are the
  server's local time.

**d) Check the copy**, once the transfer is over:

``` { .text .rclone-tpl }
rclone check "{{exo}}" {{remote}}:{{bucket}}/EXERCISE --one-way --fast-list \
  --combined "{{logs}}/exercise_check.txt"
rclone size "{{exo}}"
rclone size {{remote}}:{{bucket}}/EXERCISE
```

The exercise succeeded if the check ends with `0 differences found` and if the last
two commands print the same number of objects and the same total size.

**e) Clean up**, then send the files of the
<code><span data-rclone="logs">$HOME/rclone-logs</span></code> directory to whoever
supports you.

!!! danger "`rclone purge` permanently deletes everything under the given path, without confirmation and without a trash"
    Without the final `/EXERCISE`, the whole bucket disappears, with all the data it
    holds. Read the path three times before pressing Enter. If in doubt, first run the
    same command followed by `--dry-run`: it shows what would be deleted without
    deleting anything.

``` { .text .rclone-tpl }
rclone purge {{remote}}:{{bucket}}/EXERCISE
rm -rf "{{exo}}"
```

### Step 5 — Transfer the real data

The same commands, with the directory of the real data
(<span data-rclone="src">/path/to/data</span>) and the destination directory in the
bucket (<span data-rclone="dst">my-dataset</span>). First look for symbolic links,
which rclone ignores by default: if the following command prints anything, ask for
advice before going further.

``` { .text .rclone-tpl }
find "{{src}}" -type l | head
```

Open a `screen` session, then start the transfer inside it:

``` { .text .rclone-tpl }
screen -S transfer
```

``` { .text .rclone-tpl }
BW="{{bw}}"
rclone copy "{{src}}" "{{remote}}:{{bucket}}/{{dst}}" \
  --fast-list --transfers 4 \
  --bwlimit "$BW" -P \
  --log-file "{{logs}}/transfer_copy.log" --log-level INFO
```

Once the transfer is over, check:

``` { .text .rclone-tpl }
rclone check "{{src}}" "{{remote}}:{{bucket}}/{{dst}}" --one-way --fast-list \
  --combined "{{logs}}/transfer_check.txt"
rclone size "{{src}}"
rclone size "{{remote}}:{{bucket}}/{{dst}}"
```

If the transfer is interrupted, run exactly the same command again: rclone only sends
what is missing.

!!! warning "Always `rclone copy`, never `rclone sync`"
    `sync` deletes from the bucket whatever no longer exists on the server.

### Troubleshooting

- **`AccessDenied`**: the keys were mistyped in the configuration (step 2).
- **The site network or the server slows down**: stop rclone (++ctrl+c++ in
  `screen`), lower both bandwidths in "Your values", copy the command again and
  restart it. Once all transfers are over, this command deletes interrupted uploads
  older than 24 hours, which would otherwise still be billed in the bucket:

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

## Appendix — The generation script

The same script as in step 4 a), without the command that creates the file.

??? example "Show the script"

    ``` bash
    --8<-- "docs/guides/rclone/generate_test_data.sh"
    ```
