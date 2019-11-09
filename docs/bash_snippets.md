### diff
##### differences of 2 stdout outputs
```
diff <(egrep --color "<regexp>" <pathtofile1>) <(egrep --color "<regexp>" <pathtofile2>)
```
### egrep
```
egrep --color "<regexp>" <pathtofile>
```
### Remove lines before and after some strings in a file
##### remove before somestring, retaining somestring
`awk '/somestring/,0' <pathtofile> > <pathtodeleted_file>`
##### remove after someotherstring, including someotherstring
`sed -n '/someotherstring/q;p' <pathtofile > <pathtodeleted_file>`
### Check which ports are busy and which ports are free
```
sudo netstat -tulpn
sudo netstat -antup
sudo lsof -i -n -P
# You can verify process using port /proc:
ls -l /proc/<PID>/exe
```
### List Open Files For Process
First you need to find out PID of process. Simply use any one of the following command to obtain process id:
`ps aux | grep {program-name}` OR `ps -C {program-name} -o pid=`

To list open files for firefox process, enter:
`ls -l /proc/<PID>/fd`

##### lsof 
lsof command list open files under all Linux distributions or UNIX like operating system.

Type the following command to list open file for process ID 351:
`lsof -p 351`
### Empty the content of a file without closing it
```
truncate -s 0 filename
```
### Calculate total used disk space by files older than 180 days using find
From [stackoverflow](http://stackoverflow.com/questions/17419337/calculate-total-used-disk-space-by-files-older-than-180-days-using-find),
here an example for space occupied by files older than 5 years (1825 days)
```
find . -type f -mtime +1825 -printf '%s\n' | awk '{a+=$1;} END {printf "%.1f GB\n", a/2**30;}'
```
### generate a private/public ssh key pair and share it for ssh connection



