```
# manage slurm job manager
watch "squeue --format '%.18i %.9P %.40j %.8u %.2t %.10M %.4C %.6D %R' && sinfo"
scontrol update NodeName=ansible-upgrade-2 State=Resume

# see galaxy server logs
tail -f /home/galaxy/galaxy/uwsgi.log



# In case of conda issue with libssl / openssh...
cd /lib/x86_64-linux-gnu
ln -s libssl.so.1.1 libssl.so.1.0.0
ln -s libcrypto.so.1.1 libcrypto.so.1.0.0

# In case of problems with some conda packages (only for conda geeks)
conda install -c bioconda samtools=1.9 --force-reinstall
conda install -c bioconda --force-reinstall ucsc-genepredtobed ucsc-gtftogenepred
conda clean --dry-run --all

# In case of global issue with conda tool depencies
tar cf - _conda | pigz -9 -p 8 > _conda.pigz.tar.gz
rm -rf _conda.pigz.tar.gz && tar -I pigz -xf _conda.pigz.tar.gz 
```
