export LD_LIBRARY_PATH=/usr/lib64:$LD_LIBRARY_PATH
echo "slurm:x:151:151:Slurm:/:/sbin/nologin" >> /etc/passwd
echo "slurm:x:151:" >> /etc/group
