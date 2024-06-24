# test.sh
#!/bin/sh
# PBS -N test
# PBS -l nodes=master_vir1+master_vir2

pssh -h $PBS_NODEFILE mkdir -p /home/s2211995 1>&2
scp master:/home/s2113619/test /home/s2211995
pscp -h $PBS_NODEFILE master:/home/s2211995/test /home/s2211995 1>&2
mpiexec -np4 -machinefile $PBS_NODEFILE /home/s2211995/test
