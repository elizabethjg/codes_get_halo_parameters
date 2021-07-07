#name of the executable script
executable = get_halo.sh

# Requests 4 cores into 1 machine
machine_count = 1
request_cpus = 4

#HTCondor output files
output = condor-$(ClusterId).$(ProcId).out
error = condor-$(ClusterId).$(ProcId).err
log = condor-$(ClusterId).$(ProcId).log

# Allows only 3 hours executions
+flavour="short"

queue
