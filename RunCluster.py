import os
import subprocess
import time
import os.path

cmd = 'qstat -u steinrob'

def wait():
	proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True)
	tmp = str(proc.stdout.read())
	i = 31
	j = 0
	while tmp != "":
		if i > 3:
			print time.asctime(time.localtime()), len(tmp.split('\n')) - 3, "entries in queue"
			print time.asctime(time.localtime()), "Waiting for Cluster"
			i = 1
			j += 1
		if j > 0:
			print tmp
			j = 0
		time.sleep(30)
		i += 1
		proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True)
		tmp = str(proc.stdout.read())

cluster_command = "qsub -t 1-50:1 SubmitOne.sh"
print time.asctime(time.localtime()), cluster_command, "\n"
os.system(cluster_command)

wait()