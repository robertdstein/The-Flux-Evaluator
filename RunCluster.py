"""Script to run stacking scripts on the DESY cluster.

Through use of argparse, a given configuration for the code can be selected.
This can be given from the command line, in the form:

python RunCluster.py -c Desired_Configuration_Name -n Number_Of_Tasks -s

Each available configuration must be listed in "config.ini", and controls
options for fitting, such as which catalogue is to be used, and which seasons
of data should be included. If -s is included, then a new job is submitted
to the cluster. Having submitted the job to the cluster it will be run in
parallel Number_of_Tasks times. The shell script SubmitOne.sh is called for
each task, which in turn calls RunLocal.py with the given configuration setting.

The Wait function will periodically query the cluster
to check on the status of the job, and will output the job status occasionally.

Once all sub-tasks are completed, the script will proceed to call
MergeFiles.run() for the given configuration, combining results.

"""
import subprocess
import time
import os.path
import argparse
import MergeFiles as Mf

parser = argparse.ArgumentParser()
parser.add_argument("-n", "--n", default=50, help="Number of tasks to run")
parser.add_argument("-s", "--submit", action="store_true",
                    help="Toggle to submit job to cluster")
parser.add_argument("-c", "--config", default="Full_with_DaiFang_TDE",
                    help="Sets configuration for LLh minimisation")
cfg = parser.parse_args()

cmd = 'qstat -u steinrob'


def wait_for_cluster():
    """Runs the command cmd, which queries the status of the job on the
    cluster, and reads the output. While the output is not an empty
    string (indicating job completion), the cluster is re-queried
    every 30 seconds. Occasionally outputs the number of remaining sub-tasks
    on cluster, and outputs full table result every ~ 8 minutes. On
    completion of job, terminates function process and allows the script to
    continue.
    """
    process = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True)
    tmp = str(process.stdout.read())
    i = 31
    j = 6
    while tmp != "":
        if i > 3:
            print time.asctime(time.localtime()), len(tmp.split('\n')) - 3, \
                "entries in queue"
            print time.asctime(time.localtime()), "Waiting for Cluster"
            i = 0
            j += 1
        if j > 5:
            print tmp
            j = 0
        time.sleep(30)
        i += 1
        process = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True)
        tmp = str(process.stdout.read())


# If the job should be submitted
if cfg.submit:
    import ConfigParser

    # Sets root directory for config.ini file
    root = "/afs/ifh.de/user/s/steinrob/Desktop/python/stacking/"

    conf = ConfigParser.ConfigParser()
    conf.read(root + "config.ini")

    if cfg.config not in conf.sections():
        print "Searching for config section", cfg.config, "in", conf.sections()
        raise Exception("Config file not found.")

    else:

        # Submits job to the cluster, with a command in the form of:
        # qsub -t 1-50:1 SubmitOne.sh Full_with_DaiFang_TDE
        submit_cmd = "qsub -t 1-" + str(cfg.n) + ":1 SubmitOne.sh " + cfg.config
        print time.asctime(time.localtime()), submit_cmd, "\n"
        os.system(submit_cmd)

# In any case, check for the cluster
wait_for_cluster()

Mf.run(cfg.config)
