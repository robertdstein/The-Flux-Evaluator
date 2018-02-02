"""Script to run stacking scripts on the DESY cluster.

Through use of argparse, a given configuration for the code can be selected.
This can be given from the command line, in the form:

python RunCluster.py -c Desired_Configuration_Name -n Number_Of_Tasks -s

Each available configuration must be listed in "config.ini", and controls
options for fitting, such as which catalogue is to be used, and which seasons
of data should be included. If -x is included, then a new job is submitted
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

            n_total = len(tmp.split('\n')) - 3

            running_process = subprocess.Popen(
                cmd + " -s r", stdout=subprocess.PIPE, shell=True)
            running_tmp = str(running_process.stdout.read())

            if running_tmp != "":
                n_running = len(running_tmp.split('\n')) - 3
            else:
                n_running = 0

            print time.asctime(time.localtime()), n_total, "entries in queue. ",
            print "Of these,", n_running, "are running tasks, and",
            print n_total-n_running, "are jobs still waiting to be executed."
            print time.asctime(time.localtime()), "Waiting for Cluster"
            i = 0
            j += 1
        # if j > 5:
        #     print tmp
        #     j = 0
        time.sleep(30)
        i += 1
        process = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True)
        tmp = str(process.stdout.read())


def submit_to_cluster(
        tasks, config, conf_file="analysis_config/config.ini", ntrials=50,
        steps=15, sh_file="SubmitOne.sh"):
    # Submits job to the cluster, with a command in the form of:
    # qsub -t 1-50:1 SubmitOne.sh Full_with_DaiFang_TDE
    submit_cmd = "qsub -t 1-" + str(tasks) + ":1 " + sh_file + " " + config \
                 + " " + conf_file + " " + str(ntrials) + " " + str(steps)
    print time.asctime(time.localtime()), submit_cmd, "\n"
    os.system(submit_cmd)

# If script is run from command line, automatically uses run()
# If imported into another script, run() must be explicitly called
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-t", "--tasks", default=50,
                        help="Number of tasks to run")
    parser.add_argument("-x", "--execute", action="store_true",
                        help="Toggle to submit job to cluster")
    parser.add_argument("-c", "--config", default="Full_with_DaiFang_TDE",
                        help="Sets configuration for LLh minimisation")
    parser.add_argument("-f", "--conf_file", default="config.ini")
    parser.add_argument("-n", "--ntrials", default=15,
                        help="Number of trials per flux step per task")
    parser.add_argument("-s", "--step", default=15,
                        help="Number of flux steps per task")
    cfg = parser.parse_args()

    # If the job should be submitted
    if cfg.execute:
        import ConfigParser

        # Sets root directory for config.ini file
        root = "/afs/ifh.de/user/s/steinrob/Desktop/python/The-Flux-Evaluator/"

        conf = ConfigParser.ConfigParser()
        config_path = root + "analysis_config/" + cfg.conf_file
        conf.read(config_path)

        if cfg.config not in conf.sections():
            print "Searching for config section", cfg.config,
            print "in", conf.sections()
            raise Exception("Config file not found.")
        else:
            submit_to_cluster(cfg.tasks, cfg.config, cfg.conf_file,
                              cfg.ntrials, cfg.step)


    # In any case, check for the cluster
    wait_for_cluster()

    Mf.run(cfg.config)
