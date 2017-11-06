#!/bin/zsh
##
##(otherwise the default shell would be used)
#$ -S /bin/zsh
##
##(the running time for this job)

#$ -l h_cpu=00:29:00
#$ -l h_rss=8G
##
##Force OS = SL5/6
#$ -l os=sl6
##
##(send mail on job's end and abort)
#$ -m a
##
##(stderr and stdout are merged together to stdout)
#$ -j y
##

## name of the job
## -N TDE Stacking Analysis
##
##(redirect output to:)
#$ -o /dev/null
##

sleep $(( ( RANDOM % 60 )  + 1 ))

exec > "$TMPDIR"/${JOB_ID}_stdout.txt 2>"$TMPDIR"/${JOB_ID}_stderr.txt

cp -R /afs/ifh.de/user/a/astasik/scratch/PS_Data/FinalSample $TMPDIR
cp -R /afs/ifh.de/user/a/astasik/scratch/PS_Data/Catalog $TMPDIR
cp -R /afs/ifh.de/user/a/astasik/scratch/PS_Data/DeclinationAcceptance $TMPDIR

eval $(/cvmfs/icecube.opensciencegrid.org/py2-v2/setup.sh)

$SROOT/metaprojects/offline-software/V16-10-00/env-shell.sh python /afs/ifh.de/user/s/steinrob/Desktop/python/stacking/RunLocal.py -c $1 -f $2 -n $3 -s $4

cp $TMPDIR/${JOB_ID}_stdout.txt /afs/ifh.de/user/s/steinrob/Desktop/python/stacking/logs
cp $TMPDIR/${JOB_ID}_stderr.txt /afs/ifh.de/user/s/steinrob/Desktop/python/stacking/logs
