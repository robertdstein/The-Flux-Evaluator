import os
import sys

from common import root_path, tfe_path, input_path, storage_path, \
    output_path, log_path
import prepare_catalogue

print "\n \n"
print "************************************************************************"
print "*                                                                      *"
print "*             Initialising setup for The Flux Evaluator                *"
print "*                                                                      *"
print "************************************************************************"
print "\n"
print "Initialising directory for data storage. This could be a scratch space" \
      " or local directory."
print "\n"

print "settings.py has found the following source_path directory: \n"
print root_path
print
print "Is this correct? (y/n)"

x = ""

while x not in ["y", "n"]:
    x = raw_input("")

if x == "n":
    print "\n"
    print "Please edit settings.py to include the correct directory!"
    sys.exit()

for path in [tfe_path, input_path, storage_path, output_path, log_path]:
    if not os.path.isdir(path):
        print "Making results_path", path
        os.makedirs(path)
    else:
        print "Found directory", path

print "\n"
print "************************************************************************"
print "*                                                                      *"
print "*                   Initialising catalogue creation                    *"
print "*                                                                      *"
print "************************************************************************"
print "\n"
prepare_catalogue.make_single_sources()
