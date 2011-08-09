"""This code looks in the specified directory and determines which types of
files are located there."""

import sys
import os
import string

neo_flag = False
gyro_flag = False
tgyro_flag = False
profiles_gen_flag = False
dirtree = os.walk(sys.argv[1])

for root, dirs, files in dirtree:
    for f in files:
        if string.find(f, 'out.neo') >= 0:
            neo_flag = True
        if string.find(f, 'out.tgyro') >= 0:
            tgyro_flag = True
        if string.find(f, 'out.gyro') >= 0:
            gyro_flag = True
        if string.find(f, 'input.profiles') >= 0:
            profiles_gen_flag = True

if neo_flag or gyro_flag or tgyro_flag or profiles_gen_flag:
    print "It looks like this directory contains the following file types:"
    if neo_flag:
        print "neo"
    if gyro_flag:
        print "gyro"
    if tgyro_flag:
        print "tgyro"
    if profiles_gen_flag:
        print "profiles_gen"
    print "Try executing pyrats_help with one of these codes as the codename."
else:
    print "It looks like this directory does not contain any NEO, GYRO, TGYRO,"
    print "or profiles_gen files.  Please check the directory and try again."
