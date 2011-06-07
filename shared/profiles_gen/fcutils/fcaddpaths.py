#!/usr/bin/env python
#
# $Id: fcaddpaths.py 2601 2009-07-22 21:10:19Z cary $

import sys, os

def fcaddpaths(fcutilsexecdir):
  facetsdir = os.path.dirname(fcutilsexecdir)
  # print "fcutilsexecdir =", fcutilsexecdir
  # print "facetsdir =", facetsdir
  for i in ("bin", "share", "fcutils", "txutils"):
    canddir = os.path.join(facetsdir, i)
    # print "canddir =", canddir
    if canddir not in sys.path and os.path.exists(canddir):
      sys.path.append(canddir)
    canddir2 = os.path.join(canddir, "fileconvertors")
    if canddir2 not in sys.path and os.path.exists(canddir2):
      sys.path.append(canddir2)

