#!/usr/bin/env python
"""
  Class with members and methods associated with parsing the "ITER DB
  Files"
"""
import re, sys, os, glob, filecmp, datetime, shutil
import numpy
from string_to_number import *

class ipVariable:
     def __init__(self,varSize):
        self.values=numpy.zeros(varSize)
        self.units=""

class inputProfileData(object):
     def __init__(self,ifile):
       """ Parse and put into a class the data in an iterDB file.
       Note that that it's called iterDB, but there really isn't a DB"""
       #-------------------------------------------------
       # Set various initial members and params
       #-------------------------------------------------
       self.filename=ifile
       self.listVars=[]
       self.listUnits=[]
       self.header=[]
       idb=open(ifile,"r")

       # Read header
       while 1:
         line=idb.readline()
         if line.startswith("#"):
            self.header.append(line)
            continue
         else:
            break
       # The order is hard-coded here.  Should generalize
       self.n_exp=numpy.int(line.split("=")[1])
       line=idb.readline()
       self.bt_exp=numpy.float(line.split("=")[1])
       line=idb.readline()
       self.arho_exp=numpy.float(line.split("=")[1])
       print self.n_exp

       # Parse the file line by line.  It's a fairly simple format
       while 1:
         line=idb.readline().lstrip("#")
         if not line: break
         # Header that contains the variable lists.  Easier to remove
         # spaces between varname and first parenthesis for parsing
         headlist=line.replace(" (","(").split()
         if len(headlist)==0: break
         print "=> ", headlist
         curList=[]
         for ivar in range(5):
           varName=headlist[ivar]
           if len(varName.strip(')').split('('))>1:
              varUnits=varName.strip(')').split('(')[1]
              varName=varName.strip(')').split('(')[0]
           self.listVars.append(varName)
           self.listUnits.append(varUnits)
           curList.append(varName)
           # This is a clever way of naming the class member
           # Somewhat ugly but it gives tab completion which is sweet
           self.__dict__[varName]=ipVariable(self.n_exp)
           self.__dict__[varName].units=varUnits

         for iline in range(self.n_exp):
           values=idb.readline().split()
           for ivar in range(5):
             varName=curList[ivar]
             self.__dict__[varName].values[iline]=values[ivar]

       return

     def writeInputFile(self,ofile):
       idb=open(ofile,"w")
       nrows=len(self.listVars)/5
       # write header
       for iline in range(len(self.header)):
          idb.write(self.header[iline])
       # write scalar data
       idb.write("N_EXP = "+str(self.n_exp)+"\n")
       numsci="%5.6e" % self.bt_exp
       idb.write("BT_EXP = "+numsci +"\n")
       numsci="%5.6e" % self.arho_exp
       idb.write("ARHO_EXP = "+numsci +"\n")
       # Write out the profiles
       ivar=0
       nrows=len(self.listVars)/5
       for irow in range(nrows):
          # Write column headers
          idb.write("# ")
          ivarstart=ivar
          for icol in range(5):
            label=self.listVars[ivar]+"("+self.listUnits[ivar]+")"
            idb.write(label.center(12)+"   ")
            ivar=ivar+1
          idb.write("\n ")
          for ivalrow in range(self.n_exp):
            ivar=ivarstart
            for icol in range(5):
              varName=self.listVars[ivar]
              value="%5.6e" % self.__dict__[varName].values[ivalrow]
              if icol<4:
                 idb.write(value.upper()+"   ")
              else:
                 idb.write(value.upper())
              ivar=ivar+1
            idb.write("\n")
            if not ivalrow == self.n_exp-1: idb.write(" ")

