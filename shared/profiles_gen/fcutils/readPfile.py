#!/usr/bin/env python
"""
  Class with members and methods associated with parsing the "ITER DB
  Files"
"""

import re, sys, os, glob, filecmp, datetime, shutil
import numpy
from string_to_number import *

class variable:
     def __init__(self,varSize):
        self.values=numpy.zeros(varSize)
        self.derivs=numpy.zeros(varSize)
        self.grid=numpy.zeros(varSize)
        self.units=""

class pfileData(object):
     def __init__(self,pfile):
       """ Parse and put into a class the data in an iterDB file.
       Note that that it's called iterDB, but there really isn't a DB"""
       #-------------------------------------------------
       # Set various initial members and params
       #-------------------------------------------------
       self.filename=pfile
       self.listVars=[]
       self.listGrid=[]
       self.listDeriv=[]
       self.listLists=[]
       self.listUnits=[]
       idb=open(pfile,"r")


       # Parse the file line by line.  It's a fairly simple format
       while 1:
         line=idb.readline()
         if not line: break
         # First line
         headlist=line.split()
         arrSize=numpy.int(headlist[0])
         gridName=headlist[1]
         varName=headlist[2]
         varDerivName=headlist[3]
         varUnits=""
         if len(varName.strip(')').split('('))>1:
            varUnits=varName.strip(')').split('(')[1]
            varName=varName.strip(')').split('(')[0]
         self.listVars.append(varName)
         self.listUnits.append(varUnits)

         # This is a clever way of naming the class member
         # Somewhat ugly but it gives tab completion which is sweet
         self.__dict__[varName]=variable(arrSize)
         self.__dict__[varName].units=varUnits

         for iline in range(arrSize):
           gridVal,val,deriv=idb.readline().split()
           self.__dict__[varName].values[iline]=val
           self.__dict__[varName].grid[iline]=gridVal

       return


