#!/usr/bin/env python
"""
  Class with members and methods associated with parsing the "ITER DB
  Files"
"""

import re, sys, os, glob, filecmp, datetime, shutil
import numpy
from string_to_number import *
from iterDBsemanticMap import *

class iterDB:
     def __init__(self,iterDBfile):
          """ Parse and put into a class the data in an iterDB file.
          Note that that it's called iterDB, but there really isn't a DB"""
          #-------------------------------------------------
          # Set various initial members and params
          #-------------------------------------------------
          self.filename=iterDBfile
          self.varNames=[]
          self.variables={}
          idb=open(iterDBfile,"r")
          line=idb.readline()
          # First line
          if line[0:2]=="**":
            self.header=line
            line=idb.readline()
            if not line: return

          # Parse the file line by line.  It's a fairly simple format
          allValues=""
          while 1:
            # Lines with *'s mean's that we have new variable
            # If we have values need to stick them into the 
            # current var before carrying on to the new var
            if line[0]=="*": 
              if len(allValues)>0 and parseStatus != "Not Found": 
               # Some values are strings; e.g., species
               try:
                  a=numpy.array(allValues.split()).astype(numpy.double)
               except:
                  a=allValues.split()
               self.variables[varname]["values"]=a
               allValues=""
              
              # Start new variable.  Begin by determining line
              line=line.strip("*").strip()
              varname,parseStatus=self.determineVarname(idb,line)

            # Occasionally we parse values right after parsing name
            if parseStatus=="Done": 
              line=idb.readline()
              continue

            # Now parse the values.  If it's not variable, it's values
            line=idb.readline()
            if not line: break
            if not line[0]=="*": 
              allValues=allValues+" "+line
          return

     def determineVarname(self,idb,line):
       """ The possibility of having a single word represent a varname
       is often but not always possible.  Generally, if the line has a
       colon, then the word varname is present.   Other cases are harder
       """
       testline=line.split("species:")[0]
       if testline.find(":")>=0:
         # These are the easiest to parse
         varname,parseStatus=self.determineColonVarnames(line)
       else:
         varname,parseStatus=self.determineOtherVarnames(idb,line)
       return varname,parseStatus


     def determineColonVarnames(self,line):
       """ Fairly straightforward parsing here 
       Example: sbion : beam thermal ion source, #/(meter**3*second)
       Species do add a wrinkle.  Example:
         sion : source due to ionization, #/(meter**3*second), species: d       
       """
       varname=line.split(":")[0].strip()
       varname, line=self.fixSpecies(varname,line)
       self.varNames.append(varname)
       self.variables[varname]={}
       desc, units=self.getDescriptionUnits(line[1:].split(":")[1].strip())
       self.variables[varname]["description"]=desc
       self.variables[varname]["units"]=units
       return varname,"Not Done"

     def determineOtherVarnames(self,idb,line):
       """ The tricky stuff """
       # Try to autodetermine name for case of finding a single word
       # EXAMPLE: qrfe, RF electron heating, watts/meter**3
       if len(line.split(",")[0].split()) == 1:
         varname=line.split(",")[0].split()[0].strip()
         if varname=="total" or varname=="wdot":
           varname, status=self.areYouKiddingMe(varname,idb,line)
           return varname, status
         varname, line=self.fixSpecies(varname,line)
         self.varNames.append(varname)
         self.variables[varname]={}
         self.variables[varname]["description"]=line
         desc, units=self.getDescriptionUnits(line.lstrip(line.split(",")[0]+","))
         self.variables[varname]["description"]=desc
         self.variables[varname]["units"]=units
         return varname, "Not Done"

       # Now we are hosed and and need to do brute force pattern
       # matching using the iterDBmap to find the names
       else:
         for key in iterDBmap.keys():
           if line.find(key)>=0:
             varname=iterDBmap[key]
             varname, line=self.fixSpecies(varname,line)
             self.varNames.append(varname)
             self.variables[varname]={}
             desc, units=self.getDescriptionUnits(line)
             self.variables[varname]["description"]=desc
             self.variables[varname]["units"]=units
             return varname, "Not Done"
         print  "WARNING: Cannot determine variable name for: "
         print  line
         return "NA", "Not Found"


     def areYouKiddingMe(self,varname,idb,line):
       """ Yes, there are strings that are impossible to parse in any sane
       way because fusion scientists love to hard code 
       Currently I am handling the lines:
       *  wdot, electrons, watts/meter**3
       *  wdot, ions, watts/meter**3
       *  total, ohmic, bootstrap, beam and RF currents, amps
       For total, ohmic, bootstrap, beam and RF currents, amps
       special parsing is needed as well
       """
       if varname=="wdot":
         lspl=line.split(",")
         varname=lspl[0].strip()+lspl[1].strip()
         self.varNames.append(varname)
         self.variables[varname]={}
         self.variables[varname]["description"]=varname
         self.variables[varname]["units"]=lspl[2]
         return varname, "Not Done"
       elif varname=="total":
         if line.find("ohmic")>=0: 
           self.parseCurrents(idb,line)
           return "NA", "Done"
         else:
           print  "WARNING: Cannot determine variable name for: "
           print  line
           return "NA", "Not Found"


     def getDescriptionUnits(self,desc_units):
       """ Given the tail end of a string, get the description and
       units:  examples:
         qrfe, RF electron heating, watts/meter**3
         qrfe: RF electron heating, watts/meter**3
       Pass in: RF electron heating, watts/meter**3 and it'll split it """
       description=desc_units.split(",")[0].strip()
       if desc_units.find(",")>0:
         units=desc_units.split(",")[1].strip()
       else:
         units=" "
       return description, units

     def parseCurrents(self,idb,line):
       """ I mean really, we couldn't put them on separate lines? """
       line=idb.readline()
       #NEED TO DO SOME MORE PARSING
       return

     def fixSpecies(self,varname,line):
       """ If species: is on the line, then fix the varname and 
           strip the species part from the line.  E.g.,
             sion : source due to ionization, #/(meter**3*second), species: d       
            will end up with varname=sion_d and  line equal to:
             sion : source due to ionization, #/(meter**3*second)
       """
       if line.find("species:")>=0:
          species=line.split("species:")[1].strip()
          varname=varname+"_"+species
          line=line.split("species:")[0].strip().rstrip(",")
       return varname, line

     def printDbInfo(self):
       """ Print information on the class in a nice format """
       print "Filename: ", self.filename
       print "Header:   ", self.header
       for var in self.variables.keys():
           print "  Variable: ", var
           if  self.variables[var].has_key('description'):
             if  len(self.variables[var]['description'].strip())>0:
               print "     Description: ", self.variables[var]['description']
           if  self.variables[var].has_key('units'):
             if  len(self.variables[var]['units'].strip())>0:
               print "     Units:       ", self.variables[var]['units']

              
 
