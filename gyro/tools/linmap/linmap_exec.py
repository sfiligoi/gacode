#!/usr/bin/env python

# USAGE: 
#
#    python linmap_exec.py <directory> <nProc> <start_over>

import sys
import os
import shutil
import subprocess
import numpy as np

#----------------------------------------
# Capture command-line arguments
#
#directory=$(pwd)

directory = sys.argv[1]

nProc = int(sys.argv[2])
start_over = int(sys.argv[3])
start_over = bool(start_over)


#To solve python rounding issue, i.e. floats are not represented by exact values (important for comparisons).
#MinFloat = sys.float_info.epsilon
MinFloat = pow(10, -sys.float_info.dig)

#print "Min float: %E" % MinFloat
#sys.exit(0)

#print directory

#print nProc

#print start_over

#Go to directory
os.chdir(directory)
#directory=os.getcwd()

print "Starting gyro_linmap in directory : %s" % os.getcwd()
print "**********************************************************"

#cpu=${1}
#start_over=${2}
#-------------------------------------------------------

#file=`cat input.linmap.gen` 

#inputParams = np.loadtxt('input.linmap.gen')
inputParams = np.genfromtxt('input.linmap.gen', dtype='str')
#inputParams = np.genfromtxt('input.linmap.gen')

VariableValue = inputParams[:,0].astype(np.float)
VariableName = inputParams[:,1]

NumVariables = len(VariableValue)

i=0
#for x in $file
#do
#   p=$(($p+1))
#   if [ $p -eq 1 ] ; then DET_TOLERANCE=$x ; fi
#   if [ $p -eq 3 ] ; then ERROR_TOLERANCE=$x ; fi
#   if [ $p -eq 5 ] ; then WI_STABLE_LIMIT=$x ; fi
#if [ $p -eq 7 ] ; then startWR=$x ; fi
#if [ $p -eq 9 ] ; then startWI=$x ; fi
#if [ $p -eq 11 ] ; then RADIUSstart=$x ; fi
#if [ $p -eq 13 ] ; then RADIUSmin=$x ; fi
#if [ $p -eq 15 ] ; then RADIUSmax=$x ; fi
#if [ $p -eq 17 ] ; then RADIUSstep=$x ; fi
#if [ $p -eq 19 ] ; then L_Ystart=$x ; fi
#if [ $p -eq 21 ] ; then L_Ymin=$x ; fi
#if [ $p -eq 23 ] ; then L_Ymax=$x ; fi
#if [ $p -eq 25 ] ; then L_Ystep=$x ; fi
#done

while i < NumVariables:
	if VariableName[i] == 'DET_TOLERANCE':
		DET_TOLERANCE = VariableValue[i]
	elif VariableName[i] == 'ERROR_TOLERANCE':
		ERROR_TOLERANCE = VariableValue[i]
	elif VariableName[i] == 'WI_STABLE_LIMIT':
		WI_STABLE_LIMIT = VariableValue[i]
	elif VariableName[i] == 'startWR':
		startWR = VariableValue[i]
	elif VariableName[i] == 'startWI':
		startWI = VariableValue[i]
	elif VariableName[i] == 'RADIUSstart':
		RADIUSstart = VariableValue[i]
	elif VariableName[i] == 'RADIUSmin':
		RADIUSmin = VariableValue[i]
	elif VariableName[i] == 'RADIUSmax':
		RADIUSmax = VariableValue[i]
	elif VariableName[i] == 'RADIUSstep':
		RADIUSstep = VariableValue[i]
	elif VariableName[i] == 'L_Ystart':
		L_Ystart = VariableValue[i]
	elif VariableName[i] == 'L_Ymin':
		L_Ymin = VariableValue[i]
	elif VariableName[i] == 'L_Ymax':
		L_Ymax = VariableValue[i]
	elif VariableName[i] == 'L_Ystep':
		L_Ystep = VariableValue[i]
	else:
		print "Got false variable" 
	i += 1

print "Parameters:"
print "DET_TOLERANCE: %E" % DET_TOLERANCE
print "ERROR_TOLERANCE: %E" % ERROR_TOLERANCE
print "WI_STABLE_LIMIT: %E" % WI_STABLE_LIMIT
print "startWR: %E" % startWR
print "startWI: %E" % startWI
print "RADIUSstart: %E" % RADIUSstart
print "RADIUSmin: %E" % RADIUSmin
print "RADIUSmax: %E" % RADIUSmax
print "RADIUSstep: %E" % RADIUSstep
print "L_Ystart: %E" % L_Ystart
print "L_Ymin: %E" % L_Ymin
print "L_Ymax: %E" % L_Ymax
print "L_Ystep: %E" % L_Ystep
print "**********************************************************"


#-------------------------------------------------------




#----------------------------------------
# Adjustable: Initial guess
#
#USER INPUT
#startWR=-3.04E-01
#startWI=2.55E-01

#RADIUSstart=0.8
#RADIUSmin=0.20
#RADIUSmax=0.90
#RADIUSstep=0.010

#L_Ystart=0.50
#L_Ymin=0.050
#L_Ymax=2.00
#L_Ystep=0.010

#END OF USER INPUT
#----------------------------------------
#NstepsR=$(echo "scale=0; ($RADIUSmax - $RADIUSmin) / $RADIUSstep;" | bc)
#NstepsL_Y=$(echo "scale=0; ($L_Ymax - $L_Ymin) / $L_Ystep;" | bc)
#echo "Steps in RADIUS = $NstepsR"
#echo "Steps in L_Y = $NstepsL_Y"
NstepsR = int((RADIUSmax - RADIUSmin) / RADIUSstep)
NstepsL_Y= int((L_Ymax - L_Ymin) / L_Ystep)

print "Steps in RADIUS = %i" % NstepsR
print "Steps in L_Y = %i" % NstepsL_Y
#----------------------------------------
#cd $dir
os.chdir(directory)
#if [ ! -f input.gyro.orig ] 
#then
#    echo "No input.gyro.orig found in $dir"
#    exit 1
#fi
if not os.path.exists(directory + '/input.gyro.orig'):
	print "No input.gyro.orig found in %s" % directory
	sys.exit(0)

#cp input.gyro.orig input.gyro
shutil.copyfile('input.gyro.orig', 'input.gyro')

#echo "FIELDEIGEN_WR=$startWR" >> input.gyro
#echo "FIELDEIGEN_WI=$startWI" >> input.gyro



input_gyro = open('input.gyro', 'a')
input_gyro.write("\n")
input_gyro.write("FIELDEIGEN_WR=%E\n" % startWR)
input_gyro.write("FIELDEIGEN_WI=%E\n" % startWI)
input_gyro.close()
#sys.exit(0)
#----------------------------------------

#if $start_over ; then
#echo "Starting from scratch!"
#rm -f radius_krho_fieldeigen_omega.out
#rm -f fieldeigen_gbflux.out
#rm -f rho.out ; touch rho.out

#else
#echo "Restarting!"
#fi

if start_over :
	print "Starting from scratch!"
	if os.path.exists(directory + '/radius_krho_fieldeigen_omega.out'):
		os.remove("radius_krho_fieldeigen_omega.out")
	if os.path.exists(directory + '/fieldeigen_gbflux.out'):
		os.remove("fieldeigen_gbflux.out")
	if os.path.exists(directory + '/rho.out'):
		os.remove("rho.out")
else :
	print "Restarting!"

#radius_krho_fieldeigen_omega = open('radius_krho_fieldeigen_omega.out', 'a')
#fieldeigen_gbflux = open('fieldeigen_gbflux.out', 'a')
#rho = open('rho.out', 'a')



#----------------------------------------
#bc can't use scientific notation so transform everything to regular notation
#DET_TOLERANCE=$(echo $DET_TOLERANCE | sed 's/\([0-9]*\(\.[0-9]*\)\?\)[eE]+\?\(-\?[0-9]*\)/(\1*10^\3)/g;s/^/scale=12;/'| bc)
#ERROR_TOLERANCE=$(echo $ERROR_TOLERANCE | sed 's/\([0-9]*\(\.[0-9]*\)\?\)[eE]+\?\(-\?[0-9]*\)/(\1*10^\3)/g;s/^/scale=12;/'| bc)
#WI_STABLE_LIMIT=$(echo $WI_STABLE_LIMIT | sed 's/\([0-9]*\(\.[0-9]*\)\?\)[eE]+\?\(-\?[0-9]*\)/(\1*10^\3)/g;s/^/scale=12;/'| bc)

#echo "Your tolerance is |det| < $DET_TOLERANCE		error < $ERROR_TOLERANCE"
#echo "An eigenmode is considered stable if gamma < $WI_STABLE_LIMIT"

print "Your tolerance is |det| < %E	error < %E" % (DET_TOLERANCE,  ERROR_TOLERANCE)
print "An eigenmode is considered stable if gamma < %F" % WI_STABLE_LIMIT
#print "**********************************************************"
#----------------------------------------

#restartRadiusWR=$startWR
#restartRadiusWI=$startWI

restartRadiusWR=startWR
restartRadiusWI=startWI

resultsWR=1
resultsWI=1

#currentRADIUS=$RADIUSstart
currentRADIUS=RADIUSstart
#Start by scanning downwards in RADIUS
#downwardsRadius=true
#radiusbreak=false
downwardsRadius=True
radiusbreak=False

#while [ $(echo "$currentRADIUS >= $RADIUSmin && $currentRADIUS <= $RADIUSmax;" | bc) -eq 1 ]; do
#while currentRADIUS >= RADIUSmin and currentRADIUS <= RADIUSmax:
while not (currentRADIUS < RADIUSmin - MinFloat or currentRADIUS > RADIUSmax + MinFloat):

	#radiusbreak=false
	radiusbreak=False
	
	#currentL_Y=$L_Ystart
	currentL_Y=L_Ystart
	#start by scanning downwards in L_Y
	#downwardsL_Y=true
	downwardsL_Y=True
	
	#while [ $(echo "$currentL_Y >= $L_Ymin && $currentL_Y <= $L_Ymax;" | bc) -eq 1 ]; do
	#while currentL_Y >= L_Ymin and currentL_Y <= L_Ymax:
	while not (currentL_Y < (L_Ymin - MinFloat) or currentL_Y > (L_Ymax + MinFloat)):
		
		#converged=false
		converged=False
		trials=0
		#Loop to find a solution
		#while ! $converged  ; do
		while not converged :
			#echo "RADIUS=$currentRADIUS"
			print "**********************************************************"
			print "RADIUS = %F" % currentRADIUS
			#echo "L_Y=$currentL_Y"
			print "L_Y = %F" % currentL_Y
			#echo "RADIUS=$currentRADIUS" >> input.gyro
			#echo "L_Y=$currentL_Y" >> input.gyro
			input_gyro = open('input.gyro', 'a')
			input_gyro.write("\n")
			input_gyro.write("RADIUS=%F\n" % currentRADIUS)
			input_gyro.write("L_Y=%F\n" % currentL_Y)
			input_gyro.close()
			#sys.exit(0)
			#gyro -e . -n $cpu
			#print directory
			subprocess.call(["gyro", "-e", ".", "-n", "%i" %nProc])
			#sys.exit(0)
			#results=`tail -n 1 fieldeigen.out`
			fieldeigen_out = open('fieldeigen.out', 'r')
			lines = fieldeigen_out.readlines()
			results = map(float, lines[-1].split()) #Last line
			fieldeigen_out.close()
			#bar=( $results )
			#resultsWR=${bar[0]}
			#resultsWI=${bar[1]}
			#print "resultsWR : %s" % results
			#test = map(float, results.split())
			#print "resultsWR : %E" % test[1]
			#resultsWR=np.genfromtxt(results[0], dtype='float')
			#resultsWI=results[1].astype(np.float)
			resultsWR=results[0]
			resultsWI=results[1]
			#print "resultsWR : %E" % resultsWR
			#print "resultsWI : %E" % resultsWI
			#det=$(echo "${bar[2]}")
			#error=${bar[3]}
			det=results[2]
			error=results[3]
			#print "det : %E" % det
			#print "error : %E" % error
			#print "size : %i" % len(results)
			#sys.exit(0)

			#bc can't use scientific notation so transform everything to regular notation
			#resultsWR=$(echo $resultsWR | sed 's/\([0-9]*\(\.[0-9]*\)\?\)[eE]+\?\(-\?[0-9]*\)/(\1*10^\3)/g;s/^/scale=12;/'| bc)
			#resultsWI=$(echo $resultsWI | sed 's/\([0-9]*\(\.[0-9]*\)\?\)[eE]+\?\(-\?[0-9]*\)/(\1*10^\3)/g;s/^/scale=12;/'| bc)
			#det=$(echo $det | sed 's/\([0-9]*\(\.[0-9]*\)\?\)[eE]+\?\(-\?[0-9]*\)/(\1*10^\3)/g;s/^/scale=15;/'| bc)
			#error=$(echo $error | sed 's/\([0-9]*\(\.[0-9]*\)\?\)[eE]+\?\(-\?[0-9]*\)/(\1*10^\3)/g;s/^/scale=15;/'| bc)


			#if [ $(echo "$det <= $DET_TOLERANCE && $error <= $ERROR_TOLERANCE;" | bc) -eq 1 ]; then
			if det <= (DET_TOLERANCE + MinFloat) and error <= (ERROR_TOLERANCE + MinFloat):
				print "Converged!"
				#converged=true
				converged=True
				break
			#fi
			#sys.exit(0)
			#trials=$(echo "$trials + 1;" | bc)
			trials += 1
			#if [ $(echo "$trials >= 3;" | bc) -eq 1 ]; then #Give up searching
			if trials >= 3 : #Give up searching
				#if [ $(echo "$currentL_Y == $L_Ystart;" | bc) -eq 1 ]; then
				if currentL_Y == L_Ystart :
					#radiusbreak=true
					radiusbreak=True
				#fi
				break
			#fi
			
			#if [ $(echo "$currentL_Y == $L_Ystart;" | bc) -eq 1 ]; then
			if currentL_Y == L_Ystart :

				#if [ $(echo "$trials == 1;" | bc) -eq 1 ]; then #Try reducing step in RADIUS to half
				if trials == 1 : #Try reducing step in RADIUS to half
					#echo "Try reducing step length in RADIUS to half"
					print "Try reducing step length in RADIUS to half"
					#if $downwardsRadius ; then
					if downwardsRadius :
						#currentRADIUS=$(echo "$currentRADIUS + 0.5000*$RADIUSstep;" | bc)
						currentRADIUS = currentRADIUS + 0.5000*RADIUSstep
					else :
						#currentRADIUS=$(echo "$currentRADIUS - 0.5000*$RADIUSstep;" | bc)
						currentRADIUS = currentRADIUS - 0.5000*RADIUSstep
					#fi
				#fi
				#if [ $(echo "$trials == 2;" | bc) -eq 1 ]; then #Try reducing step in RADIUS to 1/4th
				if trials == 2 : #Try reducing step in RADIUS to 1/4th
					#echo "Try reducing step length in RADIUS to 1/4th"
					print "Try reducing step length in RADIUS to 1/4th"
					#if $downwardsRadius ; then
					if downwardsRadius :
						#currentRADIUS=$(echo "$currentRADIUS + 0.2500*$RADIUSstep;" | bc)
						currentRADIUS = currentRADIUS + 0.2500*RADIUSstep
					else :
						#currentRADIUS=$(echo "$currentRADIUS - 0.2500*$RADIUSstep;" | bc)
						currentRADIUS = currentRADIUS - 0.2500*RADIUSstep
					#fi
				#fi

			else :

				#if [ $(echo "$trials == 1;" | bc) -eq 1 ]; then #Try reducing step in L_Y to half
				if trials == 1 : #Try reducing step in L_Y to half
					#echo "Try reducing step length in L_Y to half"
					print "Try reducing step length in L_Y to half"
					#if $downwardsL_Y ; then
					if downwardsL_Y :
						#currentL_Y=$(echo "$currentL_Y + 0.5000*$L_Ystep;" | bc)
						currentL_Y = currentL_Y + 0.5000*L_Ystep
					else :
						#currentL_Y=$(echo "$currentL_Y - 0.5000*$L_Ystep;" | bc)
						currentL_Y = currentL_Y - 0.5000*L_Ystep
					#fi
				#fi
				#if [ $(echo "$trials == 2;" | bc) -eq 1 ]; then #Try reducing step in L_Y to 1/4th
				if trials == 2 : #Try reducing step in L_Y to 1/4th
					#echo "Try reducing step length in L_Y to 1/4th"
					print "Try reducing step length in L_Y to 1/4th"
					#if $downwardsL_Y ; then
					if downwardsL_Y :
						#currentL_Y=$(echo "$currentL_Y + 0.2500*$L_Ystep;" | bc)
						currentL_Y = currentL_Y + 0.2500*L_Ystep
					else :
						#currentL_Y=$(echo "$currentL_Y - 0.2500*$L_Ystep;" | bc)
						currentL_Y = currentL_Y - 0.2500*L_Ystep
					#fi
				#fi

			#fi
			
		#done
		#input_gyro.close()

		#if  $converged  ; then #Write eigenvalue to output files if converged
		if converged : #Write eigenvalue to output files if converged
			#if [ $(echo "$resultsWI >= $WI_STABLE_LIMIT;" | bc) -eq 1 ]; then
			if resultsWI >= (WI_STABLE_LIMIT - MinFloat) :

				#echo "Writing to output"
				print "Writing to output"
				#echo "$currentRADIUS	$currentL_Y	$results" >> radius_krho_fieldeigen_omega.out
				#radius_krho_fieldeigen_omega.write("%F\t%F\t%S\n" % currentRADIUS % currentL_Y % results)
				radius_krho_fieldeigen_omega = open('radius_krho_fieldeigen_omega.out', 'a')
				radius_krho_fieldeigen_omega.write("%F\t%F\t%E\t%E\t%E\t%E\n" % (currentRADIUS, currentL_Y, results[0], results[1], results[2], results[3]))
				radius_krho_fieldeigen_omega.close()
				#sys.exit(0)

				#tail -1 out.gyro.gbflux >> fieldeigen_gbflux.out
				out_gyro_gbflux = open('out.gyro.gbflux', 'r')
				lines = out_gyro_gbflux.readlines()
				fieldeigen_gbflux = open('fieldeigen_gbflux.out', 'a')
				fieldeigen_gbflux.write("%s\n" % lines[-1]) #Last line
				fieldeigen_gbflux.close()
				out_gyro_gbflux.close()
				
				#####CHECK
				#profiles_gen -i input.profiles -loc_rad "$currentRADIUS" | grep "rho " >> rho.out
				#rho
			#fi
		else :
			#echo "Couldn't find solution! Stop following branch"
			print "Couldn't find solution! Stop following branch"

			#if $radiusbreak ; then
			if radiusbreak :
				break
			#fi
			#if $downwardsL_Y ; then #Start scanning upwards
			if downwardsL_Y : #Start scanning upwards
				#downwardsL_Y=false
				downwardsL_Y=False
				#currentL_Y=$(echo "$L_Ystart + $L_Ystep;" | bc)
				currentL_Y = L_Ystart + L_Ystep
				#cp input.gyro.orig input.gyro
				shutil.copyfile('input.gyro.orig', 'input.gyro')
				#echo "FIELDEIGEN_WR=$restartRadiusWR" >> input.gyro
				#echo "FIELDEIGEN_WI=$restartRadiusWI" >> input.gyro
				input_gyro = open('input.gyro', 'a')
				input_gyro.write("\n")
				input_gyro.write("FIELDEIGEN_WR=%E\n" % restartRadiusWR)
				input_gyro.write("FIELDEIGEN_WI=%E\n" % restartRadiusWI)
				input_gyro.close()
				continue
			else :
				break #Can't continue on branch
			#fi
		#fi


		#if [ $(echo "$resultsWI < $WI_STABLE_LIMIT;" | bc) -eq 1 ]; then #Mode has become stable
		if resultsWI < (WI_STABLE_LIMIT - MinFloat) : #Mode has become stable
			#echo "Branch is stable! Stop following"
			print "Branch is stable! Stop following"
			#if $downwardsL_Y ; then #Start scanning upwards
			if downwardsL_Y : #Start scanning upwards
				#downwardsL_Y=false
				downwardsL_Y=False
				#currentL_Y=$(echo "$L_Ystart + $L_Ystep;" | bc)
				currentL_Y = L_Ystart + L_Ystep
				#cp input.gyro.orig input.gyro
				shutil.copyfile('input.gyro.orig', 'input.gyro')
				#echo "FIELDEIGEN_WR=$restartRadiusWR" >> input.gyro
				#echo "FIELDEIGEN_WI=$restartRadiusWI" >> input.gyro
				input_gyro = open('input.gyro', 'a')
				input_gyro.write("\n")
				input_gyro.write("FIELDEIGEN_WR=%E\n" % restartRadiusWR)
				input_gyro.write("FIELDEIGEN_WI=%E\n" % restartRadiusWI)
				input_gyro.close()
				continue
			else :
				break #Can't continue on branch
			#fi
		#fi


		#if [ $(echo "$currentL_Y == $L_Ystart;" | bc) -eq 1 ]; then #save solution for later
		if currentL_Y == L_Ystart : #save solution for later
			#restartRadiusWR=$resultsWR
			#restartRadiusWI=$resultsWI
			restartRadiusWR=resultsWR
			restartRadiusWI=resultsWI
		#fi
		
		#cp input.gyro.orig input.gyro
		shutil.copyfile('input.gyro.orig', 'input.gyro')
		#echo "FIELDEIGEN_WR=$resultsWR" >> input.gyro
		#echo "FIELDEIGEN_WI=$resultsWI" >> input.gyro
		input_gyro = open('input.gyro', 'a')
		input_gyro.write("\n")
		input_gyro.write("FIELDEIGEN_WR=%E\n" % resultsWR)
		input_gyro.write("FIELDEIGEN_WI=%E\n" % resultsWI)
		input_gyro.close()
	
		#if [ $(echo "$currentL_Y - $L_Ystep < $L_Ymin;" | bc) -eq 1 ]; then #Change to scanning upwards
		if (currentL_Y - L_Ystep) < (L_Ymin - MinFloat) : #Change to scanning upwards
			#downwardsL_Y=false
			downwardsL_Y = False
			#currentL_Y=$L_Ystart
			currentL_Y = L_Ystart
			#echo "FIELDEIGEN_WR=$restartRadiusWR" >> input.gyro
			#echo "FIELDEIGEN_WI=$restartRadiusWI" >> input.gyro
			input_gyro = open('input.gyro', 'a')
			input_gyro.write("FIELDEIGEN_WR=%E\n" % restartRadiusWR)
			input_gyro.write("FIELDEIGEN_WI=%E\n" % restartRadiusWI)
			input_gyro.close()
		#fi

		#if $downwardsL_Y ; then
		if downwardsL_Y :
			#currentL_Y=$(echo "$currentL_Y - $L_Ystep;" | bc)
			currentL_Y = currentL_Y - L_Ystep
		else :
			#currentL_Y=$(echo "$currentL_Y + $L_Ystep;" | bc)
			currentL_Y = currentL_Y + L_Ystep
		#fi
	#done
	#input_gyro.close()
		#print "currentL_Y: %E" % currentL_Y
		#print "L_Ymax: %E" % L_Ymax
		#print "currentL_Y <= L_Ymax: %i" % (not(currentL_Y > L_Ymax))

	#if $radiusbreak ; then
	if radiusbreak :
		#if $downwardsRadius ; then
		if downwardsRadius :
			#downwardsRadius=false
			downwardsRadius = False
			#currentRADIUS=$(echo "$RADIUSstart + $RADIUSstep;" | bc)
			currentRADIUS = RADIUSstart + RADIUSstep
			#cp input.gyro.orig input.gyro
			shutil.copyfile('input.gyro.orig', 'input.gyro')
			#echo "FIELDEIGEN_WR=$startWR" >> input.gyro
			#echo "FIELDEIGEN_WI=$startWI" >> input.gyro
			input_gyro = open('input.gyro', 'a')
			input_gyro.write("\n")
			input_gyro.write("FIELDEIGEN_WR=%E\n" % startWR)
			input_gyro.write("FIELDEIGEN_WI=%E\n" % startWI)
			input_gyro.close()
			continue
		else :
			break #Can't continue on branch, stop execution
		#fi
	#fi

	#cp input.gyro.orig input.gyro
	shutil.copyfile('input.gyro.orig', 'input.gyro')

	#echo "FIELDEIGEN_WR=$restartRadiusWR" >> input.gyro
	#echo "FIELDEIGEN_WI=$restartRadiusWI" >> input.gyro
	input_gyro = open('input.gyro', 'a')
	input_gyro.write("\n")
	input_gyro.write("FIELDEIGEN_WR=%E\n" % restartRadiusWR)
	input_gyro.write("FIELDEIGEN_WI=%E\n" % restartRadiusWI)
	input_gyro.close()

	#if [ $(echo "$currentRADIUS - $RADIUSstep < $RADIUSmin;" | bc) -eq 1 ]; then #Change to scanning upwards
	if (currentRADIUS - RADIUSstep) < (RADIUSmin - MinFloat) : #Change to scanning upwards
		#downwardsRadius=false
		downwardsRadius = False
		#currentRADIUS=$RADIUSstart
		currentRADIUS = RADIUSstart
		#echo "FIELDEIGEN_WR=$startWR" >> input.gyro
		#echo "FIELDEIGEN_WI=$startWI" >> input.gyro
		input_gyro = open('input.gyro', 'a')
		input_gyro.write("FIELDEIGEN_WR=%E\n" % startWR)
		input_gyro.write("FIELDEIGEN_WI=%E\n" % startWI)
		input_gyro.close()
	#fi

	#if $downwardsRadius ; then
	if downwardsRadius :
		#currentRADIUS=$(echo "$currentRADIUS - $RADIUSstep;" | bc)
		currentRADIUS = currentRADIUS - RADIUSstep 
	else :
		#currentRADIUS=$(echo "$currentRADIUS + $RADIUSstep;" | bc)
		currentRADIUS = currentRADIUS + RADIUSstep
	#fi
#done

#Close output files
#radius_krho_fieldeigen_omega.close()
#fieldeigen_gbflux.close()
#rho.close()

#exit
sys.exit(0)

