
all:
	cd shared ; make
	cd gyro ; make 
	cd neo ; make 
	cd gato ; make 
	cd tglf ; make
	cd tgyro ; make

clean:
	cd shared ; make clean
	cd gyro ; make clean 
	cd neo ; make clean 
	cd gato ; make clean 
	cd tglf ; make clean
	cd tgyro ; make clean
	cd gks ; make clean
	cd xptor ; make clean 
	cd gkcoll ; make clean
	cd le3 ; make clean
	rm -f shared/python/pyrats/*/*.pyc
	rm -f modules/*genmod*
