all:
	cd shared ; make
	cd le3 ; make
	cd gato ; make 
	cd neo ; make 
	cd glf23 ; make
	cd tglf ; make
	cd gyro ; make 
	cd tgyro ; make

clean:
	cd shared ; make clean
	cd le3 ; make clean
	cd gato ; make clean 
	cd neo ; make clean 
	cd glf23 ; make clean
	cd tglf ; make clean
	cd gyro ; make clean 
	cd cgyro ; make clean 
	cd tgyro ; make clean
	rm -f shared/python/pyrats/*/*.pyc
	rm -f shared/python/*/*.pyc
	rm -f modules/*genmod*
