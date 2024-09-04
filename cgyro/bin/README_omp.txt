Running CGYRO with high nomp
============================

CGYRO allocates some of its data on the stack, and the amount is 
proportional to the number of omp (OpenMP) threads you use. The 
default stack limit of 8M is only sufficient for up to about 4 omp 
threads. Thus, when using CGYRO with openMP, it is recommended 
that you put the following line into your .bashrc 

ulimit -s 32768 

(This must be set on all the nodes running CGYRO). Moreover, for 
platforms that do not use the system stack limit for OpenMP 
threads, also set

export OMP_STACKSIZE=32M

This should be sufficient for at least 128 omp threads.
