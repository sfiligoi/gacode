CC = g++ -O3 -DNDEBUG
LIB=/usr/local/lib/
INCLUDE=/usr/local/include/vtk-5.0

EXEC2 = vtkgyro
OBJECTS2 = $(EXEC2).o

.SUFFIXES : .o .C

.C.o :
	$(CC) -c -I${INCLUDE} $<

all : vtkgyro

vtkgyro : $(OBJECTS2) $(EXEC2).o
	$(CC) -o $(EXEC2) $(OBJECTS2) -L$(LIB) \
        -lvtkCommon -lvtkFiltering -lvtkGenericFiltering \
        -lvtkGraphics -lvtkHybrid -lvtkImaging -lvtkIO \
        -lvtkRendering -lvtkVolumeRendering \
        -lvtkftgl -lvtkGenericFiltering -lvtksys

clean :
	rm -f $(EXEC2) *.o *~ 

jpg:
	touch a.jpg
	rm *.jpg
