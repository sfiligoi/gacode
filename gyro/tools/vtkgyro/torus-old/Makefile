CC = g++ -O3 -DNDEBUG
LIB=/usr/local/lib/
INCLUDE=/usr/local/include/vtk-5.0

EXEC1 = inclgen
OBJECTS1 = $(EXEC1).o

EXEC2 = vtkgyro
OBJECTS2 = $(EXEC2).o

.SUFFIXES : .o .C

.C.o :
	$(CC) -c -I${INCLUDE} $<

all : vtkgyro inclgen

inclgen : $(OBJECTS1) $(EXEC1).o
	$(CC) -o $(EXEC1) $(OBJECTS1)

vtkgyro : $(OBJECTS2) $(EXEC2).o
	$(CC) -o $(EXEC2) $(OBJECTS2) -L$(LIB) \
        -lvtkCommon -lvtkFiltering -lvtkGenericFiltering \
        -lvtkGraphics -lvtkHybrid -lvtkImaging -lvtkIO \
        -lvtkjpeg -lvtkRendering -lvtkVolumeRendering \
        -lvtkftgl -lvtkGenericFiltering -lvtksys

clean :
	rm -f $(EXEC1) $(EXEC2) *.o *~ 

jpg:
	touch a.jpg incl
	rm *.jpg incl
