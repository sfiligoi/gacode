import tables
a=tables.openFile("gyro05800.h5","a")
b=tables.openFile("gyroMesh.h5","r")
a.removeNode(a.root.coarseMesh)
a.copyNode(b.root.coarseMesh,a.root,recursive=True)
a.close()
b.close()
