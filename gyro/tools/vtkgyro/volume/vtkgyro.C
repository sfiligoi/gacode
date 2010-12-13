/*==========================================================
  module: vtkgyro.C
 
  Top-level program controlling cap AND torus
  rendering.

  ==========================================================*/

#include "vtkMath.h"
#include "vtkRenderer.h"
#include "vtkRenderWindow.h"
#include "vtkStructuredPoints.h"
#include "vtkStructuredGrid.h"
#include "vtkStructuredGridGeometryFilter.h"
#include "vtkPoints.h"
#include "vtkDataSet.h"
#include "vtkPointData.h"
#include "vtkActor.h"
#include "vtkLookupTable.h"
#include "vtkPolyDataMapper.h"
#include "vtkFloatArray.h"
#include "vtkUnsignedShortArray.h"
#include "vtkRendererSource.h"
#include "vtkJPEGWriter.h"
#include "vtkCamera.h"
#include "vtkProperty.h"
#include "vtkLight.h"
#include "vtkRenderLargeImage.h"
#include "vtkPiecewiseFunction.h"
#include "vtkColorTransferFunction.h"
#include "vtkVolumeProperty.h"
#include "vtkVolumeMapper.h"
#include "vtkVolume.h"
#include "vtkVolumeRayCastCompositeFunction.h"
#include "vtkVolumeRayCastMapper.h"

#include <iostream>
#include <fstream>
#include <string>
#include <complex>
using namespace std;

#include "globals.h"
#include "parser.h"
#include "func.h"
#include "func_sq.h"
#include "c_func.h"
#include "make_name.h"
#include "read_all.h"
#include "do_square.h"
#include "do_fancy.h"
#include "do_fancy_2.h"
#include "read_VTK_parameters.h"

int main( int argc, char *argv[] )
{

  string simdir;
  string path="/home/candy/gyro/sim/";

  // Parse command line and set sim directory:

  if (argc != 2) {
    cout << "Usage: vtkgyro <simdir>" << endl;
    exit(1);
  }

  simdir = argv[1];
  cout << "Parsing information in " << simdir << endl;

  i_field = 0;

  read_VTK_parameters();
  read_all(path+simdir);

  switch (vis_method) {

  case 1:

    cout << "===========" << endl;
    cout << "SQUARE MODE" << endl;
    cout << "===========" << endl;

    n_alpha = n_phi;

    do_square();
    break;

  case 2:

    cout << "==========" << endl;
    cout << "FANCY MODE" << endl;
    cout << "==========" << endl;

    do_fancy_2();
    break;

  default :
 
    cout << "Bad vis_method." << endl; 
    exit(1);
    break;

  }

}

