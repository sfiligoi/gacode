/*==========================================================
  module: vtkgyro.C
 
  Top-level program controlling cap AND torus
  rendering.

  Geometry: (r,theta,phi)

  Caps:

  dim[0] = n_r;
  dim[1] = n_theta;
  dim[2] = 1;

  Ends:

  dim[0] = 1;
  dim[1] = n_theta;
  dim[2] = n_phi;

  ==========================================================*/

#include "vtkMath.h"
#include "vtkRenderer.h"
#include "vtkRenderWindow.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkStructuredGrid.h"
#include "vtkStructuredGridGeometryFilter.h"
#include "vtkFloatArray.h"
#include "vtkPoints.h"
#include "vtkPointData.h"
#include "vtkActor.h"
#include "vtkLookupTable.h"
#include "vtkPolyDataMapper.h"
#include "vtkFloatArray.h"
#include "vtkRendererSource.h"
#include "vtkJPEGWriter.h"
#include "vtkCamera.h"
#include "vtkProperty.h"
#include "vtkLight.h"
#include "vtkRenderLargeImage.h"

#include <iostream>
#include <fstream>
#include <string>
#include <complex>
using namespace std;

#include "globals.h"
#include "../common/parser.h"
#include "../common/make_name.h"
#include "../common/read_all.h"
#include "func.h"
#include "c_func.h"
#include "do_torus.h"
#include "do_supertorus.h"
#include "do_cap.h"
#include "read_VTK_parameters.h"

int main( int argc, char *argv[] )
{

  path="/home/candy/";

  // Parse command line and set sim directory:

  if (argc != 2) {
    cout << "Usage: vtkgyro <simdir>" << endl;
    exit(1);
  }

  simdir = argv[1];
  cout << "Parsing information in " << simdir << endl;

  i_field = 0;

  read_VTK_parameters();
  read_all();

  switch (vis_method) {

  case 1:

    cout << "========" << endl;
    cout << "CAP MODE" << endl;
    cout << "========" << endl;

    do_cap();
    break;
   
  case 2:

    cout << "==========" << endl;
    cout << "TORUS MODE" << endl;
    cout << "==========" << endl;

    do_torus();
    //do_supertorus();
    break;

  default :
 
    cout << "Bad vis_method." << endl; 
    exit(1);
    break;

  }

}

