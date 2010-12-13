/*==========================================================
  module: vtkgyro.C
 
  Top-level program controlling cap AND torus
  rendering.

  ==========================================================*/

#include "vtkMath.h"
#include "vtkRenderer.h"
#include "vtkRenderWindow.h"
#include "vtkStructuredGrid.h"
#include "vtkStructuredGridGeometryFilter.h"
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
#include "vtkRenderLargeImage.h"
#include "vtkPiecewiseFunction.h"

#include <iostream>
#include <fstream>
#include <string>
#include <complex>
using namespace std;

#include "globals.h"
#include "../common/parser.h"
#include "../common/make_name.h"
#include "../common/read_all.h"
#include "func_sq.h"
#include "func_sq2.h"
#include "do_square.h"
#include "read_VTK_parameters.h"

int main( int argc, char *argv[] )
{
  path="/home/candy/SIM/etg-ki/";

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

  cout << "===========" << endl;
  cout << "SQUARE MODE" << endl;
  cout << "===========" << endl;

  do_square();

  exit(0);

}

