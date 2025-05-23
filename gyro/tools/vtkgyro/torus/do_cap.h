/*==========================================================
  module: do_cap.h
 
  Cap rendering..

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

void do_cap( void ) {
  
  char filename[13];
  int i;
  int j;
  int k;
  int offset;

  float s;
  float x[3];


  d_theta = 2*pi/(n_theta-1);

  //-------------------------------------------------------------
  // Time-invariant objects:
  //
  vtkLookupTable *table = vtkLookupTable::New();
  table->SetNumberOfTableValues(n_colors);
 
  for (i=0;i<n_colors;i++) {  

    eta = (i*1.0)/(n_colors-1);

    r_i = pow((float)(1.0-fabs(eta-1.0)),color_exp);
    g_i = pow((float)(1.0-2*fabs(eta-0.5)),color_exp);
    b_i = pow((float)(1.0-fabs(eta-0.0)),color_exp);

    table->SetTableValue(i,r_i,g_i,b_i,1.0);

  }

  table->Build();

  //-------------------------------------------------------------

  ofstream out("incl");

  for (i_time=t1; i_time<=t2; i_time++) {

    //-------------------------------------------------------------
    // Cap at phi = 0.0
    //
    vtkStructuredGrid *sgrid1 = vtkStructuredGrid::New();
    sgrid1->SetDimensions(n_r,n_theta,1);

    vtkPoints *points1 = vtkPoints::New();
    points1->Allocate(n_r*n_theta);
  
    vtkFloatArray *scalars1 = vtkFloatArray::New();
    scalars1->SetNumberOfComponents(1);
    scalars1->SetNumberOfTuples(n_r*n_theta);

    offset = 0;
    phi = 0.0;

    for (i_theta=0; i_theta<n_theta; i_theta++)
      {
	theta = -pi+i_theta*d_theta;
	for (i_r=0; i_r<n_r; i_r++) 
	  {

	    func(r[i_r],theta,phi,&s,x);

	    points1->InsertPoint(offset,x);
	    scalars1->InsertTuple(offset,&s);

	    offset++;
	  }
      }

    sgrid1->SetPoints(points1);
    sgrid1->GetPointData()->SetScalars(scalars1);

    vtkStructuredGridGeometryFilter *contour1 = 
      vtkStructuredGridGeometryFilter::New();
    contour1->SetInput(sgrid1);

    vtkPolyDataMapper *mapper1 = vtkPolyDataMapper::New();
    mapper1->SetInput(contour1->GetOutput());
    mapper1->SetLookupTable(table);

    vtkActor *actor1 = vtkActor::New();
    actor1->SetMapper(mapper1);
    //
    //------------------------------------------------------------

    //---------------------------------------------------
    // Rendering stage:
    //
    vtkRenderer *ren = vtkRenderer::New();
    ren->SetBackground(1,1,1);
    ren->AddActor(actor1);

    vtkRenderWindow *win = vtkRenderWindow::New();
    win->AddRenderer(ren);
    win->SetSize(xsize,ysize);

    vtkCamera *cam = vtkCamera::New();
    cam->SetViewUp(0,0,1);
    cam->SetPosition(x_pos,y_pos,0.0);
    cam->SetFocalPoint(0.0,y_pos,0.0);

    ren->SetActiveCamera(cam);
 
    win->Render();

    // Write Image

    vtkRendererSource *source = vtkRendererSource::New();
    source->SetInput(ren);

    //vtkRenderLargeImage *source = vtkRenderLargeImage::New();
    //source->SetInput(ren);
    //source->SetMagnification(1);

    make_name(i_time,filename);

    out << filename << endl;

    vtkJPEGWriter *image = vtkJPEGWriter::New();
    image->SetFileName(filename);
    image->SetInput(source->GetOutput());
    image->Write();

    // Cleanup (IMPORTANT)

    points1 ->Delete();
    scalars1->Delete();
    contour1->Delete();
    sgrid1  ->Delete();
    mapper1 ->Delete();
    actor1  ->Delete();

    cam   ->Delete();
    source->Delete();
    ren   ->Delete();

    image->Delete();

    win   ->Delete();

  }

  table ->Delete();
}

