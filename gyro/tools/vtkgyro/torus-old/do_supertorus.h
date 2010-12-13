/*==========================================================
  module: do_supertorus.h
 
  Torus rendering with halo.

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

void do_supertorus( void ) {

  char filename[13];
  int i;
  int j;
  int k;
  int offset;

  float s;

  float x[3];

  i_buffer = 12;

  d_theta = 2*pi/(n_theta-1);
  d_phi   = phi_max*2*pi/(n_phi-1);

  n_rp = n_r-2*i_buffer;

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

  vtkLookupTable *halo_table = vtkLookupTable::New();
  halo_table->SetNumberOfTableValues(n_colors);
 
  for (i=0;i<n_colors;i++) {  

    eta = (i*1.0)/(n_colors-1);

    // Fixed light-blue color.

    halo_table->SetTableValue(i,0.25,1.0,1.0,0.35);

  }

  halo_table->Build();

  vtkCamera *cam = vtkCamera::New();

  vtkLight *light1 = vtkLight::New();
  light1->SetColor(1,1,1);
  light1->SetIntensity(intensity);
  light1->SetPosition(x_light_pos1,y_light_pos1,z_light_pos1);
  light1->SetFocalPoint(x_light_foc1,y_light_foc1,z_light_foc1);

  vtkLight *light2 = vtkLight::New();
  light2->SetColor(1,1,1);
  light2->SetIntensity(intensity);
  light2->SetPosition(x_light_pos2,y_light_pos2,z_light_pos2);
  light2->SetFocalPoint(x_light_foc2,y_light_foc2,z_light_foc2);

  vtkProperty *prop = vtkProperty::New();
  prop->SetSpecular(10.0);
  prop->SetSpecularPower(100);
  prop->SetDiffuse(0.2);
  prop->SetInterpolationToGouraud();
  //
  //-------------------------------------------------------------

  i_time  = t1;
  i_move  = 0;
  i_frame = 0;
  while (i_time < t2) {

    i_frame++;

    // Set the motion parameter, e0.

    if (i_time < t_stop) { 

      i_time++;
      e0 = 0.0;

    } else if (i_move < n_move) {

      i_move++;
      e0 = ((float) i_move/(n_move+1));

    } else {

      i_time++;
      e0 = 1.0;

    }

    //-------------------------------------------------------------
    // Cap at phi = 0.0
    //
    vtkStructuredGrid *sgrid1 = vtkStructuredGrid::New();
    sgrid1->SetDimensions(n_rp,n_theta,1);

    vtkPoints *points1 = vtkPoints::New();
    points1->Allocate(n_rp*n_theta);
  
    vtkFloatArray *scalars1 = vtkFloatArray::New();
    scalars1->SetNumberOfComponents(1);
    scalars1->SetNumberOfTuples(n_rp*n_theta);

    offset = 0;
    phi = 0.0;

    for (i_theta=0; i_theta<n_theta; i_theta++)
      {
	theta = -pi+i_theta*d_theta;
	for (i_r=i_buffer; i_r<n_r-i_buffer; i_r++) 
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
    actor1->SetProperty(prop);
    //
    //------------------------------------------------------------


    //-------------------------------------------------------------
    // Cap at phi = phi_max
    //
    vtkStructuredGrid *sgrid2 = vtkStructuredGrid::New();
    sgrid2->SetDimensions(n_rp,n_theta,1);

    vtkPoints *points2 = vtkPoints::New();
    points2->Allocate(n_rp*n_theta);
  
    vtkFloatArray *scalars2 = vtkFloatArray::New();
    scalars2->SetNumberOfComponents(1);
    scalars2->SetNumberOfTuples(n_rp*n_theta);

    offset = 0;
    phi = -2*pi*phi_max;

    for (i_theta=0; i_theta<n_theta; i_theta++)
      {
	theta = -pi+i_theta*d_theta;
	for (i_r=i_buffer; i_r<n_r-i_buffer; i_r++) 
	  {
	    func(r[i_r],theta,phi,&s,x);

	    points2->InsertPoint(offset,x);
	    scalars2->InsertTuple(offset,&s);

	    offset++;
	  }
      }

    sgrid2->SetPoints(points2);
    sgrid2->GetPointData()->SetScalars(scalars2);

    vtkStructuredGridGeometryFilter *contour2 = 
      vtkStructuredGridGeometryFilter::New();
    contour2->SetInput(sgrid2);

    vtkPolyDataMapper *mapper2 = vtkPolyDataMapper::New();
    mapper2->SetInput(contour2->GetOutput());
    mapper2->SetLookupTable(table);

    vtkActor *actor2 = vtkActor::New();
    actor2->SetMapper(mapper2);
    actor2->SetProperty(prop);
    //
    //-------------------------------------------------------------


    //-------------------------------------------------------------
    // Outer body of torus:
    //
    vtkStructuredGrid *sgrid3 = vtkStructuredGrid::New();
    sgrid3->SetDimensions(1,n_theta,n_phi);

    vtkPoints *points3 = vtkPoints::New();
    points3->Allocate(n_theta*n_phi);
  
    vtkFloatArray *scalars3 = vtkFloatArray::New();
    scalars3->SetNumberOfComponents(1);
    scalars3->SetNumberOfTuples(n_theta*n_phi);

    offset = 0;
    i_r = n_r-1-i_buffer;

    for (k=0; k<n_phi; k++)
      {
	phi = -k*d_phi;
	for (i_theta=0; i_theta<n_theta; i_theta++) 
	  {
	    theta = -pi+i_theta*d_theta;

	    func(r[i_r],theta,phi,&s,x);

	    points3->InsertPoint(offset,x);
	    scalars3->InsertTuple(offset,&s);

	    offset++;
	  }
      }

    sgrid3->SetPoints(points3);
    sgrid3->GetPointData()->SetScalars(scalars3);

    vtkStructuredGridGeometryFilter *contour3 = 
      vtkStructuredGridGeometryFilter::New();
    contour3->SetInput(sgrid3);

    vtkPolyDataMapper *mapper3 = vtkPolyDataMapper::New();
    mapper3->SetInput(contour3->GetOutput());
    mapper3->SetLookupTable(table);

    vtkActor *actor3 = vtkActor::New();
    actor3->SetMapper(mapper3);
    actor3->SetProperty(prop);
    //
    //---------------------------------------------------

    //-------------------------------------------------------------
    // Inner body of torus:
    //
    vtkStructuredGrid *sgrid4 = vtkStructuredGrid::New();
    sgrid4->SetDimensions(1,n_theta,n_phi);

    vtkPoints *points4 = vtkPoints::New();
    points4->Allocate(n_theta*n_phi);
  
    vtkFloatArray *scalars4 = vtkFloatArray::New();
    scalars4->SetNumberOfComponents(1);
    scalars4->SetNumberOfTuples(n_theta*n_phi);

    offset = 0;
    i_r = 0+i_buffer;

    for (k=0; k<n_phi; k++)
      {
	phi = -k*d_phi;
	for (i_theta=0; i_theta<n_theta; i_theta++) 
	  {
	    theta = -pi+i_theta*d_theta;

	    func(r[i_r],theta,phi,&s,x);

	    points4->InsertPoint(offset,x);
	    scalars4->InsertTuple(offset,&s);

	    offset++;
	  }
      }

    sgrid4->SetPoints(points4);
    sgrid4->GetPointData()->SetScalars(scalars4);

    vtkStructuredGridGeometryFilter *contour4 = 
      vtkStructuredGridGeometryFilter::New();
    contour4->SetInput(sgrid4);

    vtkPolyDataMapper *mapper4 = vtkPolyDataMapper::New();
    mapper4->SetInput(contour4->GetOutput());
    mapper4->SetLookupTable(table);

    vtkActor *actor4 = vtkActor::New();
    actor4->SetMapper(mapper4);
    actor4->SetProperty(prop);
    //
    //---------------------------------------------------

    //-------------------------------------------------------------
    // Torus halo:
    //
    d_theta = 2*pi/(n_theta/4-1);
    d_phi   = phi_max*2*pi/(n_phi/4-1);
    // 
    vtkStructuredGrid *halo_sgrid = vtkStructuredGrid::New();
    halo_sgrid->SetDimensions(1,n_theta/4,n_phi/4);

    vtkPoints *halo_points = vtkPoints::New();
    halo_points->Allocate(n_theta*n_phi/16);
  
    vtkFloatArray *halo_scalars = vtkFloatArray::New();
    halo_scalars->SetNumberOfComponents(1);
    halo_scalars->SetNumberOfTuples(n_theta*n_phi/16);

    offset = 0;
    i_r = n_r-1-i_buffer;

    for (k=0; k<n_phi/4; k++)
      {
	phi = -k*d_phi;
	for (i_theta=0; i_theta<n_theta/4; i_theta++) 
	  {
	    theta = -pi+i_theta*d_theta;

	    c_func(1.0,theta,phi,&s,x);

	    halo_points->InsertPoint(offset,x);
	    halo_scalars->InsertTuple(offset,&s);

	    offset++;
	  }
      }

    halo_sgrid->SetPoints(halo_points);
    halo_sgrid->GetPointData()->SetScalars(halo_scalars);

    vtkStructuredGridGeometryFilter *halo_contour = 
      vtkStructuredGridGeometryFilter::New();
    halo_contour->SetInput(halo_sgrid);

    vtkPolyDataMapper *halo_mapper = vtkPolyDataMapper::New();
    halo_mapper->SetInput(halo_contour->GetOutput());
    halo_mapper->SetLookupTable(halo_table);

    vtkActor *halo_actor = vtkActor::New();
    halo_actor->SetMapper(halo_mapper);
    halo_actor->SetProperty(prop);
    //
    //---------------------------------------------------

    make_name(i_frame,filename);

    //---------------------------------------------------
    // Rendering stage:
    //
    vtkRenderer *ren = vtkRenderer::New();

    vtkRenderWindow *win = vtkRenderWindow::New();
    win->AddRenderer(ren);
    win->SetSize(xsize,ysize);
 
    ren->SetBackground(1,1,1);
    ren->AddActor(actor1);
    ren->AddActor(actor2);
    ren->AddActor(actor3);
    ren->AddActor(actor4);
    //    ren->AddActor(halo_actor);

    e1 = 1.0-e0;
    e2 = 4.0*e0*e1;

    x_cam = x_pos*e1-y_pos*e0;
    y_cam = y_pos*e1-x_pos*e0;

    fx_cam = x_foc*e1-y_pos*e0;
    fy_cam = y_foc*e1-x_foc*e0;

    cam->SetViewUp(0,0,1);
    cam->SetPosition(x_cam,y_cam,z_pos*e2);
    cam->SetFocalPoint(fx_cam,fy_cam,z_foc*e2);
    cam->SetViewAngle(view_angle+dview_angle*e2);

    ren->SetActiveCamera(cam);
    ren->AddLight(light1);
    ren->AddLight(light2);
 
    win->Render();

    // Write Image

    //vtkRendererSource *source = vtkRendererSource::New();
    //source->SetInput(ren);

    vtkRenderLargeImage *source = vtkRenderLargeImage::New();
    source->SetInput(ren);
    source->SetMagnification(image_magnification);

    vtkJPEGWriter *image = vtkJPEGWriter::New();
    image->SetFileName(filename);
    image->SetInput(source->GetOutput());
    image->Write();
    image->Delete(); 

    // Cleanup (IMPORTANT)

    points1 ->Delete();
    scalars1->Delete();
    contour1->Delete();
    sgrid1  ->Delete();
    mapper1 ->Delete();
    actor1  ->Delete();

    points2 ->Delete();
    scalars2->Delete();
    contour2->Delete();
    sgrid2  ->Delete();
    mapper2 ->Delete();
    actor2  ->Delete();

    points3 ->Delete();
    scalars3->Delete();
    contour3->Delete();
    sgrid3  ->Delete();
    mapper3 ->Delete();
    actor3  ->Delete();
 
    points4 ->Delete();
    scalars4->Delete();
    contour4->Delete();
    sgrid4  ->Delete();
    mapper4 ->Delete();
    actor4  ->Delete();

    halo_points ->Delete();
    halo_scalars->Delete();
    halo_contour->Delete();
    halo_sgrid  ->Delete();
    halo_mapper ->Delete();
    halo_actor  ->Delete();

    source->Delete();
    ren   ->Delete();
    win   ->Delete();
  }

  table->Delete();
  halo_table->Delete();
  cam->Delete();
  prop->Delete();

}
