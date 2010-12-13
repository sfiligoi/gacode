/*==========================================================
  module: do_square.h
 
  Flux-tube square rendering.

  Geometry: (r,alpha)

  Caps:

  dim[0] = n_r;
  dim[1] = n_theta;
  dim[2] = 1;

  Ends:

  dim[0] = 1;
  dim[1] = n_theta;
  dim[2] = n_phi;

  ==========================================================*/

void do_fancy_2( void ) {
  
  char filename[13];
  int i;
  int j;
  int k;
  int offset;

  unsigned short s;
  float x[3];

  if (n_n == 1) {
    ly = 2*pi/n[0];
  } else {
    ly = 2*pi/n[1];
  }

  d_phi   = ly/(n_phi-1);
  d_theta = ly/(n_theta-1);

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
    // Square
    //
    vtkStructuredPoints *spt = vtkStructuredPoints::New();
    spt->SetScalarType(VTK_UNSIGNED_SHORT);
    spt->SetDimensions(n_r,n_phi,n_theta);
    spt->SetOrigin(0,0,0);
    spt->SetSpacing(1.0/n_r,1.0/n_phi,1.0/n_theta);
  
    vtkUnsignedShortArray *scalars = vtkUnsignedShortArray::New();
    scalars->SetNumberOfComponents(1);
    scalars->SetNumberOfTuples(n_r*n_phi*n_theta);

    offset = 0;

    for (i_theta=0; i_theta<n_theta; i_theta++)
      {
	theta = (i_theta-n_theta/2.0)*d_theta;
	for (i_phi=0; i_phi<n_phi; i_phi++)
	  {
	    phi = i_phi*d_phi;
	    for (i_r=0; i_r<n_r; i_r++) 
	      {
		func(r[i_r],phi,theta,&s,x);

		scalars->InsertTupleValue(offset,&s);
 
                //cout << s << endl;

		offset++;
	      }
	  }

      }

    spt->GetPointData()->SetScalars(scalars);
    scalars->Delete();

    vtkPiecewiseFunction *opac=vtkPiecewiseFunction::New();
    opac->AddPoint(64,1.0);
    opac->AddPoint(128,0.0);
    opac->AddPoint(192,1.0);

    vtkPiecewiseFunction *gopac=vtkPiecewiseFunction::New();
    gopac->AddPoint(0,1.0);
    gopac->AddPoint(64,0.0);
    gopac->AddPoint(128,0.0);
    gopac->AddPoint(192,0.0);
    gopac->AddPoint(255,1.0);

    vtkColorTransferFunction *tfunc=vtkColorTransferFunction::New();
    tfunc->AddRGBPoint(0,1.0,0.0,0.0);
    tfunc->AddRGBPoint(128,0.0,1.0,0.0);
    tfunc->AddRGBPoint(255,0.0,0.0,1.0);

    vtkVolumeProperty *volprop=vtkVolumeProperty::New();
    volprop->SetColor(tfunc);
    volprop->SetScalarOpacity(opac);
    //volprop->SetGradientOpacity(gopac);
    //volprop->ShadeOn();
    volprop->SetInterpolationTypeToLinear();

    vtkVolumeRayCastCompositeFunction *compo= 
      vtkVolumeRayCastCompositeFunction::New();

    vtkVolumeRayCastMapper *mapper=vtkVolumeRayCastMapper::New();
    mapper->SetVolumeRayCastFunction(compo);
    mapper->SetInput(spt);

    vtkVolume *vol=vtkVolume::New();
    vol->SetMapper(mapper);
    vol->SetProperty(volprop);
    //
    //------------------------------------------------------------

    //---------------------------------------------------
    // Rendering stage:
    //
    vtkRenderer *ren = vtkRenderer::New();
    ren->SetBackground(1,1,1);
    //ren->AddActor(actor1);
    ren->AddVolume(vol);

    vtkRenderWindow *win = vtkRenderWindow::New();
    win->AddRenderer(ren);
    win->SetSize(xsize,ysize);

    vtkCamera *cam = vtkCamera::New();
    cam->SetViewUp(0,0,1);
    cam->SetPosition(x_pos,y_pos,z_pos);
    cam->SetFocalPoint(x_foc,y_foc,z_foc);

    ren->SetActiveCamera(cam);
 
    win->Render();

    //
    //------------------------------------------------------------

   //---------------------------------------------------
    // Write Image

    vtkRendererSource *source = vtkRendererSource::New();
    source->SetInput(ren);

    make_name(i_time,filename);

    out << filename << endl;

    vtkJPEGWriter *image = vtkJPEGWriter::New();
    image->SetFileName(filename);
    image->SetInput(source->GetOutput());
    image->Write();

    // Cleanup (IMPORTANT)

    //points1 ->Delete();
    //contour1->Delete();
    //sgrid1  ->Delete();
    //mapper1 ->Delete();
    //actor1  ->Delete();

    //cam   ->Delete();
    source->Delete();
    ren   ->Delete();

    image->Delete();

    win   ->Delete();

  }

  table ->Delete();
}

