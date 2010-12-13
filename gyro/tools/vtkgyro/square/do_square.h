/*==========================================================
  module: do_square.h
  ==========================================================*/

void do_square( void ) {
  
  char filename[13];
  int i;
  int j;
  int k;
  int offset;

  float s;
  float x[3];

  if (n_n == 1) {
    d_alpha = (2*pi/n[0])/(n_alpha-1);
  } else {
    d_alpha = (2*pi/n[1])/(n_alpha-1);
  }

  d_r = (r[n_r-1]-r[0])/(n_rp-1);

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

  i_n_0 = 0;
  
  string file;
  
  //file = path+simdir+"/u.out"; 

  file = path+simdir+"/moment_n.out"; 
  n_field=1;
  i_field=1;

  ifstream in(file.c_str());
 
  if (!in) {
    cout << "Cannot open " << file << endl;
  }

  //-----------------------------------------------------------
  // fields
  //
  n_elements = n_theta_plot*n_r*n_field*n_n;

  cout << "** Reading " << file << "**" << endl;

  f_real = new float[n_elements];
  f_imag = new float[n_elements];


  //--------------------------------------------------------------

  for (i_time=0; i_time<=n_time; i_time++) {

    for (int i=0; i<n_elements; i++) {
      in >> f_real[i] >> f_imag[i];
    }

    if (i_time < t1) {

      read_progress(i_time,t1-1);

    } else {

      vtkStructuredGrid *sgrid1 = vtkStructuredGrid::New();
      sgrid1->SetDimensions(n_rp,n_alpha,1);

      vtkPoints *points1 = vtkPoints::New();
      points1->Allocate(n_rp*n_alpha);
  
      vtkFloatArray *scalars1 = vtkFloatArray::New();
      scalars1->SetNumberOfComponents(1);
      scalars1->SetNumberOfTuples(n_rp*n_alpha);

      offset = 0;

      for (i_alpha=0; i_alpha<n_alpha; i_alpha++)
	{
	  alpha = i_alpha*d_alpha;
	  for (i_r=0; i_r<n_rp; i_r++) 
	    {
              rp = r[0]+i_r*d_r;
	      func_sq(rp,alpha,&s,x);

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

      //=======================================================

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
      cam->SetPosition(x_pos,y_pos,z_pos);
      cam->SetFocalPoint(x_foc,y_foc,z_foc);

      ren->SetActiveCamera(cam);
 
      win->Render();

      // Write Image

      vtkRendererSource *source = vtkRendererSource::New();
      source->SetInput(ren);

      make_name(i_time,filename);

      out << filename << endl;

      vtkJPEGWriter *image = vtkJPEGWriter::New();
      image->ProgressiveOff();
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
  }

  table ->Delete();

  in.close();

}

