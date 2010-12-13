/*==========================================================
  module: make_name
 
  Generate a filename based on integer value i.
  ==========================================================*/

void make_name(int i, char *filename) 
{

  // Filename root:

  char root[14] = "frame"; 

  // Digits:

  char *num[] = {"0","1","2","3","4","5","6","7","8","9"};

  int i1;
  int i2;
  int i3;
  int i4;

  i4 = i/1000;
  i3 = (i-1000*i4)/100;
  i2 = (i-1000*i4-100*i3)/10;
  i1 = i-1000*i4-100*i3-10*i2;  

  strncat(root,num[i4],1); // framei
  strncat(root,num[i3],1); // frameij
  strncat(root,num[i2],1); // frameijk
  strncat(root,num[i1],1); // frameijkl
  strncat(root,".jpg",5);  // frameijkl.jpg

  for (int n=0; n<14; n++) filename[n] = root[n];

  cout << filename << endl;
}

