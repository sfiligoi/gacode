/*==========================================================
  module: make_name
 
  Generate a filename based on integer value i.
  ==========================================================*/

void make_name(int i, char *filename) 
{

  // Filename root:

  char root[13] = "frame"; 

  // Digits:

  char *num[] = {"0","1","2","3","4","5","6","7","8","9"};

  int i1;
  int i2;
  int i3;

  i3 = i/100;
  i2 = (i-100*i3)/10;
  i1 = i-100*i3-10*i2;  

  strncat(root,num[i3],1); // filei
  strncat(root,num[i2],1); // fileij
  strncat(root,num[i1],1); // fileijk
  strncat(root,".jpg",5);  // fileijk.jpg

  for (int n=0; n<13; n++) filename[n] = root[n];

  cout << filename << endl;
}

