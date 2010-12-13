#include <iostream>
#include <fstream>
#include <string>
using namespace std;

#include "make_name.h"

int main( int argc, char *argv[] )
{
  
  char filename[13];
  int i_time;

  //-----------------------
  int t1 = 1;
  int t2 = 744;
  //-----------------------

  ofstream out("incl");

  for (i_time=t1; i_time<=t2; i_time++) {

    make_name(i_time,filename);
    out << filename << endl;
  }

}

