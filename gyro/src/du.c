/* Dummy function only to get the linker include these procedures.
   Didn't find a way to to this with a linker command line. */
extern int printf(),fopen(),fclose(),fwrite();
int f()
{
  printf();
  fopen();
  fclose();
  fwrite();
}
