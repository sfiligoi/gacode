# 7.19 Input/output <stdio.h>

cdef extern from *:
    ctypedef char const_char "const char"
    ctypedef void const_void "const void"

cdef extern from "stdio.h" nogil:

    ctypedef struct FILE
    cdef FILE *stdin
    cdef FILE *stdout
    cdef FILE *stderr
    
    enum: FOPEN_MAX
    enum: FILENAME_MAX
    FILE *fopen   (const_char *FILENAME, const_char  *OPENTYPE)
    FILE *freopen (const_char *FILENAME, const_char *OPENTYPE, FILE *STREAM)
    int  fclose   (FILE *STREAM)
    int  remove   (const_char *FILENAME)
    int  rename   (const_char *OLDNAME, const_char *NEWNAME)
    FILE *tmpfile ()

    enum: _IOFBF
    enum: _IOLBF
    enum: _IONBF
    int setvbuf (FILE *STREAM, char *BUF, int MODE, size_t SIZE)
    enum: BUFSIZ
    void setbuf (FILE *STREAM, char *BUF)

    size_t fread  (void *DATA, size_t SIZE, size_t COUNT, FILE *STREAM)
    size_t fwrite (const_void *DATA, size_t SIZE, size_t COUNT, FILE *STREAM)
    int    fflush (FILE *STREAM)

    enum: EOF
    int feof   (FILE *STREAM)
    int ferror (FILE *STREAM)

    enum: SEEK_SET
    enum: SEEK_CUR
    enum: SEEK_END
    int      fseek  (FILE *STREAM, long int OFFSET, int WHENCE)
    void     rewind (FILE *STREAM)
    long int ftell  (FILE *STREAM)

    ctypedef long long int fpos_t
    ctypedef fpos_t const_fpos_t "const fpos_t"
    int fgetpos (FILE *STREAM, fpos_t *POSITION)
    int fsetpos (FILE *STREAM, const_fpos_t *POSITION)

    int scanf    (const_char *TEMPLATE, ...)
    int sscanf   (const_char *S, const_char *TEMPLATE, ...)
    int fscanf   (FILE *STREAM, const_char *TEMPLATE, ...)

    int printf   (const_char *TEMPLATE, ...)
    int sprintf  (char *S, const_char *TEMPLATE, ...)
    int snprintf (char *S, size_t SIZE, const_char *TEMPLATE, ...)
    int fprintf  (FILE *STREAM, const_char *TEMPLATE, ...)

    void perror  (const_char *MESSAGE)

    char *gets  (char *S)
    char *fgets (char *S, int COUNT, FILE *STREAM)

    int  puts   (const_char *S)
    int  fputs  (const_char *S, FILE *STREAM)

# 7.20 General utilities <stdlib.h>

cdef extern from *:
    ctypedef char const_char "const char"
    ctypedef void const_void "const void"

cdef extern from "stdlib.h" nogil:

    # 7.20.1 Numeric conversion functions
    int atoi (const_char *STRING)
    long atol (const_char *STRING)
    long long atoll (const_char *STRING)
    double atof (const_char *STRING)
    long strtol (const_char *STRING, char **TAILPTR, int BASE)
    unsigned long int strtoul (const_char *STRING, char **TAILPTR, int BASE)
    long long int strtoll (const_char *STRING, char **TAILPTR, int BASE)
    unsigned long long int strtoull (const_char *STRING, char **TAILPTR, int BASE)
    float strtof (const_char *STRING, char **TAILPTR)
    double strtod (const_char *STRING, char **TAILPTR)
    long double strtold (const_char *STRING, char **TAILPTR)
    
    # 7.20.2 Pseudo-random sequence generation functions
    enum: RAND_MAX
    int rand ()
    void srand (unsigned int SEED)

    # 7.20.3 Memory management functions
    void *calloc (size_t COUNT, size_t ELTSIZE)
    void free (void *PTR)
    void *malloc (size_t SIZE)
    void *realloc (void *PTR, size_t NEWSIZE)

    # 7.20.4 Communication with the environment
    enum: EXIT_FAILURE
    enum: EXIT_SUCCESS
    void exit (int STATUS)
    void _Exit (int STATUS)
    int atexit (void (*FUNCTION) ())
    void abort ()
    char *getenv (const_char *NAME)
    int system (const_char *COMMAND)

    #7.20.5 Searching and sorting utilities
    void *bsearch (const_void *KEY, const_void *ARRAY,
                   size_t COUNT, size_t SIZE, 
                   int (*COMPARE)(const_void *, const_void *))
    void qsort (void *ARRAY, size_t COUNT, size_t SIZE,
                int (*COMPARE)(const_void *, const_void *))

    # 7.20.6 Integer arithmetic functions
    int abs (int NUMBER)
    long int labs (long int NUMBER)
    long long int llabs (long long int NUMBER)
    ctypedef struct div_t:
        int quot
        int rem
    div_t div (int NUMERATOR, int DENOMINATOR)
    ctypedef struct ldiv_t:
        long int quot
        long int rem
    ldiv_t ldiv (long int NUMERATOR, long int DENOMINATOR)
    ctypedef struct lldiv_t:
        long long int quot
        long long int rem
    lldiv_t lldiv (long long int NUMERATOR, long long int DENOMINATOR)

import numpy as np

def loadtxt(char* f, int dSize):
    cdef FILE *pFile
    cdef long lSize
    cdef char * buffer
    cdef size_t result
    cdef int i
    i=0

    pFile = fopen(f, "rb")
    if (pFile==NULL):
        fputs ("File error\n",stderr)
        exit (1)

    #obtain file size:
    fseek (pFile , 0 , SEEK_END)
    lSize = ftell (pFile)
    rewind (pFile)

    #allocate memory to contain the whole file:
    buffer = <char*> calloc(lSize, sizeof(char))
    if (buffer == NULL):
        fputs ("Memory error\n",stderr)
        exit (2)

    #Copy the file into the buffer:
    result = fread (buffer,1,lSize,pFile)
    if (result != lSize):
        fputs ("Reading error\n",stderr)
        exit (3)

    #The whole file is now loaded in the memory buffer.
    array = np.asarray(np.frombuffer(buffer, dtype='S' + str(dSize)), dtype=np.float32)

    #terminate
    fclose (pFile)
    free (buffer)
    return array