c -*-f90-*-
c
c ----------------------------------------------------------------------
c --- dimensional parameters for TORAY/GAFIT
c ----------------------------------------------------------------------
c
      integer nxmx, nzmx, nzmx2
      integer kj, kjp
c
c... 65x65 and 51 grid pts
c
c      parameter (nxmx = 65, nzmx = 65, nzmx2 = 2*nzmx)
c      parameter (kj = 52, kjp = 3*kj + 1)
c
c... 129x129 and 201 grid pts
c
      parameter (nxmx = 129, nzmx = 129, nzmx2 = 2*nzmx)
      parameter (kj = 202, kjp = 3*kj + 1)
c
      integer    nraymx      , ncntres
      parameter (nraymx = 100, ncntres = 500)
c
      integer    noutmx      , noutmxp
      parameter (noutmx =  11, noutmxp = noutmx + 1)
c
      integer    ndiag      , nvec     , nvtok     , nspace
      parameter (ndiag =  10, nvec = 22, nvtok = 10, nspace = 29)
c
c... GAFIT stuff
c
      integer    nw, nh, nh2
      parameter (nw = nxmx, nh = nzmx, nh2 = nzmx2)
      integer    npsimax, njmax, njmaxp
      integer    ncntra, nlimtra
      parameter (ncntra = 2001, nlimtra = 201)
c
c... 65x65 and 51 grid pts
c
c      parameter (npsimax = 65)
c      parameter (njmax = 51, njmaxp = njmax + 1)
c
c... 129x129 and 201 grid pts
c
      parameter (npsimax = 129)
      parameter (njmax = 201, njmaxp = njmax + 1)

