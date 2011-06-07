      subroutine part(work,largs,strnxt)
c******************************************************
c*****part finds the next item in the string.
c*****Input:
c*****work-a character variable.
c*****largs-the number of characters to be parsed
c*****Output:
c*****strnxt-the next item in the list,
c*****       delimited by matched blanks, "", or {}
c*****Note: On return, work and largs have been modified to reflect
c*****      the removal of the first item.
c*****last revision: 12/94 s.e.attenberger, w.a.houlberg, ornl
c*******************************************************
      implicit none
c
      integer largs, ltst, lenst, ntst
      character*(*) work,  strnxt
      character*1 find,    tab,    wtst
      character*1 space,      quote,      lbr,     rbr
      data        space/' '/, quote/'"'/, lbr/'{'/, rbr/'}'/
      tab=char(9)
c
      strnxt=' '
      if(largs.eq.0) return
c*****look for start of item (non-blank character)
      ltst=0
100   ltst=ltst+1
      if    ((work(ltst:ltst).eq.space.or.work(ltst:ltst).eq.tab)
     >        .and. ltst.lt.largs) then
        go to 100
      elseif((work(ltst:ltst).eq.space.or.work(ltst:ltst).eq.tab)
     >        .and. ltst.eq.largs) then
c*****  no items in this list
        strnxt=' '
        work=' '
        largs=0
        return
      elseif(work(ltst:ltst).eq.quote) then
        find=quote
        ltst=ltst+1
      elseif(work(ltst:ltst).eq.lbr) then
        find=rbr
        ltst=ltst+1
      else
        find=space
      endif
      lenst=1
      strnxt(lenst:lenst)=work(ltst:ltst)
c*****Start of item is character ltst.  Now search for end of item.
      ntst=ltst
200   ntst=ntst+1
c*****Treat tab like space for testing, but dont replace tab in work.
      wtst=work(ntst:ntst)
      if(wtst.eq.tab) wtst=space
      if(ntst.gt.largs.and.find.ne.space) then
        write(*,*)' fatal error, missing ///', find,'/// near'
        write(*,*)work(1:largs)
        stop
      elseif(ntst.le.largs.and.wtst.ne.find) then
        lenst=lenst+1
        strnxt(lenst:lenst)=work(ntst:ntst)
        go to 200
      elseif(ntst.gt.largs.and.find.eq.space) then
c*****  successful exit, end of string.
c*****  (no space delimiter is required at the end of a string)
        work(1:largs)=' '
        largs=0
      elseif(wtst.eq.find)then
c*****  successful exit
        work(1:largs-ntst)=work(ntst+1:largs)
        work(largs-ntst+1:largs)=' '
        largs=largs-ntst
      endif
c
      return
      end
