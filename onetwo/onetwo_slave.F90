!-------------------------------------------------------------------------------
!
SUBROUTINE onetwo_slave
!
!  jmp may 2006
!

   USE mpi12
   USE param, ONLY: kj
   IMPLICIT NONE
   INTEGER iflag_mpistop(1)
   
   DO
      CALL mpi12_recv_int(iflag_mpistop,1,99)
      if(iflag_mpistop(1).eq.1) then
      	print *,'final call at slave'
      	exit
      endif
      CALL glf_slave(kj)
   END DO
  
END SUBROUTINE  