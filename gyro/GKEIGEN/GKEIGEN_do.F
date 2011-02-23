!-----------------------------------------------------
! GKEIGEN_do.F
!
! PURPOSE:
!  Determine the linear time-derivative matrix 
!  time_derivative_matrix and find the gkeigen_n_values
!  eigenvectors and eigenvalues satisfying the condition
!  specified in GKEIGEN_METHOD.
!  Only works in linear mode.
!-----------------------------------------------------

      subroutine GKEIGEN_do

      use gyro_globals
      use gyro_pointers
      use GKEIGEN_globals
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#include "finclude/petsc.h"
#include "finclude/petscvec.h"
#include "finclude/petscmat.h"
#include "finclude/petscis.h"
#include "finclude/petscsys.h"
#include "finclude/slepc.h"
#include "finclude/slepceps.h"
!--------------------------------------------------
! PETSc, SLEPc, and other local variable declarations. These variables
! are called in PETSc/SLEPc functions so are declared this way to
! conform to the specific PETSc/SLEPc build being used.
!
      Integer        n_converged,                                       &
     &               n_iter,                                            &
     &               irows(h_length_block),                             &
     &               icols(h_width_loc),                                &
     &               mrow,                                              &
     &               ncol

      Vec            dummy_vec

      PetscScalar    dummy_scalar,                                      &
     &               one_col(h_length_loc),                             &
     &               M_seq(seq_length)

      PetscReal      residual

      Mat                                                               &
     &               time_derivative_matrix

      PetscErrorCode ierr

      EPS            epsolve
!
!-------------------------------------------------

      integer :: irestart
      complex, dimension(500,h_length_loc) :: restart_cols
!      complex, dimension(0:h_length_loc-1,0:h_width_loc-1) :: M_loc
      complex, dimension(:,:), allocatable :: M_loc
      character (len=80) :: gkeigen_exit_message

10    Format(I2, ':', I2, ':', I2, ':' I3)
20    Format('Mode: ', I3, '    Omega: ', es18.9, '    Gamma: ', es18.9,&
     &       '    Residulal: ', es18.9) 
30    Format(es20.9, es20.9, es20.9) 
40    Format(I2, ',', I2, ' owns rows', I5, '-',I5, ' and columns ', I5,&
     &       '-', I5, ' with h_length_loc = ', I5)
50    Format(I2, ',', I2, ' jstate_start = ', I5, ' Local jsets: ',     &
     &       I5, '-', I5)
60    Format('Subset 1: ', I5, ' Subset 2: ', I5, ' Send from: ', I5,   &
     &    ' Send to: ', I5)
70    Format('Set: ', I2, ',', ' Proc: ', I2, ',',' Loc: ', I4, ',',    &
     &       ' Glob: ',I4,',',' Val: ',g12.5, ' + ', g12.5, ' i')
80    Format(I5, ' eigenvalues found')


      l_print = ( (gkeigen_j_set == 0) .AND. (i_proc == 0) )

      allocate(M_loc(0:h_length_loc-1,0:h_width_loc-1))     


!------------------------------------------------
! Make flags and inputs self-consistant.
!
      If (gkeigen_matrixonly == 1) gkeigen_mwrite_flag = 1

!-- Temporarily disable the seldom used matrix reading/writing.
!
      If (gkeigen_mwrite_flag == 1) call catch_error('Matrix writing tem&
     &porarily disabled in GKEIGEN')
      eigensolve_restart_flag = 0
!---------------------------------------

      If (gkeigen_n_values > h_length) Then
        call catch_error('INCONSISTENT: gkeigen_n_values must be less ',&
     &'than or equal to matrix dimension ', h_length, '.')
      Else If (gkeigen_kspace_dim < gkeigen_n_values) Then
        call catch_error('INCONSISTENT: gkeigen_kspace_dim must be      &
     &greater than or equal to gkeigen_n_values.')  
      Else If (gkeigen_kspace_dim > h_length) Then
        gkeigen_kspace_dim = h_length
        If (l_print) print *, 'Input gkeigen_kspace_dim must be ',      &
     &    'less than or equal to matrix dimension ', h_legnth, '.'
        If (l_print) print *, 'Reset gkeigen_kspace_dim = ', h_length
        call send_line('Reset gkeigen_kspace_dim = ', h_length)
      EndIf

!--------------------------------------------------


!-------------------------------------------------
! Initialize the SLEPc operating environment,
! placeholder vector dummy_vec, and matrix
! time_derivative_matrix.
!
      call SlepcInitialize(PETSC_NULL_CHARACTER,ierr)
!
      call VecCreate(MPI_COMM_WORLD,dummy_vec,ierr)
      call VecSetSizes(dummy_vec,PETSC_DECIDE,h_length,ierr)
      call VecSetFromOptions(dummy_vec,ierr)
!
!
      call MatCreate(MPI_COMM_WORLD,time_derivative_matrix,ierr)
      call MatSetSizes(time_derivative_matrix,PETSC_DECIDE,             &
     &                 PETSC_DECIDE,h_length,h_length,ierr)
      call MatSetType(time_derivative_matrix,MATDENSE,ierr)
!-------------------------------------------------


!------------------------------------------------
! Get all elements of time_derivative_matrix column
! by column with get_RHS.  Grid points owned by each
! processor are grouped together in the total state
! vector that defines time_derivative_matrix. 
!

      call MatGetOwnershipRange(time_derivative_matrix,jstate_start,    &
     &                          jstate_end,ierr)
      call MatGetSize(time_derivative_matrix,mrow,ncol,ierr) 

!      Print *, i_proc, 'PETSc owns ', jstate_start, jstate_end,         &
!     &         'size: ', mrow, ' x ', ncol

      istate_start = gkeigen_j_set * h_width_loc
      istate_end   = istate_start + h_width_loc - 1

      jstate_start = 0
      Do i_proc_e = 0, i_proc - 1
         jstate_start = jstate_start + n_nek_1 / n_proc_1
         If ( i_proc_e < Mod(n_nek_1,n_proc_1) )                        &
     &                  jstate_start = jstate_start + 1
      EndDo
      jstate_start = jstate_start * n_kinetic * n_x * n_stack
      jstate_end = jstate_start + h_length_loc - 1

!      Do jelem = 1, h_length_loc
!        irows(jelem) = jstate_start + jelem - 1
!      EndDo

      call date_and_time(values=time_array)
      If (l_print) Then
       Print 10, time_array(5),time_array(6),time_array(7),time_array(8)
       Print *, 'Starting value assignments.'
      EndIf

      If (eigensolve_restart_flag == 1) Then
        call GKEIGEN_matrix_read(n_last_col)
        If (l_print) print *, " Restarting from column ", n_last_col
      Else
        n_last_col = 0
        If ( (l_print) .AND. (gkeigen_mwrite_flag==1) ) Then
          open(unit=1,file=trim(path)//file_eigen_restart,              &
     &         status='replace')
          write (1,*) "!  Matrix elements for eigensolver."
          close(unit=1)
        EndIf
      EndIf


      restart_cols(:,:) = (0.0,0.0)
      istate = 0

      Do i_proc_e = 0, n_proc_1 - 1

       Do is1 = 1, n_kinetic

        p_nek_loc1 = 0
        Do p_nek1 = 1+i_proc_e, n_nek_1, n_proc_1

         If (i_proc_e == i_proc) p_nek_loc1 = p_nek_loc1 + 1
         Do i1 = 1, n_x

          Do m1 = 1, n_stack

           l_in_i_block = ((istate >= istate_start) .AND.               &
     &                   (istate <= istate_end))
!           If ((eigensolve_restart_flag == 1) .AND.                     &
!     &                 (istate < n_last_col)) Then
!
!                   one_col(1:h_length_loc) =                            &
!     &                       eigensolve_matrix(1:h_length_loc,istate+1)
!
           If (l_in_i_block) Then

             h(:,:,:,:) = (0.0,0.0)  
             If (i_proc_e == i_proc)                                    &
     &                h(m1,i1,p_nek_loc1,is1) = (1.0,0.0)
               call get_field_explicit
               call gyro_rhs_total
               ielem = istate - istate_start
               jstate = jstate_start
               Do is2 = 1, n_kinetic
                p_nek_loc2 = 0
                Do p_nek2 = 1+i_proc, n_nek_1, n_proc_1
                 p_nek_loc2 = p_nek_loc2 + 1
                 Do i2 = 1, n_x
                  Do m2 = 1, n_stack
                    jelem = jstate-jstate_start
                    one_col(jelem+1) = RHS(m2,i2,p_nek_loc2,is2)
                    M_loc(jelem,ielem)= RHS(m2,i2,p_nek_loc2,is2)  
!                    M_loc(jelem,ielem)=(1.,0.)*jelem+(0.,1.)*ielem
!                    if (l_print) print *, M_loc(jelem,ielem)
                    jstate = jstate + 1
                  EndDo ! m2
                 EndDo ! i2
                EndDo ! p_nek2
               EndDo ! is2

             EndIf

!             If (gkeigen_matrixonly == 0)                               &
!     &        call MatSetValues(time_derivative_matrix,h_length_loc,    &
!     &                     irows,1,istate,one_col,INSERT_VALUES,ierr)
             istate = istate + 1
!             irestart = istate-n_last_col
!             If (irestart >= 1) restart_cols(irestart,:)=one_col(:)
!             If (irestart == 500) Then
!               If (gkeigen_mwrite_flag==1)                              &
!     &           call GKEIGEN_matrix_write(restart_cols,n_last_col)
!               n_last_col = istate
!               restart_cols(:,:) = (0.0,0.0)
!             EndIf

          EndDo ! m1

         EndDo ! i1

        EndDo ! p_nek1

       EndDo ! is1

      EndDo ! i_proc_e
!
!---------------------------------------------

      If (gkeigen_mwrite_flag==1)                                       &
     &       call GKEIGEN_matrix_write(restart_cols,n_last_col)

      If (gkeigen_matrixonly == 1) Then
        call date_and_time(values=time_array)
        If (l_print) Then
          Print 10, time_array(5),time_array(6),time_array(7),          &
     &              time_array(8)
          Print *, 'Finished matrix construction without solving.'
        EndIf
        Goto 1000
      EndIf

      call date_and_time(values=time_array)
      If (l_print) Then
       Print 10, time_array(5),time_array(6),time_array(7),time_array(8)
       Print *, 'Starting message passing phase.'
      EndIf


!--------------------------------------------------------------------
! Matrix elements are passed from the processor calculating GKEIGEN
! processor to the processor that owns that element in PETSC, then
! assigned to the PETSc object time_derivative_matrix.

      Do j_sr = 0, gkeigen_proc_mult-1

        gkeigen_j_set_e0 = Mod(j_sr-gkeigen_j_set,gkeigen_proc_mult)
        If (gkeigen_j_set_e0 < 0)                                       &
     &         gkeigen_j_set_e0 = gkeigen_j_set_e0 + gkeigen_proc_mult 

        Do i_sr = 0, 1

          select case (i_sr)

            case (0)
              gkeigen_j_set_e1=Min(gkeigen_j_set_e0,gkeigen_j_set)
              gkeigen_j_set_e2=Max(gkeigen_j_set_e0,gkeigen_j_set)
            case (1)
              gkeigen_j_set_e1=Max(gkeigen_j_set_e0,gkeigen_j_set)
              gkeigen_j_set_e2=Min(gkeigen_j_set_e0,gkeigen_j_set)

          end select

          Do ielem = 0, h_width_loc-1
            icols(ielem+1) = gkeigen_j_set_e1*h_width_loc + ielem
          EndDo

          jelem_start = gkeigen_j_set_e2*h_length_block
          jelem_end = jelem_start + h_length_block - 1
          Do jelem = jelem_start, jelem_end
            i_index = jelem - jelem_start + 1
            irows(i_index) = jstate_start + jelem
          EndDo

          j_send_proc = i_proc*gkeigen_proc_mult + gkeigen_j_set_e1
          j_recv_proc = i_proc*gkeigen_proc_mult + gkeigen_j_set_e2

          l_send = (j_proc_tot==j_send_proc)
          l_recv = (j_proc_tot==j_recv_proc)
          l_set =  (l_send .AND. l_recv)

          If ((l_set==.TRUE.).AND.(i_sr==1)) Then
            l_set  = .FALSE.
            l_send = .FALSE.
            l_recv = .FALSE.
          EndIf

         If (l_send) Then
          Do ielem = 0, h_width_loc-1
            Do jelem = jelem_start, jelem_end
              iseq = (jelem-jelem_start)*h_width_loc + ielem + 1
              M_seq(iseq) = M_loc(jelem,ielem)
!              if (l_print) print *, iseq, ielem, M_loc(jelem,ielem)
            EndDo
          EndDo
         EndIf


         If (l_set) Then

!          print *, 'Setting on ', j_proc_tot
          call MatSetValues(time_derivative_matrix,h_length_block,irows,&
     &        h_width_loc,icols,M_seq,INSERT_VALUES,ierr)

         Else If (l_send) Then

!          print *, 'Sending from ', j_proc_tot, ' to', j_recv_proc,     &
!     &              size(M_seq)
          call MPI_SEND(M_seq,size(M_seq),MPIU_SCALAR,                  &
     &                  j_recv_proc,j_send_proc,MPI_COMM_WORLD,ierr)
!          print *, 'Send confirmed', j_proc_tot, ' to', j_recv_proc

         Else If (l_recv) Then

!          print *, 'Receiving on ', j_proc_tot, ' from', j_send_proc,   &
!     &              size(M_seq)
          call MPI_RECV(M_seq,size(M_seq),MPIU_SCALAR,                  &
     &                  j_send_proc,j_send_proc,MPI_COMM_WORLD,ierr)
!          print *, 'Receipt confirmed', j_proc_tot, ' from', j_send_proc
          call MatSetValues(time_derivative_matrix,h_length_block,irows,&
     &        h_width_loc,icols,M_seq,INSERT_VALUES,ierr)

         EndIf

         call MPI_barrier(MPI_COMM_WORLD,ierr)
!         if (l_print) print *, '*************************************'
        EndDo  ! gkeigen_j_set_e2

      EndDo  ! gkeigen_j_set_e1
!
!------------------------------------------------------------------------

      deallocate(M_loc)

      call date_and_time(values=time_array)
      If (l_print) Then
       Print 10, time_array(5),time_array(6),time_array(7),time_array(8)
       Print *, 'Starting matrix assembly.'
      EndIf
      call MatAssemblyBegin(time_derivative_matrix,MAT_FINAL_ASSEMBLY,  &
     &                      ierr)
      call MatAssemblyEnd(time_derivative_matrix,MAT_FINAL_ASSEMBLY,    &
     &                    ierr)

!      call MatView(time_derivative_matrix,PETSC_VIEWER_STDOUT_WORLD,    &
!     &             ierr)

!--------------------------------------------------------------------
! This creates the EigenProbelmSolver (EPS) framework and solves for the 
! gkeigen_n_values eigenvalues and vectors with the largest real part using a 
! Krylov subspace of maximum dimension gkeigen_kspace_dim.

!      gkeigen_kspace_dim = Max ( gkeigen_n_values+1,gkeigen_kspace_dim )
!      gkeigen_kspace_dim = Min ( gkeigen_kspace_dim,Int(h_length)   )

      call date_and_time(values=time_array)
      If (l_print) Then
       Print 10, time_array(5),time_array(6),time_array(7),time_array(8)
       Print *, 'Starting eigensolver.'
       Print *, "Using Krylov subspace of dimension ",gkeigen_kspace_dim
      EndIf
      call EPSCreate(MPI_COMM_WORLD,epsolve,ierr)
      call EPSSetOperators(epsolve,time_derivative_matrix,              &
     &                     PETSC_NULL_OBJECT,ierr)

      call EPSSetDimensions(epsolve,gkeigen_n_values,                   &
     &                      gkeigen_kspace_dim,h_length,ierr)

!--------------------------------------------------------------
! Search for eigenvalues in the region of the complex plane
! specified by gkeigen_method. GA frequency sign convention
! is opposite to sign of the imagainary part returned by SLEPc.
! GA convention is used in descriptions.
!
! (1) Largest real part / Largest growth rate
! (2) Largest imaginary part / Largest frequency in ion direction
! (3) Smallest imaginary part / Largest frequency in electron direction
! (4) Smallest real part / Most stable modes
!
!
      select case (gkeigen_method)

      case (1)

        call send_line(' Seeking largest growth rates.')
        If (l_print) print *, ' Seeking largest growth rates.'
        call EPSSetWhichEigenpairs(epsolve,EPS_LARGEST_REAL,ierr)

      case (2)

        call send_line(' Seeking largest negative frequency.')
        If (l_print) print *, ' Seeking largest negative frequency.'
        call EPSSetWhichEigenpairs(epsolve,EPS_LARGEST_IMAGINARY,ierr)

      case (3)

        call send_line(' Seeking largest positive frequency.')
        If (l_print) print *, ' Seeking largest positive frequency.'
        call EPSSetWhichEigenpairs(epsolve,EPS_SMALLEST_IMAGINARY,ierr)

      case (4)

        call send_line(' Seeking smallest growth rates.')
        If (l_print) print *, ' Seeking smallest growth rates.'
        call EPSSetWhichEigenpairs(epsolve,EPS_SMALLEST_REAL,ierr)

      case default

        call catch_error('INVALID: gkeigen_method')

      end select
!--------------------------------------------------------------------

      call EPSSetTolerances(epsolve,gkeigen_tol,gkeigen_iter,ierr)
      call EPSSetFromOptions(epsolve,ierr)
      call EPSSolve(epsolve,ierr)

      call date_and_time(values=time_array)
      If (l_print) Then
       Print 10, time_array(5),time_array(6),time_array(7),time_array(8)
       Print *, 'Finished eigensolve.'
      EndIf
!-------------------------------------------------------------------

!----------------------------------------------------------------
! All converged solutions are returned. Eigenvalues are stored
! in eigensolve_lambda, vectors in eigensolve_vector.

      call EPSGetConverged(epsolve,n_converged,ierr)
      If (l_print) Print *, n_converged, ' converged solution(s)'
      open(unit=1,file=trim(path)//file_eigen_freq,status='replace')
      step = 0
      Do iev = 0, n_converged-1
        call EPSGetEigenpair(epsolve,iev,dummy_scalar,PETSC_NULL,       &
     &          dummy_vec,PETSC_NULL,ierr)
        call EPSComputeResidualNorm(epsolve,iev,residual,ierr)
        ev_omega = -aimag(dummy_scalar)
        ev_gamma = real(dummy_scalar)
        If (l_print) Then
          print 20, iev+1, ev_omega, ev_gamma, residual
          write (1,30) ev_omega, ev_gamma, residual
          call gyro_write_precision(10,abs(ev_omega))
          step = step + 1
          call gyro_write_precision(10,abs(ev_gamma))
        EndIf
      EndDo
      close(unit=1)
!---------------------------------------------------------------------

      If (l_print) Then
        call EPSGetIterationNumber(epsolve,n_iter,ierr)
        Print *, n_iter, ' iterations performed'
      EndIf

      write(gkeigen_exit_message,80) n_converged
      call set_exit_status(gkeigen_exit_message,2)

      call VecDestroy(dummy_vec,ierr)
      call MatDestroy(time_derivative_matrix,ierr)
      call EPSDestroy(epsolve,ierr)

      call SlepcFinalize(ierr)

!      deallocate(M_loc) 

1000  end subroutine

