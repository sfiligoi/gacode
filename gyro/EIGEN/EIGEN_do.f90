!-----------------------------------------------------
! EIGEN_do.F
!
! PURPOSE:
!  Determine the linear time-derivative matrix 
!  time_derivative_matrix and find the gkeigen_n_values
!  most unstable eigenvectors and eigenvalues.
!  Only works in linear mode.
!-----------------------------------------------------

      subroutine EIGEN_do

      use gyro_globals
      use gyro_pointers
      use EIGEN_globals
      use math_constants
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#include "finclude/petsc.h"
#include "finclude/petscvec.h"
#include "finclude/petscmat.h"
#include "finclude/petscis.h"
#include "finclude/petscsys.h"
#include "finclude/slepc.h"
#include "finclude/slepceps.h"

!--------------------------------------------------
! PETSc, SLEPc, and other local variable declarations.
!
      Integer        n_converged,jstate_start,jstate_end,n_iter,        &
     &               time_array(8),irestart,irows(h_length_loc),jelem,  &
     &               mrow,ncol
      Vec            dummy_vec
      PetscScalar    dummy_scalar,                                      &
     &               one_col(h_length_loc),                             &
     &               restart_cols(500,h_length_loc)
      PetscReal      residual,ev_omega,ev_gamma
      Mat                                                               &
     &               time_derivative_matrix
      PetscErrorCode ierr
      EPS            epsolve
!
!-------------------------------------------------

10    Format(I2, ':', I2, ':', I2, ':' I3)
20    Format('Mode: ', I3, '    Omega: ', es18.9, '    Gamma: ', es18.9,&
     &       '    Residulal: ', es18.9) 
30    Format(es20.9, es20.9, es20.9) 

!------------------------------------------------
! Make flags and inputs self-consistant.
!
      If (gkeigen_matrixonly == 1) gkeigen_mwrite_flag = 1

      If (gkeigen_n_values > h_length) Then
        call catch_error('INCONSISTENT: gkeigen_n_values must be less ',&
     &'than or equal to matrix dimension ', h_length, '.')
      Else If (gkeigen_kspace_dim < gkeigen_n_values) Then
        call catch_error('INCONSISTENT: gkeigen_kspace_dim must be      &
     &greater than or equal to gkeigen_n_values.')  
      Else If (gkeigen_kspace_dim > h_length) Then
        gkeigen_kspace_dim = h_length
        If (i_proc==0) print *, 'Input gkeigen_kspace_dim must be ',    &
     &    'less than or equal to matrix dimension ', h_legnth, '.'
        If (i_proc==0) print *, 'Reset gkeigen_kspace_dim = ', h_length
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
      call VecCreate(GYRO_COMM_WORLD,dummy_vec,ierr)
      call VecSetSizes(dummy_vec,PETSC_DECIDE,h_length,ierr)
      call VecSetFromOptions(dummy_vec,ierr)
!
!
      call MatCreate(GYRO_COMM_WORLD,time_derivative_matrix,ierr)
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

      jstate_start = 0
      Do i_proc_e = 0, i_proc - 1
         jstate_start = jstate_start + n_nek_1 / n_proc_1
         If ( i_proc_e < Mod(n_nek_1,n_proc_1) )                        &
     &                  jstate_start = jstate_start + 1
      EndDo
      jstate_start = jstate_start * n_kinetic * n_x * n_stack

      Do jelem = 1, h_length_loc
        irows(jelem) = jstate_start + jelem - 1
      EndDo

      call date_and_time(values=time_array)
      If (i_proc == 0) Then
       Print 10, time_array(5),time_array(6),time_array(7),time_array(8)
       Print *, 'Starting value assignments.'
      EndIf

      allocate(eigensolve_matrix(h_length_loc,h_length))
      eigensolve_matrix(:,:) = (0.0,0.0)

      If (eigensolve_restart_flag == 1) Then
        call EIGEN_matrix_read(n_last_col)
        If (i_proc==0) print *, " Restarting from column ", n_last_col
      Else
        n_last_col = 0
        If ( (i_proc==0) .AND. (gkeigen_mwrite_flag==1) ) Then
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

                 If ((eigensolve_restart_flag == 1) .AND.               &
     &                                   (istate < n_last_col)) Then

                   one_col(1:h_length_loc) =                            &
     &                       eigensolve_matrix(1:h_length_loc,istate+1)

                 Else

                  h(:,:,:,:) = (0.0,0.0)
                  If (i_proc_e == i_proc)                               &
     &                     h(m1,i1,p_nek_loc1,is1) = (1.0,0.0)
                  call get_field_explicit
                  call get_RHS
                  jstate = jstate_start
                  Do is2 = 1, n_kinetic
                     p_nek_loc2 = 0
                     Do p_nek2 = 1+i_proc, n_nek_1, n_proc_1
                        p_nek_loc2 = p_nek_loc2 + 1
                        Do i2 = 1, n_x
                           Do m2 = 1, n_stack
                              jelem = jstate-jstate_start+1
                              one_col(jelem) = RHS(m2,i2,p_nek_loc2,is2)
                              jstate = jstate + 1
                           EndDo ! m2
                        EndDo ! i2
                     EndDo ! p_nek2
                  EndDo ! is2

                 EndIf

                 If (gkeigen_matrixonly == 0)                           &
     &            call MatSetValues(time_derivative_matrix,h_length_loc,&
     &                     irows,1,istate,one_col,INSERT_VALUES,ierr)
                 istate = istate + 1
                 irestart = istate-n_last_col
                 If (irestart >= 1) restart_cols(irestart,:)=one_col(:)
                 If (irestart == 500) Then
                   If (gkeigen_mwrite_flag==1)                          &
     &               call EIGEN_matrix_write(restart_cols,n_last_col)
                   n_last_col = istate
                   restart_cols(:,:) = (0.0,0.0)
                 EndIf

               EndDo ! m1

            EndDo ! i1

         EndDo ! p_nek1

       EndDo ! is1

      EndDo ! i_proc_e
!
!---------------------------------------------

      If (gkeigen_mwrite_flag==1)                                       &
     &       call EIGEN_matrix_write(restart_cols,n_last_col)

      If (gkeigen_matrixonly == 1) Then
        call date_and_time(values=time_array)
        If (i_proc == 0) Then
          Print 10, time_array(5),time_array(6),time_array(7),          &
     &              time_array(8)
          Print *, 'Finished matrix construction without solving.'
        EndIf
        Goto 1000
      EndIf

      call date_and_time(values=time_array)
      If (i_proc == 0) Then
       Print 10, time_array(5),time_array(6),time_array(7),time_array(8)
       Print *, 'Starting matrix assembly.'
      EndIf
      call MatAssemblyBegin(time_derivative_matrix,MAT_FINAL_ASSEMBLY,  &
     &                      ierr)
      call MatAssemblyEnd(time_derivative_matrix,MAT_FINAL_ASSEMBLY,    &
     &                    ierr)


!--------------------------------------------------------------------
! This creates the EigenProbelmSolver (EPS) framework and solves for the 
! gkeigen_n_values eigenvalues and vectors with the largest real part using a 
! Krylov subspace of maximum dimension gkeigen_kspace_dim.

!      gkeigen_kspace_dim = Max ( gkeigen_n_values+1,gkeigen_kspace_dim )
!      gkeigen_kspace_dim = Min ( gkeigen_kspace_dim,Int(h_length)   )

      call date_and_time(values=time_array)
      If (i_proc == 0) Then
       Print 10, time_array(5),time_array(6),time_array(7),time_array(8)
       Print *, 'Starting eigensolver.'
       Print *, "Using Krylov subspace of dimension ",gkeigen_kspace_dim
      EndIf
      call EPSCreate(GYRO_COMM_WORLD,epsolve,ierr)
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
        If (i_proc==0) print *, ' Seeking largest growth rates.'
        call EPSSetWhichEigenpairs(epsolve,EPS_LARGEST_REAL,ierr)

      case (2)

        call send_line(' Seeking largest negative frequency.')
        If (i_proc==0) print *, ' Seeking largest negative frequency.'
        call EPSSetWhichEigenpairs(epsolve,EPS_LARGEST_IMAGINARY,ierr)

      case (3)

        call send_line(' Seeking largest positive frequency.')
        If (i_proc==0) print *, ' Seeking largest positive frequency.'
        call EPSSetWhichEigenpairs(epsolve,EPS_SMALLEST_IMAGINARY,ierr)

      case (4)

        call send_line(' Seeking smallest growth rates.')
        If (i_proc==0) print *, ' Seeking smallest growth rates.'
        call EPSSetWhichEigenpairs(epsolve,EPS_SMALLEST_REAL,ierr)

      case default

        call catch_error('INVALID: gkeigen_method')

      end select
!--------------------------------------------------------------------

      call EPSSetTolerances(epsolve,gkeigen_tol,gkeigen_iter,ierr)
      call EPSSetFromOptions(epsolve,ierr)
      call EPSSolve(epsolve,ierr)

      call date_and_time(values=time_array)
      If (i_proc == 0) Then
       Print 10, time_array(5),time_array(6),time_array(7),time_array(8)
       Print *, 'Finished eigensolve.'
      EndIf
!-------------------------------------------------------------------

!----------------------------------------------------------------
! All converged solutions are returned. Eigenvalues are stored
! in eigensolve_lambda, vectors in eigensolve_vector.

      call EPSGetConverged(epsolve,n_converged,ierr)
      If (i_proc ==0) Print *, n_converged, ' converged solution(s)'
      open(unit=1,file=trim(path)//file_eigen_freq,status='replace')
      step = 0
      Do iev = 0, n_converged-1
        call EPSGetEigenpair(epsolve,iev,dummy_scalar,PETSC_NULL,       &
     &          dummy_vec,PETSC_NULL,ierr)
        call EPSComputeResidualNorm(epsolve,iev,residual,ierr)
        ev_omega = -aimag(dummy_scalar)
        ev_gamma = real(dummy_scalar)
        If (i_proc == 0) Then
          print 20, iev+1, ev_omega, ev_gamma, residual
          write (1,30) ev_omega, ev_gamma, residual
          call write_prec(10,abs(ev_omega))
          step = step + 1
          call write_prec(10,abs(ev_gamma))
        EndIf
      EndDo
      close(unit=1)
!---------------------------------------------------------------------

      If (i_proc == 0) Then
        call EPSGetIterationNumber(epsolve,n_iter,ierr)
        Print *, n_iter, ' iterations performed'
      EndIf
      call VecDestroy(dummy_vec,ierr)
      call MatDestroy(time_derivative_matrix,ierr)
      call EPSDestroy(epsolve,ierr)

      call SlepcFinalize(ierr)

 
1000  end subroutine

