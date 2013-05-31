subroutine tgyro_iteration_pppl

  use tgyro_globals
  use tgyro_iteration_variables

  select case (tgyro_iteration_method)

  case (2)  ! Levenberg-Marquardt

     allocate(jfg_transpose(p_max,p_max))
     nu = relax(1)

  case (3)  ! Backtracking

     lambda = 1.0

  end select

  allocate(bfull(p_max))
  allocate(x_best(p_max))
  allocate(f_best(p_max))
  allocate(res_best(p_max))

  !---------------------------------------------------
  ! ZEROTH ITERATION
  !
  ! One pass to get fluxes.
  !
  ! We *do* want the iteration number to increase in 
  ! the case where we restart without iterating (this
  ! can be used to compare transport models).
  !  
  if (tgyro_relax_iterations == 0 .and. loc_restart_flag == 1) i_tran=i_tran+1
  !
  call tgyro_target_vector(x_vec,g_vec)
  if (loc_restart_flag == 0 .or. tgyro_relax_iterations == 0) then
     ! Need to determine initial fluxes
     gyro_restart_method = 1
     call tgyro_flux_vector(x_vec,f_vec,0.0,0)
     gyro_restart_method = 2
  else
     ! Initial fluxes already computed
     p = 0
     do i=2,n_r
        if (loc_ti_feedback_flag == 1) then
           p = p+1
           f_vec(p) = eflux_i_tot(i)
        endif
        if (loc_te_feedback_flag == 1) then
           p = p+1
           f_vec(p) = eflux_e_tot(i)
        endif
        if (loc_ne_feedback_flag == 1) then
           p = p+1
           f_vec(p) = pflux_e_tot(i)
        endif
        if (loc_er_feedback_flag == 1) then
           p = p+1
           f_vec(p) = mflux_tot(i)
        endif
     enddo
     ! GYRO restart data available
     gyro_restart_method = 2
  endif
  res0 = 0.0
  call tgyro_residual(f_vec,g_vec,res,p_max,loc_residual_method)

  if (loc_restart_flag == 0 .or. tgyro_relax_iterations == 0) then
     call tgyro_write_data(1)
  endif
  !----------------------------------------------------

  do i_tran_loop=1,tgyro_relax_iterations

     i_tran = i_tran+1

     ! Initialize gradients

     x_vec0 = x_vec
     g_vec0 = g_vec
     f_vec0 = f_vec
     res0   = res

     !----------------------------------------------
     ! Build dQ^T/dz (dense matrix)
     !
     do p=1,p_max

        x_vec(:) = x_vec0(:)
        x_vec(p) = x_vec0(p)+dx

        call tgyro_target_vector(x_vec,g_vec)
        jg(:,p) = (g_vec(:)-g_vec0(:))/dx

     enddo

     ! Reset gradients
     x_vec = x_vec0

     ! Reset profiles to be consistent with gradient.
     ! NOTE: This corrects an error (pre July 2010).

     p = 0
     do i=2,n_r
        if (loc_ti_feedback_flag == 1) then
           p = p + 1
           dlntidr(therm_vec(:),i) = x_vec(p)
        endif
        if (loc_te_feedback_flag == 1) then
           p = p + 1
           dlntedr(i) = x_vec(p)
        endif
        if (loc_ne_feedback_flag == 1) then
           p = p + 1
           dlnnedr(i) = x_vec(p)
           ! Set dlnnidr(1,i) according to quasineutrality
           call tgyro_quasigrad(ne(i),dlnnedr(i),ni(:,i),dlnnidr(:,i),zi_vec(:),loc_n_ion)
        endif
        if (loc_er_feedback_flag == 1) then
           p = p + 1
           f_rot(i) = x_vec(p)
        endif
     enddo

     call tgyro_profile_functions 
     !----------------------------------------------

     if (tgyro_global_newton_flag == 1) then

        !----------------------------------------------
        ! Build dQ/dz (dense matrix)
        !
        do p=1,p_max

           x_vec(:) = x_vec0(:)
           x_vec(p) = x_vec0(p)+dx

           call tgyro_flux_vector_dense(x_vec,f_vec)
           jf(:,p) = (f_vec(:)-f_vec0(:))/dx

        enddo
        x_vec = x_vec0
        !----------------------------------------------

     else

        !----------------------------------------------
        ! Build dQ/dz (block diagonal matrix)
        !
        ! (p  ,p) (p  ,p+1) (p  ,p+2)
        ! (p+1,p) (p+1,p+1) (p+1,p+2)
        ! (p+2,p) (p+2,p+1) (p+2,p+2)

        jf(:,:) = 0.0

        ip = -1

        if (loc_ti_feedback_flag == 1) then

           ip = ip+1
           call tgyro_flux_vector(x_vec,f_vec,dx,1)
           do p=1,p_max,n_evolve
              do pp=0,n_evolve-1
                 jf(p+pp,p+ip) = (f_vec(p+pp)-f_vec0(p+pp))/dx
              enddo
           enddo

        endif

        if (loc_te_feedback_flag == 1) then

           ip = ip+1
           call tgyro_flux_vector(x_vec,f_vec,dx,2)
           do p=1,p_max,n_evolve
              do pp=0,n_evolve-1
                 jf(p+pp,p+ip) = (f_vec(p+pp)-f_vec0(p+pp))/dx
              enddo
           enddo

        endif

        if (loc_ne_feedback_flag == 1) then

           ip = ip+1
           call tgyro_flux_vector(x_vec,f_vec,dx,3)
           do p=1,p_max,n_evolve
              do pp=0,n_evolve-1
                 jf(p+pp,p+ip) = (f_vec(p+pp)-f_vec0(p+pp))/dx
              enddo
           enddo

        endif

        if (loc_er_feedback_flag == 1) then

           ip = ip+1
           call tgyro_flux_vector(x_vec,f_vec,dx,4)
           do p=1,p_max,n_evolve
              do pp=0,n_evolve-1
                 jf(p+pp,p+ip) = (f_vec(p+pp)-f_vec0(p+pp))/dx
              enddo
           enddo

        endif
        !
        !----------------------------------------------

     endif

     !----------------------------------------------
     ! Total Jacobian: (dQ/dz-dQ^T/dz)
     !
     jfg(:,:) = jf(:,:)-jg(:,:)
     !----------------------------------------------

     sum_res0 = sum(res0(1:p_max))

     b(:) = -(f_vec0(:)-g_vec0(:))
     !---------------------

     ! Levenberg-Marquart Algorithm
     if (tgyro_iteration_method == 2) then


        ! J^TxJ is approximately Hessian of Square Residual
        ! J^Txb is gradient descent of Square Residual
        ! LM solves (J^TJ + nu*Diag(J^TJ))dx=J^T(g-f)
        jfg_transpose(:,:) = transpose(jfg(:,:))
        jfg(:,:) = matmul(jfg_transpose, jfg)
        b(:) = matmul(jfg_transpose, b)

        ! Normalize matrices
        ! Store non-normal matrix for scaling
        jf(:,:) = jfg(:,:)
        do p=1,p_max
           do ip=1,p_max
              jfg(ip,p) = jfg(ip,p)/sqrt(jf(p,p)*jf(ip,ip))
           enddo
           b(p) = b(p) / sqrt(jf(p,p))
        end do

        ! Scale diagonal component by damping factor
        do p=1,p_max
           jfg(p,p) = jfg(p,p) + nu*jfg(p,p)
        enddo

     endif
     !-------------------

     ! LAPACK matrix factorization into L/U components
     call DGETRF(p_max,p_max,jfg,p_max,ipiv,ierr) 

     if (ierr > 0) then
        call tgyro_catch_error('L/U decomposition failure: singular matrix.')
     elseif (ierr < 0) then
        call tgyro_catch_error('L/U decomposition failure: illegal value.')
     endif

     ! LAPACK matrix solve (jfg)x=b (pass jfg and b, return x in b).
     call DGETRS('N',p_max,1,jfg,p_max,ipiv,b,p_max,ierr)

     if (ierr < 0) then
        call tgyro_catch_error('Matrix solve failure: illegal value.')
     endif

     ! Rescale step if doing L-M
     if (tgyro_iteration_method == 2) then
        do p=1,p_max
           b(p) = b(p) / sqrt(jf(p,p))
        enddo
     endif

     ! Check to see if step length exceeds maximum 
     ! Determine direction with greatest overshoot
     lambda = 1.0
     b_flag(:) = ' '
     do p=1,p_max
        if (abs(b(p)) > loc_dx_max/r_min) then
           lambda = min(lambda, loc_dx_max/r_min/abs(b(p)))
           b_flag(p) = '*'
        endif
     enddo
     ! Scale each step by same amount
     b(:) = b(:)*lambda
     bfull(:) = b(:) ! full Newton step, limited by step length

     ! Update gradient using Newton-step.
     x_vec(:) = x_vec0(:)+b(:)

     ! Check to see if gradient is negative
     ! Check for convergence
     ! Just flag for now, could be ok
     do p=1,p_max
        if (x_vec(p) < 0.0) then
           !      x_vec(p) = 0.0001
           b_flag(p) = '#'
        endif
        if (abs(b(p)) < newton_tol) then
           b_flag(p) = '?'
        endif
     enddo

     !----------------------------------------------

     call tgyro_target_vector(x_vec,g_vec)
     call tgyro_flux_vector(x_vec,f_vec,0.0,0)
     !
     ! Compute initial residual
     call tgyro_residual(f_vec,g_vec,res,p_max,loc_residual_method)
     call tgyro_write_intermediate(0,res)
     sum_res = sum(res(1:p_max))     

     ! Levenberg-Marquardt Update
     if (tgyro_iteration_method == 2) then

        if (i_proc_global == 0) then
           write (*,*)
           write (*,*) "Iteration:", i_tran
           write (*,*) "   Residual:", sum_res, "L-M damp factor:", nu
        endif

        ! Residual increased, increase damping factor by lm_boost if not backtracking
        ! If backtracking, increase damping only if backtrack fails
        if (sum_res > sum_res0) then

           if (tgyro_backtrack_method == 0) then
              ! Increase damping
              nu = nu * lm_boost
              nu = max(0.01, nu)

              ! Reset Solution
              x_vec(:) = x_vec0(:)
              call tgyro_target_vector(x_vec,g_vec)
              gyro_restart_method = 1
              call tgyro_flux_vector(x_vec,f_vec,0.0,0)
              gyro_restart_method = 2

              ! Recompute residual
              call tgyro_residual(f_vec,g_vec,res,p_max,loc_residual_method)
              call tgyro_write_intermediate(1,res)
           endif

        else
           ! Accept step, decrease damping factor by multiplying by lm_drop
           nu = nu * lm_drop
           nu = min(0.01, nu)
        endif

        relax(1) = nu

     endif

     ! Backtracking Method
     if (tgyro_backtrack_method > 0) then

        if (i_proc_global == 0) then
           write (*,*) "----------"
           write (*,*) "Iteration:", i_tran
           write (*,*) "Sum_Res0:", sum_res0
        endif

        correct_num = 0 ! number of correction steps
        lambda = 1.0    ! step size length

        ! Store best step so far
        if (sum_res <= sum_res0) then
           sum_res_best = sum_res
           lambda_best = lambda
           x_best = x_vec
           f_best = f_vec
           res_best = res
        else
           sum_res_best = sum_res0
           lambda_best = 0.0
           x_best = x_vec0
           f_best = f_vec0
           res_best = res0
        endif

        do ! Backtracking loop
           lambda_old = lambda

           if(i_proc_global == 0) then
              write (*,*) "Correct Num:", correct_num
              write (*,*) "Lambda:", lambda
              write (*,*) "Sum_Res:", sum_res
              write (*,*)
           endif

           if (sum_res < newton_tol) then
              call tgyro_write_data(1)
              ! This exit call is broken.
              call tgyro_catch_error('Local or global minimum found.')
              exit ! (we have converged either to the real root or a local minimum)
           endif

           ! Golden Ratio Section Search for Minimum of lambda
           if (tgyro_backtrack_method == 1) then
              if (correct_num == 0) then ! initialize

                 lambda_a = 0.0 ! Left and right boundary points
                 lambda_b = 1.0

                 if (sum_res > sum_res0*(1.0-0.29*lambda_old) .and. loc_relax > 1.0) then
                    lambda_c = sum_res0/(sum_res + sum_res0)
                    !                 lambda_c = (1.0 - golden_ratio)**2

                    ! First guess of center value, c
                    lambda = lambda_c
                    correct_flag = 1
                 else ! Full Newton Step was ok, so we don't monkey with line search
                    correct_flag = 0
                 endif

              else

                 ! Pick new optimal center point
                 ! If ends of interval (a,b) converge, we are there and exit loop
                 if ((abs(lambda_a - lambda_b) < 1.0e-2) .or. (sum_res < sum_res0*(1.0-1.0e-4*lambda_old))) then
                    correct_flag = 0
                 else
                    lambda = lambda_c + (1.0-golden_ratio)**2*(lambda_b - lambda_c)
                    correct_flag = 1
                 endif

              endif
           endif

           if (tgyro_backtrack_method == 2) then
              ! Scale back Newton Step if residual has not decreased enough and if we
              ! expect a further error reduction of at least a factor of two
              if (sum_res > sum_res0*(1-0.29*lambda_old) .and. loc_relax > 1.0) then
                 lambda = lambda_old**2 *sum_res0/(sum_res - sum_res0*(1.0-2.0*lambda_old))
                 lambda = max(0.1*lambda_old, lambda)
                 b_flag(:) = '>'
                 correct_flag = 1
              endif
           endif

           ! If the backtracking step gets too small, we may be stuck in a valley and would want
           ! to actually move off of the Newton Direction. To do so, we scale only those directions
           ! that would have their residual increase by the previous step.
           ! If lambda is not too small, scale everything by the same amount.
           ! It would be smarter to move more in a steepest descent direction.
           if (tgyro_backtrack_method == 3) then
              if(lambda < 1.0e-1) then
                 do p=1,p_max
                    if (res(p) > res0(p)) then
                       b_flag(p) = ':'
                       b(p) = bfull(p)  * lambda_old / loc_relax
                       correct_flag = 1
                    else
                       b(p) = bfull(p) * lambda_old
                    endif
                 enddo
                 lambda = lambda_old / loc_relax
              else
                 b(:) = bfull(:) * lambda
              endif
           endif

           b(:) = bfull(:) * lambda
           ! If we've corrected the Newton step, recalculate solution
           ! otherwise, leave correction loop
           ! if we've tried too many corrections, move to next iteration
           ! If we've rotated the propagation direction, move to next iteration
           if ((correct_flag == 0) .or. (correct_num > correct_max)) then
              x_vec = x_best
              f_vec = f_best
              res = res_best
              call tgyro_target_vector(x_vec,g_vec)
              call tgyro_write_intermediate(1,res)
              sum_res = sum(res(:))

              ! If doing L-M and line search failed, increase damping factor
              if ((tgyro_iteration_method == 2) .and. (correct_num > correct_max) .or. (sum_res >= sum_res0)) then
                nu = nu * lm_boost
                nu = max(0.01, nu)
              endif

              ! Diagnostics
              relax(:) = (x_vec(:)-x_vec0(:))/bfull(:)
              relax(1) = nu

              exit
           else
              relax(:) = b(:)/bfull(:) ! For diagnostic purposes
              x_vec(:) = x_vec0(:) + b(:)
              call tgyro_target_vector(x_vec,g_vec)
              call tgyro_flux_vector(x_vec,f_vec,0.0,0)
              call tgyro_residual(f_vec,g_vec,res,p_max,loc_residual_method)
              call tgyro_write_intermediate(1,res)
              sum_res = sum(res(:))

              if (correct_flag == 2) then ! We've rotated off Newton direction
                 correct_flag = 0
                 exit
              endif

              if ((tgyro_backtrack_method == 1)) then
                 if (correct_num == 0) then ! we've taken step size lambda_c
                    sum_res_c = sum_res
                 else

                    ! Change intervals
                    if (sum_res < sum_res_c) then ! min between c and b
                       lambda_a = lambda_c
                       lambda_c = lambda
                       sum_res_c = sum_res
                    else ! min between a and lambda
                       lambda_b = lambda_a
                       lambda_a = lambda
                    endif

                 endif

              endif
              correct_flag = 0
              correct_num = correct_num + 1

              ! Keep track of best residual
              if (sum_res < sum_res_best) then
                 lambda_best = lambda
                 res_best(:) = res(:)
                 x_best(:) = x_vec(:)
                 f_best(:) = f_vec(:)
              endif
           endif

        end do ! step correction loop

     endif

     call tgyro_write_data(1)

  enddo

end subroutine tgyro_iteration_pppl
