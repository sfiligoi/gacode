! This is called after cgyro_kernel to perform clean-ups from CPU and GPU (deallocate)

#define cgyro_del(x) \
  if (allocated(x)) deallocate(x)

#if defined(OMPGPU)

#define cgyro_del_gpu(x) \
  if (allocated(xl)) then \
!$omp target exit data map(release:x) \
     deallocate(x) \
  endif

#define cgyro_del_fixed_gpu(x) \
!$omp target exit data map(release:x)

#elif defined(_OPENACC)

#define cgyro_del_gpu(x) \
  if (allocated(xl)) then \
!$acc exit data delete(x) \
     deallocate(x) \
  endif

#define cgyro_del_fixed_gpu(x) \
!$acc exit data delete(x)

#else

#define cgyro_del_gpu(x) \
  if (allocated(x)) deallocate(x)

  ! nothing to do
#define cgyro_del_fixed_gpu(x)

#endif

subroutine cgyro_cleanup
  use cgyro_globals
  use parallel_lib

#if defined(_OPENACC) || defined(OMPGPU)
#define CGYRO_GPU_FFT
#endif

  !!!!!!!!!!!!!!!!!!!!!!!!!!
  ! From cgyro_init_manager
  !!!!!!!!!!!!!!!!!!!!!!!!!!
  
  cgyro_del(energy)
  cgyro_de_gpul(vel)
  cgyro_del(w_e)
  cgyro_del(e_deriv1_mat)
  cgyro_del(e_deriv1_rot_mat)
  cgyro_del_gpu(xi)
  cgyro_del(w_xi)
  cgyro_del(xi_lor_mat)
  cgyro_del(xi_deriv_mat)
  
  cgyro_del(theta)
  cgyro_del(thetab)
  cgyro_del(w_theta)
  cgyro_del(g_theta)
  cgyro_del(g_theta_geo)
  cgyro_del(bmag)
  cgyro_del(btor)
  cgyro_del(bpol)
  cgyro_del(k_perp)
  cgyro_del(k_x)
  cgyro_del(bigr)
  cgyro_del(bigr_r)
  cgyro_del(itp)
  cgyro_del_gpu(omega_stream)
  cgyro_del(omega_trap)
  cgyro_del(omega_rdrift)
  cgyro_del(omega_adrift)
  cgyro_del(omega_aprdrift)
  cgyro_del(omega_cdrift)
  cgyro_del(omega_cdrift_r)
  cgyro_del(omega_gammap)
  cgyro_del(lambda_rot)
  cgyro_del(dlambda_rot)
  cgyro_del(dens_rot)
  cgyro_del(dens_ele_rot)
  cgyro_del(dens_avg_rot)
  cgyro_del(dlnndr_avg_rot)
  cgyro_del(omega_rot_trap)
  cgyro_del(omega_rot_u)
  cgyro_del(omega_rot_drift)
  cgyro_del(omega_rot_drift_r)
  cgyro_del(omega_rot_edrift)
  cgyro_del(omega_rot_edrift_r)
  cgyro_del(omega_rot_star)
  cgyro_del(gtime)
  cgyro_del(freq)
  cgyro_del(freq_err)
  cgyro_del_gpu(fcoef)
  cgyro_del_gpu(gcoef)
  cgyro_del_gpu(field)
  cgyro_del_gpu(field_loc)
  cgyro_del(field_dot)
  cgyro_del(field_old)
  cgyro_del(field_old2)
  cgyro_del(field_old3)
  cgyro_del(moment)
  cgyro_del(moment_loc)
  cgyro_del(cflux)
  cgyro_del(cflux_loc)
  cgyro_del(gflux)
  cgyro_del(gflux_loc)
  cgyro_del(cflux_tave)
  cgyro_del(gflux_tave)
  cgyro_del(recv_status)
  cgyro_del_gpu(icd_c)
  cgyro_del_gpu(dtheta)
  cgyro_del_gpu(dtheta_up)
  cgyro_del_gpu(source)
  cgyro_del_gpu(h0_old)
  cgyro_del_gpu(rhs)
  cgyro_del_gpu(h_x)
  cgyro_del_gpu(g_x)
  cgyro_del_gpu(h0_x)
  cgyro_del_gpu(cap_h_c)
  cgyro_del_gpu(cap_h_ct)
  cgyro_del_gpu(cap_h_c_dot)
  cgyro_del_gpu(cap_h_c_old)
  cgyro_del_gpu(cap_h_c_old2)
  cgyro_del_gpu(omega_cap_h)
  cgyro_del_gpu(omega_h)
  cgyro_del_gpu(omega_s)
  cgyro_del_gpu(omega_ss)
  cgyro_del_gpu(jvec_c)
  cgyro_del_gpu(jvec_c_nl)
  cgyro_del(jvec_v)
  cgyro_del_gpu(dvjvec_c)
  cgyro_del_gpu(dvjvec_v)
  cgyro_del_gpu(jxvec_c)
  cgyro_del_gpu(upfac1)
  cgyro_del_gpu(upfac2)
  cgyro_del_gpu(cap_h_v)
  cgyro_del_gpu(upwind_res_loc)
  cgyro_del_gpu(upwind_res)
  cgyro_del_gpu(upwind32_res_loc)
  cgyro_del_gpu(upwind32_res)
  cgyro_del_gpu(f_nl)
  cgyro_del_gpu(g_nl)
  cgyro_del_gpu(fpack)
  cgyro_del_gpu(gpack)
  cgyro_del_gpu(cmat)
  cgyro_del_gpu(cmat_fp32)
  cgyro_del_gpu(cmat_stripes)
  cgyro_del_gpu(cmat_simple)

#ifndef CGYRO_GPU_FFT
  cgyro_del(fx)
  cgyro_del(gx)
  cgyro_del(fy)
  cgyro_del(gy)
  cgyro_del(uxmany)
  cgyro_del(uymany)
  cgyro_del(vx)
  cgyro_del(vy)
  cgyro_del(uv)
#endif  

#ifdef CGYRO_GPU_FFT
  cgyro_del_gpu(fxmany)
  cgyro_del_gpu(gxmany)
  cgyro_del_gpu(fymany)
  cgyro_del_gpu(gymany)
  cgyro_del_gpu(uxmany)
  cgyro_del_gpu(uymany)
  cgyro_del_gpu(vxmany)
  cgyro_del_gpu(vymany)
  cgyro_del_gpu(uvmany)
#endif    

  !!!!!!!!!!!!!!!!!!!!!!!!!!
  ! From cgyro_make_profiles
  !!!!!!!!!!!!!!!!!!!!!!!!!!
  
  cgyro_del_gpu(px)
  cgyro_del(geo_yin)
 
  cgyro_del_fixed_gpu(z)
  cgyro_del_fixed_gpu(temp) 
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!
  ! From cgyro_init_arrays
  !!!!!!!!!!!!!!!!!!!!!!!!!!

  cgyro_del(vfac)
  cgyro_del(sum_den_h)
  cgyro_del(cderiv)
  cgyro_del(uderiv)
  cgyro_del_gpu(c_wave)
  cgyro_del(hzf)
  cgyro_del(xzf)

  !!!!!!!!!!!!!!!!!!!!!!!!!!
  ! From cgyro_mpi_grid
  !!!!!!!!!!!!!!!!!!!!!!!!!!

  cgyro_del_gpu(ie_v)
  cgyro_del_gpu(ix_v)
  cgyro_del_gpu(is_v)
  cgyro_del_gpu(iv_v)
  cgyro_del_gpu(ir_c)
  cgyro_del_gpu(it_c)
  cgyro_del_gpu(ic_c)
  cgyro_del(ica_c)
  cgyro_del(icb_c)

  !!!!!!!!!!!!!!!!!!!!!!!!!!
  ! From cgyro_parallel_lib
  !!!!!!!!!!!!!!!!!!!!!!!!!!

  call parallel_lib_clean
  
end subroutine cgyro_cleanup
 
