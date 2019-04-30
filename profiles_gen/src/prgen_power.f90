subroutine prgen_power

  use expro

  implicit none
  real, dimension(expro_n_exp) :: temp
  
  ! Total auxiliary electron power  
  temp = expro_qohme+expro_qbeame+expro_qrfe+expro_qione
  call volint(temp,expro_pow_e_aux)
  ! Total electron power 
  temp = temp+expro_qbrem+expro_qsync+expro_qline-expro_qei+expro_qfuse
  call volint(temp,expro_pow_e)

  ! Total auxiliary ion power 
  temp = expro_qbeami+expro_qrfi+expro_qioni+expro_qcxi
  call volint(temp,expro_pow_i_aux)
  ! Total ion power 
  temp = temp+expro_qei+expro_qfusi
  call volint(temp,expro_pow_i)

  ! Exchange power
  call volint(expro_qei,expro_pow_ei)

  ! Fusion power
  call volint(expro_qfuse,expro_pow_e_fus)
  call volint(expro_qfusi,expro_pow_i_fus)

  ! Radiated power (sink/negative)
  call volint(expro_qbrem,expro_pow_e_brem)
  call volint(expro_qsync,expro_pow_e_sync)
  call volint(expro_qline,expro_pow_e_line)

  ! Particle/momentum
  call volint(expro_qpar,expro_flow_beam)
  call volint(expro_qmom,expro_flow_mom)
 
end subroutine prgen_power

subroutine volint(f,fdv)

  use expro
  
  implicit none

  integer :: i
  real, intent(in) :: f(expro_n_exp)
  real, intent(out) :: fdv(expro_n_exp)

  fdv(1) = 0.0

  ! Integration is exact for constant f (density)
  do i=2,expro_n_exp
     fdv(i) = fdv(i-1)+0.5*(f(i)+f(i-1))*(expro_vol(i)-expro_vol(i-1))
  enddo

end subroutine volint
