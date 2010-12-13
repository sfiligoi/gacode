subroutine tgyro_qfm_map

  use tgyro_globals
  use qfm_interface

  implicit none

  qfm_in(:)  = 0.0

  qfm_in(1) = q(i_r)
  qfm_in(2) = s(i_r)
  qfm_in(3) = delta(i_r)
  qfm_in(4) = s_delta(i_r)
  qfm_in(5) = kappa(i_r)
  qfm_in(6) = s_kappa(i_r)
  qfm_in(7) = shift(i_r)
  qfm_in(8) = r_maj(i_r)/r_min
  qfm_in(9) = r(i_r)/r_min
  qfm_in(10) = r_min*dlnnedr(i_r)
  qfm_in(11) = r_min*dlnnidr(1,i_r)
  qfm_in(12) = r_min*dlntedr(i_r)
  qfm_in(13) = r_min*dlntidr(1,i_r)
  qfm_in(14) = ti(1,i_r)/te(i_r) 
  qfm_in(15) = betae_unit(i_r)*loc_betae_scale
  qfm_in(16) = sqrt(mi(1)/(me*loc_me_multiplier))
  qfm_in(17) = nue(i_r)*r_min/c_s(i_r)*loc_nu_scale
  qfm_in(18) = gamma_eb(i_r)*r_min/c_s(i_r)

  ! Default impurities (to prevent div-by-zero)
  qfm_in(19) = 1.0
  qfm_in(21) = 1.0
  qfm_in(23) = 1.0
  qfm_in(25) = 1.0

  ! Proper settings:
  if (loc_n_ion >= 2) then
     ! Ion #2
     qfm_in(19) = zi_vec(2)
     qfm_in(20) = ni(2,i_r)/ne(i_r)
     qfm_in(21) = sqrt(mi(1)/mi(2))
     qfm_in(22) = r_min*dlnnidr(2,i_r)
  endif
  if (loc_n_ion >= 3) then
     ! Ion #3
     qfm_in(23) = zi_vec(3)
     qfm_in(24) = ni(3,i_r)/ne(i_r)
     qfm_in(25) = sqrt(mi(1)/mi(3))
     qfm_in(26) = r_min*dlnnidr(3,i_r)
  endif
  qfm_in(27) = ni(1,i_r)/ne(i_r)

end subroutine tgyro_qfm_map
