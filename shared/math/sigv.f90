real function sigv(ti)

  real, intent(in) :: ti
  real :: r0,a1,a2,a3,a4,a5,a6

  r0  = 0.2935
  a1 = -21.378
  a2 = -25.204
  a3 = -7.101e-2
  a4 = 1.938e-4
  a5 = 4.925e-6
  a6 = -3.984e-8

  sigv = exp(a1/ti**r0+a2+ti*(a3+ti*(a4+ti*(a5+ti*a6))))

end function sigv
