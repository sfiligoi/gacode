program test

  use geo

  double precision, dimension(1) :: theta

  theta(1) = 0.1
  
  call geo_interp(1,theta,.true.)
  
  print *,geo_b
  
end program test
