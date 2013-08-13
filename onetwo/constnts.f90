  MODULE CONSTNTS 
      IMPLICIT NONE
!
! --- store all conversion factors, mathematical constants,
! --- and physical constants here
!
      real *8                                                   &
             pi, pisq, twopi, piov4, piov2, tpiov2, tpiov4,     &
             fpiov4, spiov4,u0, fourpioc, fourpi, psimks,       &
             psikgaus,rootpi,root2,charg4,cee,xmasse,xmassp,    &
             charge,rot2pi, rotpi2,joupkev, jouperg,kevperg
!
! --- initialize constants in /constnts/ common block (INCLUDE file constnts.i)
!
      data                                                            &
           pi, pisq, twopi, piov4, piov2, tpiov2, tpiov4, fpiov4,     &
           spiov4, fourpioc, fourpi, psimks, psikgaus,            &
           rootpi,root2,rot2pi, rotpi2                                &
         / 3.14159265358978, 9.86960440108935, 6.28318530717958,      &
           0.785398163397448, 1.57079632679489, 4.712388904,          &
           2.356194491, 3.926990818, 5.497787145,                     &
           4.1916900439e-10, 12.5663706144, 1.0e-5, 1.0e5 ,           &
           1.772453851,1.414213562,2.506628275, 1.253314137 /
!
!
      data  u0    /1.2566370614e-6/     ! H/m
      data xmasse / 9.10938970e-28 /    ! electron mass   (g)
      data xmassp / 1.67262310e-24 /    ! proton   mass   (g)
      data charge / 4.80320670e-10 /    ! electron charge (esu)
      data charg4 / 5.32261610e-38 /    ! charge**4       (esu**4)
      data cee    / 2.99792458e+10 /    ! speed of light  (cm/s)
      data       joupkev, jouperg /1.60217733e-16, 1.0e-07/
      data       kevperg /6.2415064e+08/
!
!
!
  END   MODULE CONSTNTS 
