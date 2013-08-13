  MODULE rad_loss

  USE nrtype,             ONLY :I4B,DP

  SAVE
  REAL(DP),SAVE,ALLOCATABLE,DIMENSION(:,:) :: brems_nions
  REAL(DP),SAVE ::  brems_prim,brems_imp
  !brems_nions is dynamically sized to (nj,nprim+nimp=nions)




  END MODULE rad_loss
