       module aid_newton
         
            character *8,save, dimension(:),allocatable :: nameu_gp
            real *8 , save, dimension(:),allocatable :: eqtn_scale
            integer *4 ,save,dimension(:),allocatable ::grid_point
            integer *4 ,save,dimension(:),allocatable ::itranrd
            integer *4 freeze_xte,freeze_xti,freeze_xrbp,freeze_xwr
            integer *4 freeze_xni, freeze_xnue,freeze_alpha_exp
            integer *4 freeze_source,freeze_xsn,tot_iters,jac_skip
            integer *4 jelc_clamp, random_pert,conv_skip,jion_clamp
            integer *4 ,save ::  wrt_nwt,tipts,tipts2,itncount
            real *8 ,save,dimension(:),allocatable ::typx
            real *8    SSQR,gradmax,SSQR_HIST
            logical xvset,jac_known
            data xvset,jac_known /.false., .false./ 
            data tipts /0/      !total iters per time step counter,used
                                !when wrt_nwt = 1
            data tipts2 /0/     !counts number of time steps allowed in
                                !file nwt.test before file is rewritten
                                !tipts is used only if wrt_nwt > 1
       end module aid_newton
