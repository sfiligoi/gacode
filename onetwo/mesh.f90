      MODULE MESH
      USE param, only : kj,kjm1,kbctim
      implicit none
!
    
      integer*4,public ::                                            &
           te_var_edge,ti_var_edge,rot_var_edge,ni_var_edge,         & !jmp.den
           te_index,ti_index,rot_index,ni_index,                     & !jmp.den
           te_index_save,ti_index_save,rot_index_save,ni_index_save, & !jmp.den
           imesh
      logical,public  ::                                             &
           freeze_te_index,freeze_ti_index,freeze_rot_index,         & !jmp.den
           freeze_ni_index                                             !jmp.den
      real*8,public ::                                               &
                   r(kj), ra(kjm1), dr(kjm1), drr(kj), rrm(kj),      &
                   rrp(kjm1), roa(kj), r_mesh(kj),                   &
                   grho1_mesh(kj), grho2_mesh(kj), rho_edge,         &
                   r_mesh_mod_edge(kj),r_at_t0(kj),reqdsk_box_edge,  &
                   fix_edge_te(kbctim), fix_edge_ti(kbctim),         &
                   fix_edge_rot(kbctim),fix_edge_ni(kbctim)
!                  NOTE: fix_edge_ni not yet used 8/9/05 HSJ
!
!                   grho1_mesh = <grad rho>,grho2_mesh = < (grad rho )**2 >
!                   are defined on the r grid in rhoset.
 
!  fix_edge_te,ti,rot can be in (0.0,1.0) normalized flux or in grid points
!  (grid points stored  as floating point numbers )
      END MODULE MESH
