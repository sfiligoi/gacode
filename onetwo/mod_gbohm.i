c
c --- INCLUDE file mod_gbohm.i
c     modified bohm, gyrobohm, internal transport
c     barrier confinement models.
c

c
      common /mod_gbohm/
     .                  chie_mgb(kj),chii_mgb(kj),chie_bohm(kj),
     .                  chie_gb(kj),chii_bohm(kj),chii_bitb(kj),
     .                  chie_bgb(kj),chii_bgbitb(kj),chiif_mgb(kj),
     .                  ce0_mgb,alpha_mgb,ci1_mgb,chief_mgb(kj),
     .                  xke_mgb(kj),xki_mgb(kj),chii_gb(kj),
     .                  ci2_mgb,ce_bgb,cfe_mgb,cfe_bgb,
     .                  cfi_mgb,cfi_bgbc1_g,c2_g,c_theta,
     .                  include_itb
