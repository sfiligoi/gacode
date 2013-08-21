c
c --- INCLUDE file rebut.i
c
      common /rebut/ crit_grad(kj), grad_te(kj), xchierl(kj),
     .               xchiirl(kj), xeffrl(kj), xeffsf(kj), eeneav(kj),
     .               enpav(kj), wrebut, drebut(kj), xkerebut(kj),
     .               xkirebut(kj), timerebut, relaxrebut, xchierlsv(kj),
     .               xchiirlsv(kj), xkangrl(kj), xkangrlsv(kj),
     .               itercorct, qrebsmth, tirlw, rlw_model
      integer        itercorct, qrebsmth, tirlw
      character*3                                rlw_model   ! 22 Aug 95
c
c --- xeffsf is single fluid chi effective, which is determined from
c --- analysis or simulation mode, depending on how the code is run.
c --- xeffrl is single fluid chi effective based on rebut model of
c --- chie and chi. drebut is particle diffusion coefficient.
c --- relaxrebut is (under) relaxation parameter for diffusivity.
c --- itercorct keeps track of nuymber of corrector iterations done
c --- qrebsmth = 0 no effect = 1, smooth the q profile used in
c --- Rebut-Lallia model (q is not smoothed in other parts of the code)
c
