      module precision_m
      implicit none
      integer, parameter :: sp = selected_real_kind(6,30)
      integer, parameter :: dp = selected_real_kind(14,60)
      integer, parameter :: wp = dp

      integer, parameter :: singlePrecision = sp
      integer, parameter :: doublePrecision = dp
      integer, parameter :: fp_kind = wp

      integer, parameter :: i4 = selected_int_kind(9)
      integer, parameter :: i8 = selected_int_kind(14)
      end module precision_m

