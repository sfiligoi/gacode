    MODULE limiter
     integer,parameter :: maxlimpt = 200
!
!  note that the extremes of the limiter are stored as
!    xlimiter(nlimiter+1) = xlimmin
!                      2    xlimmax
!    ylimiter(nlimiter+1) = ylimmin
!                      2    ylimmax
!  so nlimiter must be .le. maxlimpt-2
!
     real*8, dimension(:) :: xlimiter(maxlimpt), &
                             ylimiter(maxlimpt)
     integer *4 nlimiter
    END MODULE limiter
