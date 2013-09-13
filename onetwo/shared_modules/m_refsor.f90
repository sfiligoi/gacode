MODULE m_refsor
  USE nrtype,                   ONLY : DP,SP,I4B
  !INTEGER, PARAMETER :: kdp = SELECTED_REAL_KIND(15)
   INTEGER, PARAMETER :: kdp = DP
  PUBLIC :: refsor
  PRIVATE :: kdp
  !PRIVATE :: R_refsor, I_refsor, D_refsor
  !PRIVATE :: R_inssor, I_inssor, D_inssor
  !PRIVATE :: R_subsor, I_subsor, D_subsor
  INTERFACE refsor
     MODULE PROCEDURE D_refsor, R_refsor, I_refsor
  END INTERFACE
CONTAINS

  SUBROUTINE D_refsor (XDONT)
    !  Sorts XDONT into ascending order - Quicksort
    ! __________________________________________________________
    !  Quicksort chooses a "pivot" in the set, and explores the
    !  array from both ends, looking for a value > pivot with the
    !  increasing index, for a value <= pivot with the decreasing
    !  index, and swapping them when it has found one of each.
    !  The array is then subdivided in 2 ([3]) subsets:
    !  { values <= pivot} {pivot} {values > pivot}
    !  One then call recursively the program to sort each subset.
    !  When the size of the subarray is small enough, one uses an
    !  insertion sort that is faster for very small sets.
    !  Michel Olagnon - Apr. 2000
    ! __________________________________________________________
    ! __________________________________________________________
!    REAL (kind=kdp), DIMENSION (:), INTENT (InOut) :: XDONT
     REAL(DP), DIMENSION (:), INTENT (InOut) :: XDONT
    ! __________________________________________________________
    !
    !
    CALL D_subsor (XDONT, 1, SIZE (XDONT))
    CALL D_inssor (XDONT)
    RETURN
  END SUBROUTINE D_refsor
  RECURSIVE SUBROUTINE D_subsor (XDONT, IDEB1, IFIN1)
    !  Sorts XDONT from IDEB1 to IFIN1
    ! __________________________________________________________
    REAL(kind=kdp), DIMENSION (:), INTENT (InOut) :: XDONT
    INTEGER, INTENT (In) :: IDEB1, IFIN1
    ! __________________________________________________________
    INTEGER, PARAMETER :: NINS = 16 ! Max for insertion sort
    INTEGER :: ICRS, IDEB, IDCR, IFIN, IMIL
    REAL(kind=kdp) :: XPIV, XWRK
    !
    IDEB = IDEB1
    IFIN = IFIN1
    !
    !  If we don't have enough values to make it worth while, we leave
    !  them unsorted, and the final insertion sort will take care of them
    !
    IF ((IFIN - IDEB) > NINS) THEN
       IMIL = (IDEB+IFIN) / 2
       !
       !  One chooses a pivot, median of 1st, last, and middle values
       !
       IF (XDONT(IMIL) < XDONT(IDEB)) THEN
          XWRK = XDONT (IDEB)
          XDONT (IDEB) = XDONT (IMIL)
          XDONT (IMIL) = XWRK
       END IF
       IF (XDONT(IMIL) > XDONT(IFIN)) THEN
          XWRK = XDONT (IFIN)
          XDONT (IFIN) = XDONT (IMIL)
          XDONT (IMIL) = XWRK
          IF (XDONT(IMIL) < XDONT(IDEB)) THEN
             XWRK = XDONT (IDEB)
             XDONT (IDEB) = XDONT (IMIL)
             XDONT (IMIL) = XWRK
          END IF
       END IF
       XPIV = XDONT (IMIL)
       !
       !  One exchanges values to put those > pivot in the end and
       !  those <= pivot at the beginning
       !
       ICRS = IDEB
       IDCR = IFIN
       ECH2: DO
          DO
             ICRS = ICRS + 1
             IF (ICRS >= IDCR) THEN
                !
                !  the first  >  pivot is IDCR
                !  the last   <= pivot is ICRS-1
                !  Note: If one arrives here on the first iteration, then
                !        the pivot is the maximum of the set, the last value is equal
                !        to it, and one can reduce by one the size of the set to process,
                !        as if XDONT (IFIN) > XPIV
                !
                EXIT ECH2
                !
             END IF
             IF (XDONT(ICRS) > XPIV) EXIT
          END DO
          DO
             IF (XDONT(IDCR) <= XPIV) EXIT
             IDCR = IDCR - 1
             IF (ICRS >= IDCR) THEN
                !
                !  The last value < pivot is always ICRS-1
                !
                EXIT ECH2
             END IF
          END DO
          !
          XWRK = XDONT (IDCR)
          XDONT (IDCR) = XDONT (ICRS)
          XDONT (ICRS) = XWRK
       END DO ECH2
       !
       !  One now sorts each of the two sub-intervals
       !
       CALL D_subsor (XDONT, IDEB1, ICRS-1)
       CALL D_subsor (XDONT, IDCR, IFIN1)
    END IF
    RETURN
  END SUBROUTINE D_subsor
  SUBROUTINE D_inssor (XDONT)
    !  Sorts XDONT into increasing order (Insertion sort)
    ! __________________________________________________________
    REAL(kind=kdp), DIMENSION (:), INTENT (InOut) :: XDONT
    ! __________________________________________________________
    INTEGER :: ICRS, IDCR
    REAL(kind=kdp) :: XWRK
    !
    DO ICRS = 2, SIZE (XDONT)
       XWRK = XDONT (ICRS)
       IF (XWRK >= XDONT(ICRS-1)) CYCLE
       XDONT (ICRS) = XDONT (ICRS-1)
       DO IDCR = ICRS - 2, 1, - 1
          IF (XWRK >= XDONT(IDCR)) EXIT
          XDONT (IDCR+1) = XDONT (IDCR)
       END DO
       XDONT (IDCR+1) = XWRK
    END DO
    !
    RETURN
    !
  END SUBROUTINE D_inssor
  !
  SUBROUTINE R_refsor (RXDONT)
    !  Sorts XDONT into ascending order - Quicksort
    ! __________________________________________________________
    !  Quicksort chooses a "pivot" in the set, and explores the
    !  array from both ends, looking for a value > pivot with the
    !  increasing index, for a value <= pivot with the decreasing
    !  index, and swapping them when it has found one of each.
    !  The array is then subdivided in 2 ([3]) subsets:
    !  { values <= pivot} {pivot} {values > pivot}
    !  One then call recursively the program to sort each subset.
    !  When the size of the subarray is small enough, one uses an
    !  insertion sort that is faster for very small sets.
    !  Michel Olagnon - Apr. 2000
    ! __________________________________________________________
    ! _________________________________________________________
    REAL(SP), DIMENSION (:), INTENT (InOut) :: RXDONT
    ! __________________________________________________________
    !
    !
    CALL R_subsor (RXDONT, 1, SIZE (RXDONT))
    CALL R_inssor (RXDONT)
    RETURN
  END SUBROUTINE R_refsor
  RECURSIVE SUBROUTINE R_subsor (XDONT, IDEB1, IFIN1)
    !  Sorts XDONT from IDEB1 to IFIN1
    ! __________________________________________________________
    REAL, DIMENSION (:), INTENT (InOut) :: XDONT
    INTEGER, INTENT (In) :: IDEB1, IFIN1
    ! __________________________________________________________
    INTEGER, PARAMETER :: NINS = 16 ! Max for insertion sort
    INTEGER :: ICRS, IDEB, IDCR, IFIN, IMIL
    REAL :: XPIV, XWRK
    !
    IDEB = IDEB1
    IFIN = IFIN1
    !
    !  If we don't have enough values to make it worth while, we leave
    !  them unsorted, and the final insertion sort will take care of them
    !
    IF ((IFIN - IDEB) > NINS) THEN
       IMIL = (IDEB+IFIN) / 2
       !
       !  One chooses a pivot, median of 1st, last, and middle values
       !
       IF (XDONT(IMIL) < XDONT(IDEB)) THEN
          XWRK = XDONT (IDEB)
          XDONT (IDEB) = XDONT (IMIL)
          XDONT (IMIL) = XWRK
       END IF
       IF (XDONT(IMIL) > XDONT(IFIN)) THEN
          XWRK = XDONT (IFIN)
          XDONT (IFIN) = XDONT (IMIL)
          XDONT (IMIL) = XWRK
          IF (XDONT(IMIL) < XDONT(IDEB)) THEN
             XWRK = XDONT (IDEB)
             XDONT (IDEB) = XDONT (IMIL)
             XDONT (IMIL) = XWRK
          END IF
       END IF
       XPIV = XDONT (IMIL)
       !
       !  One exchanges values to put those > pivot in the end and
       !  those <= pivot at the beginning
       !
       ICRS = IDEB
       IDCR = IFIN
       ECH2: DO
          DO
             ICRS = ICRS + 1
             IF (ICRS >= IDCR) THEN
                !
                !  the first  >  pivot is IDCR
                !  the last   <= pivot is ICRS-1
                !  Note: If one arrives here on the first iteration, then
                !        the pivot is the maximum of the set, the last value is equal
                !        to it, and one can reduce by one the size of the set to process,
                !        as if XDONT (IFIN) > XPIV
                !
                EXIT ECH2
                !
             END IF
             IF (XDONT(ICRS) > XPIV) EXIT
          END DO
          DO
             IF (XDONT(IDCR) <= XPIV) EXIT
             IDCR = IDCR - 1
             IF (ICRS >= IDCR) THEN
                !
                !  The last value < pivot is always ICRS-1
                !
                EXIT ECH2
             END IF
          END DO
          !
          XWRK = XDONT (IDCR)
          XDONT (IDCR) = XDONT (ICRS)
          XDONT (ICRS) = XWRK
       END DO ECH2
       !
       !  One now sorts each of the two sub-intervals
       !
       CALL R_subsor (XDONT, IDEB1, ICRS-1)
       CALL R_subsor (XDONT, IDCR, IFIN1)
    END IF
    RETURN
  END SUBROUTINE R_subsor
  SUBROUTINE R_inssor (XDONT)
    !  Sorts XDONT into increasing order (Insertion sort)
    ! __________________________________________________________
    REAL, DIMENSION (:), INTENT (InOut) :: XDONT
    ! __________________________________________________________
    INTEGER :: ICRS, IDCR
    REAL :: XWRK
    !
    DO ICRS = 2, SIZE (XDONT)
       XWRK = XDONT (ICRS)
       IF (XWRK >= XDONT(ICRS-1)) CYCLE
       XDONT (ICRS) = XDONT (ICRS-1)
       DO IDCR = ICRS - 2, 1, - 1
          IF (XWRK >= XDONT(IDCR)) EXIT
          XDONT (IDCR+1) = XDONT (IDCR)
       END DO
       XDONT (IDCR+1) = XWRK
    END DO
    !
    RETURN
    !
  END SUBROUTINE R_inssor
  !
  SUBROUTINE I_refsor (IXDONT)
    !  Sorts XDONT into ascending order - Quicksort
    ! __________________________________________________________
    !  Quicksort chooses a "pivot" in the set, and explores the
    !  array from both ends, looking for a value > pivot with the
    !  increasing index, for a value <= pivot with the decreasing
    !  index, and swapping them when it has found one of each.
    !  The array is then subdivided in 2 ([3]) subsets:
    !  { values <= pivot} {pivot} {values > pivot}
    !  One then call recursively the program to sort each subset.
    !  When the size of the subarray is small enough, one uses an
    !  insertion sort that is faster for very small sets.
    !  Michel Olagnon - Apr. 2000
    ! __________________________________________________________
    ! __________________________________________________________
    INTEGER(I4B), DIMENSION (:), INTENT (InOut)  :: IXDONT
    ! __________________________________________________________
    !
    !
    CALL I_subsor (IXDONT, 1, SIZE (IXDONT))
    CALL I_inssor (IXDONT)
    RETURN
  END SUBROUTINE I_refsor
  RECURSIVE SUBROUTINE I_subsor (XDONT, IDEB1, IFIN1)
    !  Sorts XDONT from IDEB1 to IFIN1
    ! __________________________________________________________
    INTEGER, DIMENSION (:), INTENT (InOut) :: XDONT
    INTEGER, INTENT (In) :: IDEB1, IFIN1
    ! __________________________________________________________
    INTEGER, PARAMETER :: NINS = 16 ! Max for insertion sort
    INTEGER :: ICRS, IDEB, IDCR, IFIN, IMIL
    INTEGER :: XPIV, XWRK
    !
    IDEB = IDEB1
    IFIN = IFIN1
    !
    !  If we don't have enough values to make it worth while, we leave
    !  them unsorted, and the final insertion sort will take care of them
    !
    IF ((IFIN - IDEB) > NINS) THEN
       IMIL = (IDEB+IFIN) / 2
       !
       !  One chooses a pivot, median of 1st, last, and middle values
       !
       IF (XDONT(IMIL) < XDONT(IDEB)) THEN
          XWRK = XDONT (IDEB)
          XDONT (IDEB) = XDONT (IMIL)
          XDONT (IMIL) = XWRK
       END IF
       IF (XDONT(IMIL) > XDONT(IFIN)) THEN
          XWRK = XDONT (IFIN)
          XDONT (IFIN) = XDONT (IMIL)
          XDONT (IMIL) = XWRK
          IF (XDONT(IMIL) < XDONT(IDEB)) THEN
             XWRK = XDONT (IDEB)
             XDONT (IDEB) = XDONT (IMIL)
             XDONT (IMIL) = XWRK
          END IF
       END IF
       XPIV = XDONT (IMIL)
       !
       !  One exchanges values to put those > pivot in the end and
       !  those <= pivot at the beginning
       !
       ICRS = IDEB
       IDCR = IFIN
       ECH2: DO
          DO
             ICRS = ICRS + 1
             IF (ICRS >= IDCR) THEN
                !
                !  the first  >  pivot is IDCR
                !  the last   <= pivot is ICRS-1
                !  Note: If one arrives here on the first iteration, then
                !        the pivot is the maximum of the set, the last value is equal
                !        to it, and one can reduce by one the size of the set to process,
                !        as if XDONT (IFIN) > XPIV
                !
                EXIT ECH2
                !
             END IF
             IF (XDONT(ICRS) > XPIV) EXIT
          END DO
          DO
             IF (XDONT(IDCR) <= XPIV) EXIT
             IDCR = IDCR - 1
             IF (ICRS >= IDCR) THEN
                !
                !  The last value < pivot is always ICRS-1
                !
                EXIT ECH2
             END IF
          END DO
          !
          XWRK = XDONT (IDCR)
          XDONT (IDCR) = XDONT (ICRS)
          XDONT (ICRS) = XWRK
       END DO ECH2
       !
       !  One now sorts each of the two sub-intervals
       !
       CALL I_subsor (XDONT, IDEB1, ICRS-1)
       CALL I_subsor (XDONT, IDCR, IFIN1)
    END IF
    RETURN
  END SUBROUTINE I_subsor
  SUBROUTINE I_inssor (XDONT)
    !  Sorts XDONT into increasing order (Insertion sort)
    ! __________________________________________________________
    INTEGER, DIMENSION (:), INTENT (InOut) :: XDONT
    ! __________________________________________________________
    INTEGER :: ICRS, IDCR
    INTEGER :: XWRK
    !
    DO ICRS = 2, SIZE (XDONT)
       XWRK = XDONT (ICRS)
       IF (XWRK >= XDONT(ICRS-1)) CYCLE
       XDONT (ICRS) = XDONT (ICRS-1)
       DO IDCR = ICRS - 2, 1, - 1
          IF (XWRK >= XDONT(IDCR)) EXIT
          XDONT (IDCR+1) = XDONT (IDCR)
       END DO
       XDONT (IDCR+1) = XWRK
    END DO
    !
    RETURN
    !
  END SUBROUTINE I_inssor
  !
END MODULE m_refsor
