  MODULE shot_info

    USE nrtype,   ONLY : I4B,dp
    IMPLICIT NONE
    TYPE shot_attrib
       INTEGER(I4B) shot_nmbr
       REAL(DP)     shot_time
    END TYPE shot_attrib
    TYPE (shot_attrib) shot_id

  END MODULE shot_info
