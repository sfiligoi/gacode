 MODULE file_proc

    CHARACTER (len = 256)int_io,new_name


    CONTAINS 

    SUBROUTINE append_time_to_filename (time,filename)
!----------------------------------------------------------
! -- 
!----------------------------------------------------------
    CHARACTER (len=*) filename
    CHARACTER *4 id
    INTEGER j
       !find location of .nc or .cdf  (if netcdf file)
       new_name = ADJUSTL(filename)
       j = 0
       j = INDEX(new_name,'.nc',BACK = .TRUE. )
       IF( j .NE. 0)id ='.nc '
       IF(j == 0)THEN
          j = INDEX(new_name,'.cdf',BACK = .TRUE. )
          id = '.cdf'
       ENDIF
       IF(j == 0)THEN
          j = LEN_TRIM(new_name)
          id =''
       ENDIF

       IF(j == 0) RETURN 
     print *,'old  filename =',filename
       WRITE(int_io,'A,A,1pe12.6,A')new_name(1:j-1),'_',time,//id
       READ (int_i0,'A')new_name
     print *,'new filename =',new_name

    RETURN
    END     SUBROUTINE append_time_to_filename

    SUBROUTINE rename_disk_file(orig_name)
!----------------------------------------------------------------------
! -- rename file in current working directory fromorig_name to new_name
! -- Note that new_name is assumed loaded in this module by a previous
! -- call to append_time_to_filename
!-----------------------------------------------------------------------

     CHARACER(len=*) orig_name,new_name

        call system("mv " // trim(orig_name) // " " // trim(new_name))

     RETURN
     END SUBROUTINE rename_disk_file

 END MODULE file_proc
