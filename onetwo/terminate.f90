
	MODULE  terminate
	CONTAINS
		subroutine stop(message,iounit)
		   CHARACTER *(*) message
                   INTEGER iounit
                   
                   write(iounit,'(a)')message
                   Call Exit(1)
	        RETURN
	        END subroutine stop
	END MODULE terminate
