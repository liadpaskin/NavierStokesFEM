!    
      SUBROUTINE IOMNGR  ( INT )

	  use modtapes
      use modvar
	
      IMPLICIT REAL*8 (A-H,O-Z)

	  CHARACTER*70      filein1  , &
	                    filein2  , &
                        filein3  , &
                        fileerror, &
                        FILEOUT , fileplot

	logical result

      IF ( INT .EQ. 0 ) THEN

! --- DEFINICAO DAS UNIDADES LOGICAS PARA I/O

       !WRITE (*,'(A)') ' ENTER WIHT THE JOBNAME:'
       !READ  (*,'(A)')  FILEINPUT    
                                                                     
       FILEINPUT = "model01" 
        filer = FILEINPUT
 
        IBLANK = 0


	DO I = 1,80
	  IF (FILER(I:I).EQ.' ') THEN
	    IBLANK=I
	    GO TO 10 
	  ENDIF     
	ENDDO
10    CONTINUE 

        IF ( IBLANK .EQ. 0 ) THEN
          WRITE (*,'(A)') 'IOMNGR.F90 - NOME INVALIDO PARA JOBNAME ABORTANDO !'
          !PAUSE
          read (*,*)
          STOP
        ENDIF

!     MONTA OS NOMES DOS ARQUIVOS
        WRITE (FILEIN1,'(A)')   FILER(1:IBLANK-1)//'.CNTR'
        WRITE (FILEIN2,'(A)')   FILER(1:IBLANK-1)//'.GEO'
        WRITE (FILEIN3,'(A)')   FILER(1:IBLANK-1)//'.RHSV'
        
        !result = SYSTEMQQ('md OUT1')
        !result = SYSTEMQQ('md OUT2')


        WRITE (FILEERROR,'(A)') './OUT1/'//FILER(1:IBLANK-1)//'.ERRO'                
        WRITE (FILEOUT,'(A)')   './OUT1/'//FILER(1:IBLANK-1)//'.LIS'
        WRITE (FILEPLOT,'(A)')  './OUT1/'//FILER(1:IBLANK-1)//'.PLOT'

          ! write(filename,'(a)') './out/'//filer(1:iblank-1)//'.'//tipo//filename(1:4)

        
        WRITE (*,100)
        WRITE (*,101)
        WRITE (*,*) FILEIN1
        WRITE (*,102)
        WRITE (*,*) FILEIN2
        WRITE (*,103)
        WRITE (*,*) FILEIN3

        WRITE (*,200)
        WRITE (*,104)
        WRITE (*,*) FILEOUT
        WRITE (*,105)
        
        WRITE (*,*) FILEPLOT


100     FORMAT (//, 'DATA INPUT MODULE:',//) 
101     FORMAT (/,'** CONTROL DATA  ** (JOBNAME.CNTR)')        
102     FORMAT (/,'** MESH DATA     ** (JOBNAME.MESH)')        
103     FORMAT (/,'** BOUNDARY DATA ** (JOBNAME.RHSV)') 

200     FORMAT (//, 'DATA OUTPUT MODULE:',//) 
104     FORMAT (/,'** LIST DATA     ** (JOBNAME.LIS) ')      

105     FORMAT (/,'** PLOT DATA     ** (JOBNAME.PLOT) ')      
        
        IIN1    = 1
        IIN2    = 2
        IIN3    = 3
        
        IOUT    = 4
        
        IERROR  = 5
        
        IPLOT   = 61

        IPLT    = 7

      

! --------------------------------------------------------------------

        OPEN (UNIT=IIN1 ,FILE=FILEIN1)
        OPEN (UNIT=IIN2 ,FILE=FILEIN2)
        OPEN (UNIT=IIN3 ,FILE=FILEIN3)
        OPEN (UNIT=IERROR ,FILE=FILEERROR)
        OPEN (UNIT=IOUT,FILE=FILEOUT) 
        OPEN (UNIT=IPLOT,FILE=FILEPLOT) 

      ELSE

        PRINT *,' FECHA ',IIN1
        CLOSE (UNIT=IIN1   )
        PRINT *,' FECHA ',IIN2
        CLOSE (UNIT=IIN2   )
        PRINT *,' FECHA ',IIN3
        CLOSE (UNIT=IIN3   )
        PRINT *,' FECHA ',IOUT
        CLOSE (UNIT=IOUT  )
        PRINT *,' FECHA ',IERROR
        CLOSE (UNIT=IERROR)
        PRINT *,' FECHA ',IPLOT
        CLOSE (UNIT=IPLOT)
      ENDIF


      RETURN
      END
