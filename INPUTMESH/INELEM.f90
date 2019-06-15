!$====================================================================================================$

!$Id: inelem.F90 527 2011-09-05 18:00:00Z LAMCE $

!$HeadURL: --- $

!$Revision: --- $

!$Date: 2011-09-05 18:00:00 -0300 (seg, 05 set 2011) $

!$Author: LAMCE $

!$====================================================================================================$

! $Id: inelem.F90 1 2011-09-05 18:00:00 LAMCE $

!***** INPUTMESH/inelem

!  NAME

!  inelem

!

!  DESCRIPTION

!  L� O TIPO DO MATERIAL, O ELEMENTO (PARA DEFINI��O DO N� DE N�S POR ELEMENTO)E L� AS CONECTIVIDADES. AL�M DISSO FAZ A
!  CONTAGEM PARA COMPARAR: N� DE PARTES LOCAL COM GLOBAL E N� DE N�S POR ELEMENTO LOCAL COM GLOBAL E COMPARA TAMB�M SE
!  O N� DE N�S DO ELEMENTO CONDIZ COM O TIPO DO MATERIAL.
 
!

!  INPUTS

!  NUMMAT (N�MERO TOTAL DE PARTES),NUME (N�MERO TOTAL DE ELEMENTOS),IIN2 (ARQUIVO DE EXTENS�O GEO),
!  IOUT (ARQUIVO DE SA�DA), NNOEL (N� DE N�S POR ELEMENTO),IERROR (ARQUIVO DE INFORMA��O DO ERRO).

!  RESULT

!  NNOEL_MESH,NNOEL, NUME_MESH, NUME, NUMPARTS, ETYPE, NUMEPART, IOUT, IERROR

!

!  USES

!  

!

!  NOTES

!  

!

!  SOURCE        
      
     
SUBROUTINE INELEM (INCID,NUME,NNOEL,IIN2,NUMMAT,IOUT,IERROR,MTYPE,ETYPE2)!,IOUT,LM,ND,ID,MTYPE,
! Variavel criada:
! NUMPART: No. total de partes.
! ETYPE: Tipo de elemento.
! NUMEPART: No. de Elementos da parte em questao.
! NUMMAT: No total de materiais (partes).
      
       IMPLICIT NONE
     
       INTEGER             :: I,N,J,INO,NUME_MESH,NUMEPART,NUMPARTS,NNOEL_MESH,NPART, ios
       INTEGER,INTENT(in)  :: NUMMAT,NUME,IIN2,IOUT,NNOEL,IERROR
       INTEGER,INTENT(out) :: INCID(NUME,NNOEL),MTYPE (NUME)
       CHARACTER*6         :: ETYPE
       CHARACTER*5         :: AUX2
       CHARACTER*13        :: AUX,ETYPE2
      

      WRITE(IOUT,2040)

! ... LEITURA DAS CONECTIVIDADES, MATERIAL E CORPO REFERENTE AOS ELEMENTOS

      NUME_MESH = 0
      NUMPARTS = 0
      
      DO
        
       READ (IIN2,'(A)', iostat=ios)  AUX
       if (ios  /= 0) exit  ! EOF
          IF(AUX.eq.'      ') GOTO 1
          READ (AUX,'(A4,I9)') AUX2,NPART
          READ (IIN2,*)
          READ (IIN2,'(A6)') ETYPE    
          READ (IIN2,'(I8)') NUMEPART
         
          SELECT CASE (ETYPE)
              CASE ('tria6')
                  NNOEL_MESH = 6
              CASE ('quad4')
                  NNOEL_MESH = 4
              CASE ('hexa8')
                  NNOEL_MESH = 8
              CASE ('hexa27')
                  NNOEL_MESH = 27
              CASE DEFAULT
                  WRITE (IOUT,2000)
           END SELECT
             
           DO I = 1, NUMEPART
              NUME_MESH = NUME_MESH + 1
              READ  (IIN2,*) N,(INCID(N,INO),INO=1,NNOEL_MESH)
             !WRITE (IOUT,*) N,(INCID(I,INO),INO=1,NNOEL),ETYPE
              MTYPE (N)= NPART           
              
           
           ENDDO
           
           
           NUMPARTS = NUMPARTS + 1

           
           
                                                     
    ENDDO
    
   
1 continue

if (etype2 .ne. 'tdpflex') then
 IF(NNOEL_MESH.NE.NNOEL) THEN
       WRITE (*,100)
       WRITE (IERROR,100)
       read (*,*)
       STOP  
    ENDIF
endif
    
    IF (NUME_MESH.NE.NUME) THEN
        WRITE (*,200)
        WRITE (IERROR,200)
        read (*,*)
        STOP  
    ENDIF
    
    IF(NUMMAT.NE.NUMPARTS) THEN
       WRITE (*,300)
       WRITE (IERROR,300)
       read (*,*)
       STOP  
    ENDIF

      RETURN
  100 FORMAT (/,/,'***********************************************************************',/,/,'ERROR &
  (INELEM.f90): NNOEL(jobname.cntr) .NE. NNOEL_MESH(jobname.mesh)',/,/,'*************************************')
      
  200 FORMAT (/,/,'***********************************************************************',/,/,'ERROR &
  (INELEM.f90): NUME(jobname.cntr) .NE. NUME_MESH(jobname.mesh)',/,/,'*************************************')
  
  300 FORMAT (/,/,'***********************************************************************',/,/,'ERROR &
  (INELEM.f90): NUMMAT(jobname.cntr) .NE. NUMEPARTS(jobname.mesh)',/,/,'*************************************')
  
 1020 FORMAT (6I8)
 2040 FORMAT (                                                          &
          ' INFORMACOES DOS ELEMENTOS ',//,                             &
          ' ____________________________________________________ ',/)
 1030 FORMAT (6I8)
 
 2000 FORMAT ('TIPO DE ELEMENTO INVALIDO.')
      END
      
     
 !********
