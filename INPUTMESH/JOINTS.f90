!$====================================================================================================$

!$Id: joints.F90 527 2011-08-30 12:18:04Z LAMCE $

!$HeadURL: --- $

!$Revision: --- $

!$Date: 2011-08-30 19:36:00 -0300 (ter, 30 ago 2011) $

!$Author: LAMCE $

!$====================================================================================================$

! $Id: joints.F90 1 2011-08-30 19:36:00 Fabio, Vitor e Daniel $

!****x* INPUTMESH/joints

!  NAME

!  joints

!

!  DESCRIPTION

!  l� as coordenadas de cada n� da malha
 
!

!  INPUTS

!  NUMNP

!

!  RESULT

!  NUMNP_MESH, X, Y, Z

!

!  USES

!  modvar, modmesh

!

!  NOTES

!  

!

!  SOURCE      
      
      
      SUBROUTINE JOINTS (X,Y,Z,NUMNP,IIN2,IOUT,IERROR)     
	
      
      IMPLICIT none


      INTEGER :: NUMNP_MESH, I, INO        
          
      REAL*8 , INTENT(out)   ::   X  (NUMNP)          , &
                                  Y  (NUMNP)          , &
                                  Z  (NUMNP)

	  INTEGER , INTENT(in)   ::   NUMNP               , &
	                              IIN2                , &
	                              IOUT                , &
	                              IERROR	

	  CHARACTER*20  :: teste
                            

!_______________________________________________________________________
!
!     LEITURA DOS DADOS FORNECIDOS
!
!      N   ===   NUMERO DO NO
!
!      X   ===   COORDENADA X DO NO
!      Y   ===   COORDENADA Y DO NO
!      Z   ===   COORDENADA Z DO NO
!
!_______________________________________________________________________
!


      READ (IIN2,'(A)') teste
      READ (IIN2,'(A)') teste
      READ (IIN2,'(A)') teste
      READ (IIN2,'(A)') teste
      READ (IIN2,'(A)') teste
          
      
      READ (IIN2,'(I8)') NUMNP_MESH
      
      
      IF (NUMNP.NE.NUMNP_MESH) THEN
         WRITE (*,100)
         WRITE (IERROR,100)
         read (*,*)
         STOP
      ENDIF
   
   
      WRITE (IOUT,200)


      DO I = 1, NUMNP    
         READ  (IIN2,300) INO, X(INO), Y(INO), Z(INO)
      ENDDO
      

      
      RETURN
      

100  FORMAT (/,/,&
             '******************************************************************',/,&
             /,'ERROR(JOINTS.f90): NUMNP(jobname.cntr) .NE. NUMNP_MESH(jobname.mesh)',/,/,&
             '******************************************************************') 
200  FORMAT (' DADOS DOS PONTOS NODAIS',/)
300  FORMAT (I8,3E12.5)


     END
