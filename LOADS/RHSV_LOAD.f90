SUBROUTINE RHSV_LOAD (ID,NUMNP,NEQ,NEQP,NGL,IIN3,IOUT,IERROR,fno)

use modloads,  only: f, up , fbody
use modvar, only: nbody
IMPLICIT none
       
integer,intent(in)   :: ID(ngl,numnp)

integer :: Ino, N, idof, ngl, NUMNP, BC, IIN3,  IOUT, IERROR,  I, J, K, neq , neqp
real*8    :: v2 ,Fno (NUMNP,NGL+1,2)
character*70 :: line
character*4  :: ptype

allocate (F              (0:neq));                      F(:) = 0.d0   
allocate (Up            (0:neqp));                     Up(:) = 0.d0 
allocate (fbody          (0:neq));                  fbody(:) = 0.d0

DO I=1,numnp
    DO J=1,ngl
        IF (ID(J,I) .lt. 0) THEN
            Up(-ID(J,I)) = Fno(I,J,1)
        ELSEIF (ID(J,I) .gt. 0) THEN
            F(ID(J,I)) = Fno(I,J,2)
        ENDIF
    ENDDO
ENDDO

!if (nbody .gt. 0) then
!    call rigidbodyread
!endif

RETURN

ENDSUBROUTINE
