SUBROUTINE RHSV_LOAD_PHI (ID,NUMNP,NEQ,NEQP,NGL,IIN3,IOUT,IERROR,fno)

use modloads,  only: f_PHI, up_PHI
use modsolver, only: SedPHI
IMPLICIT none
       
integer,intent(in)   :: ID(1,numnp)

integer :: Ino, N, idof, ngl, NUMNP, BC, IIN3,  IOUT, IERROR,  I, J, K, neq , neqp
real*8    :: v2 ,Fno (NUMNP,NGL+1,2)
character*70 :: line
character*4  :: ptype

allocate (F_PHI              (0:neq));                      F_PHI(:) = 0.d0   
allocate (Up_PHI            (0:neqp));                     Up_PHI(:) = 0.d0 
allocate (SedPHI            (neqp,2));                     SedPHI(:,:) = 0.d0

DO I=1,numnp
    DO J=1,1
        IF (ID(J,I) .lt. 0) THEN
            Up_PHI(-ID(J,I)) = Fno(I,ngl+1,1)
        ELSEIF (ID(J,I) .gt. 0) THEN
            F_PHI(ID(J,I)) = Fno(I,ngl+1,2)
        ENDIF
    ENDDO
ENDDO 

RETURN

ENDSUBROUTINE