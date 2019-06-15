SUBROUTINE RHSV_MAKE (LM,INCID,ID,NNOEL,ND,NUME,NUMNP,NEQ,NEQP,NGL,IIN3,IOUT,IERROR,fno,found,invp)


IMPLICIT none
       
INTEGER,INTENT(out)  :: LM(nume,nd), ID(ngl,numnp)
integer,intent(in)   :: incid(nume,nnoel)

logical*1 :: BF, IC , found(NGL+1,NUMNP)
integer :: Ino, N, idof, ngl, NUMNP, BC, IIN3, ND, IOUT, IERROR, NUME, nnoel, I, J, K, icount , neq , neqp , invp(numnp)
real*8    :: v2 , Fno (NUMNP,NGL+1,2)
character*70 :: line
character*4  :: ptype

icount = 0

DO I=1,numnp
    DO J=1,ngl
        IF  (.not. found(J,I)) THEN
            neq = neq + 1
            ID (J,I) = neq
        ELSE
           IF (fno(i,j,1) .eq. 0.d0) THEN
                ID (J,I) = 0
           ELSE
                neqp = neqp + 1
                ID (J,I) = -neqp
            ENDIF
        ENDIF
    ENDDO
ENDDO

DO I=1,nume
    DO J=1,nnoel
        DO K=1,ngl
            icount = icount + 1
            LM(I,icount) = ID (K,incid(I,J))
        ENDDO
    ENDDO
    icount = 0
ENDDO

RETURN

ENDSUBROUTINE