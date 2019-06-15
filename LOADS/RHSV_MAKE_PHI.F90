SUBROUTINE RHSV_MAKE_PHI (LM,INCID,ID,NNOEL,NUME,NUMNP,NEQ,NEQP,NGL,IIN3,IOUT,IERROR,fno,found,invp)

IMPLICIT none
       
INTEGER,INTENT(out)  :: LM(nume,nnoel), ID(1,numnp)
integer,intent(in)   :: incid(nume,nnoel)

logical*1 :: BF, IC , found(NGL+1,NUMNP)
integer :: Ino, N, idof, ngl, NUMNP, BC, IIN3, IOUT, IERROR, NUME, nnoel, I, J, K, icount , neq , neqp , invp(numnp)
real*8    :: v2 , Fno (NUMNP,NGL+1,2)
character*70 :: line
character*4  :: ptype

icount = 0

DO I=1,numnp
    DO J=ngl+1,ngl+1
        IF  (.not. found(J,I)) THEN
            neq = neq + 1
            ID (1,I) = neq
        ELSE
           IF (fno(i,j,1) .eq. 0.d0) THEN
                ID (1,I) = 0
           ELSE
                neqp = neqp + 1
                ID (1,I) = -neqp
           ENDIF
        ENDIF
    ENDDO
ENDDO

DO I=1,nume
    DO J=1,nnoel
        icount = icount + 1
        LM(I,icount) = ID (1,incid(I,J))
    ENDDO
    icount = 0
ENDDO

RETURN

ENDSUBROUTINE