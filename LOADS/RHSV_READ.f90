SUBROUTINE RHSV_READ (NUMNP,NGL,IIN3,IOUT,IERROR,invp,ireordering,fno,found)
    
use modmesh, only: x, y, prop, mtype     ! Bacalho para gerar RHSV
use modvar, only: grav
    
IMPLICIT none
   
logical*1,INTENT(out)  :: found(ngl+1,numnp)
real*8   ,INTENT(out)    :: Fno (NUMNP,NGL+1,2)


logical*1 :: BF, IC
integer :: Ino, N, idof, ngl, NUMNP, BC, IIN3,  IOUT, IERROR, I, J, K,invp(numnp), ireordering,fin
real*8    :: v2 , xmax, xRho
character*70 :: line
character*4  :: ptype

if (ireordering .eq. 0) then 
    do i = 1 , numnp
        invp(i) = i
    enddo
endif

!xRho    = prop(mtype(1),4)       ! Densidade
!xmax = 0.d0
!write (666,'(A)') '*BC1'
!do i=1, numnp                                               !
!    if (x(i) .gt. xmax) xmax = x(i)
    !if (x(i) .eq. 4.0d0) then                               !
    !    write (666,*) i                                      !                           !
    !endif
!    if (y(i) .eq. -1.0d0) then                               !
!        write (666,*) i, 1, 0.0d0                           ! Bacalho para gerar RHSV
!        write (666,*) i, 2, 0.0d0                           !                                                !
!    elseif (x(i) .eq. 0.0d0) then                               !
!        write (666,*) i, 1, 1.0d0                           ! Bacalho para gerar RHSV
!        write (666,*) i, 2, 1.0d0                           !
    !elseif (x(i) .eq. 2.0d0) then                               !
    !    write (666,*) i, 1, 0.0d0                           ! Bacalho para gerar RHSV
    !    write (666,*) i, 2, 0.0d0                           !
    !elseif (y(i) .eq. 0.0d0) then                               !
    !    write (666,*) i, 1, 0.0d0                           ! Bacalho para gerar RHSV
    !    write (666,*) i, 2, 0.0d0                           !   
!    endif
!enddo
!write (666,'(A)') '*BC2'
!do i=1, numnp                                               !
!    if (x(i) .eq. xmax) then                               !
!        write (666,*) i, 1, xRho*grav*y(i)                                      !                           !
!    endif
!enddo                                                        !

DO

    READ (IIN3,'(a)', iostat = fin) line
    if (fin /= 0) exit

    IF (line(1:1) .eq. '*') THEN
        READ (line,'(a)') ptype
        SELECT CASE (ptype)
            CASE ('*BC1')
                BC = 1
                BF = .false.
                IC = .false.
            CASE ('*BC2')
                BC = 2
                BF = .false.
                IC = .false.
            CASE ('*BF ')
                BF = .true.
                IC = .false.
            CASE ('*IC ')
                IC = .true.
                BF = .false.
         ENDSELECT
    ELSEIF (line(1:1) .eq. "n") then
    write(*,*) "sem condi��o de contorno"
 
    ELSE
        
        READ (line,*) Ino, idof, v2
        
        IF (BC .eq. 1) THEN
            found(idof,invp(Ino)) = .true.
        ENDIF
        
        Fno (invp(Ino), idof, BC) = v2
                
    ENDIF
              
ENDDO

RETURN

ENDSUBROUTINE
