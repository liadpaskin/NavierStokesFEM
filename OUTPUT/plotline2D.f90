subroutine PlotLine2D(u, up, id, neq, neqp, ngl, num)

use modtapes, only: filer, iblank
use modmesh, only: x, y
use modvar, only: numnp

implicit none

! global
integer, intent(in) :: neq, ngl, neqp , id(ngl,numnp), num
real*8, intent(in)  :: up(0:neqp), u(0:neq)

!      Local
character*70 :: filename, filename1
integer      :: i, j, ieq, xcont, ino, aux, ilin
real*8, allocatable ::  ap(:)
integer, allocatable::  l_vert(:)

allocate (ap(-neqp:neq)) ; ap = 0.d0
allocate (l_vert(numnp)) ; l_vert = 0.d0

filename = filer
write (filename(1:4),'(i4.4)') num
write(filename1,'(a)') './OUT2/'//filer(1:iblank-1)//'_PLOTline.'//filename(1:4)
open (unit=202,file=filename1,form='formatted')
write (202,'(a,i5,1x,a)') 'Saida para Driven Cavity, linha de centro vertical, passo: ', 0000


do i = 1, numnp
  do j = 1, ngl
     ieq = id(j,i)
     if (ieq.ge.0) then
       ap(ieq)=u(ieq)
     else
       ap(ieq)=up(-ieq)
     endif
  enddo
enddo

xcont = 0
do ino = 1, numnp
    if ( (x(ino) .eq. 0.5d0)) then
        xcont = xcont + 1
        l_vert    (xcont) = ino
    endif
enddo

do I=1,xcont-1
        do J= I+1,xcont
            if (y(l_vert(I)) .gt. y(l_vert(J)) ) then
                    aux=l_vert(I)
                    l_vert(I)=l_vert(J)
                    l_vert(J)=aux
            endif
        enddo
enddo

do i = 1, xcont
    write (202,'(e12.4E3, A, e12.4E3)') y(l_vert(i)), ' ', ap(ID(1,l_vert(i)))
enddo

CLOSE (UNIT=202  )

end
