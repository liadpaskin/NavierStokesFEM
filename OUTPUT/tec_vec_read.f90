
      subroutine tec_vec_read (tipo,disp,up,id,ngl,numnp,neq,neqp,istep)
	  
	  use modtapes,     only: filer, iblank, iplt
      use modmesh,      only: invp
      
      implicit none
      
      integer :: i, istep, neq, neqp, ngl, numnp, id   (ngl,numnp)
      
      character*70 filename, clixo
      character*1  tipo
      real*8 :: disp (0:neq), up (0:neqp)
      real*8, allocatable ::  ap(:)
      
      allocate (ap(-neqp:neq)) ; ap = 0.d0

      filename = filer
      write (filename(1:4),'(i4.4)') istep
      write(filename,'(a)') './OUT2_init/'//filer(1:iblank-1)//'.'//tipo//filename(1:4)//'.dat'

      open (unit=iplt,file=filename)
      
      do i = 1, 9
        read (iplt,'(a)')
      enddo

           IF (ngl .eq. 3) THEN
               do i = 1, numnp
                 read (iplt,'(3e16.9)') ap(id(1,i)),ap(id(2,i)),ap(id(3,i))
                 IF (ap(id(1,i)) .ne. 0) then
                    write(*,*) 'oi'
                 endif
               enddo
           ELSE
             write (*,*) 'Error in: Subroutine ensg_vec'
             read(*,*)
             STOP
           ENDIF
       
    do i=1, neq
        disp(i) = ap(i)
    enddo
    
    deallocate (ap)

   close (unit=iplt)
      
   return
   end
