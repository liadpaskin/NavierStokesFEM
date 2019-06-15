  subroutine initialInterpolate (tipo,disp,up,id,ngl,numnp,neq,neqp,istep)
	  
	  use modtapes,     only: filer, iblank, iplt, IOUT, ierror
      use modmesh,      only: invp, x, y, z, mtype
      use modvar,    only: nume,nummat,nnoel,ncprop,etype
      
      implicit none
      
      real*8              :: Nf(8)
      REAL*8              :: x15, y15, z15,x14, y14, z14, x12, y12, z12, DetJ, DetJinv, Jinv(3,3), xp, yp, zp, r, s, t, x1, y1, z1
      real*8, allocatable :: x0(:), y0(:), z0(:),u0(:,:)
      integer, allocatable:: incid0(:,:)
      integer :: numnp0, nume0, ino, ino0, iel0, inoel, nnoel0

      integer :: istep, neq, neqp, ngl, numnp, id(ngl,numnp), marc, iin0, i
      character*70 filename, FILEIN0
      character*1  tipo
      real*8 :: disp (0:neq), up (0:neqp)
      real*8, allocatable ::  ap(:)

      numnp0 = 35301
      nume0 = 4000*8
      nnoel0 = 8

      allocate (incid0     (nume0,nnoel0));                       incid0(:,:) = 0
      allocate (x0              (numnp0));                      x0(:) = 0.d0 
      allocate (y0              (numnp0));                      y0(:) = 0.d0 
      allocate (z0              (numnp0));                      z0(:) = 0.d0  
      allocate (u0              (numnp0,ngl));                  u0(:,:)  = 0.d0
      allocate (ap(-neqp:neq)) ; ap = 0.d0
      
      IIN0 = 100
      WRITE (FILEIN0,'(A)')   FILER(1:IBLANK-1)//'_init.GEO'
      OPEN (UNIT=IIN0 ,FILE=FILEIN0)

      filename = filer
      write (filename(1:4),'(i4.4)') istep
      write(filename,'(a)') './OUT2_init/'//filer(1:iblank-1)//'.'//tipo//filename(1:4)//'.dat'
      open (unit=iplt,file=filename)
      do i = 1, 12
        read (iplt,'(a)')
      enddo
           IF (ngl .eq. 3) THEN
               do i = 1, numnp0
                 read (iplt,'(6e16.9)') x0(i),y0(i),z0(i),u0(i,1),u0(i,2),u0(i,3)
               enddo
               do i = 1, nume0
                 read (iplt,*) incid0(i,1),incid0(i,2),incid0(i,3),incid0(i,4),incid0(i,5),incid0(i,6),incid0(i,7),incid0(i,8)
               enddo
           ELSE
             write (*,*) 'Error in subroutine initialInterpolate'
             read(*,*)
             STOP
           ENDIF
           
      do ino = 1,numnp
      	  if (id(1,ino) < 1) then
      	  		cycle
      	  endif
	      xp = x(ino)
	      yp = y(ino)
	      zp = z(ino)
	      !marc = 0
          !do ino0 = 1,numnp0
          !      if (( abs(x(ino)-x0(ino0)) < 0.00001d0) .and. (abs(y(ino)-y0(ino0)) < 0.00001d0) .and. (abs(z(ino)-z0(ino0)) < 0.00001d0)) then
          !            ap(id(1,ino)) = u0(ino0,1)
          !            ap(id(2,ino)) = u0(ino0,2)
          !            ap(id(3,ino)) = u0(ino0,3)
          !            marc = 1
          !            exit
          !      endif
          !enddo
          !if (marc .eq. 1) then
          !      cycle
          !endif
          do iel0 = 1,nume0
		      	
               x1 = x0(incid0(iel0,1))
               y1 = y0(incid0(iel0,1))
               z1 = z0(incid0(iel0,1))

               if ((abs(x1-xp)>1E-01) .or. (abs(y1-yp)>1E-01) .or. (abs(z1-zp)>1E-01)) then
                    cycle
               endif

		       x12 = x0(incid0(iel0,2)) - x0(incid0(iel0,1))  !
		       x14 = x0(incid0(iel0,4)) - x0(incid0(iel0,1))  ! Vetores que definem o hexaedro
		       x15 = x0(incid0(iel0,5)) - x0(incid0(iel0,1))  !
		       y12 = y0(incid0(iel0,2)) - y0(incid0(iel0,1))  !
		       y14 = y0(incid0(iel0,4)) - y0(incid0(iel0,1))  ! Vetores que definem o hexaedro
		       y15 = y0(incid0(iel0,5)) - y0(incid0(iel0,1))  !
		       z12 = z0(incid0(iel0,2)) - z0(incid0(iel0,1))  !
		       z14 = z0(incid0(iel0,4)) - z0(incid0(iel0,1))  ! Vetores que definem o hexaedro
		       z15 = z0(incid0(iel0,5)) - z0(incid0(iel0,1))  !
           
		      detJ = z12 * x14 * y15 - z12 * x15 * y14 - y12 * x14 * z15 + y12 * x15 * z14 + x12 * y14 * z15 - x12 * y15 * z14
		      detJinv = 1.d0/detJ
		      
		      Jinv(1,1) = (y14 * z15 - y15 * z14)*detJinv
		      Jinv(1,2) = (-x14 * z15 + x15 * z14)*detJinv
		      Jinv(1,3) = (x14 * y15 - x15 * y14)*detJinv
		      Jinv(2,1) = (-y12 * z15 + y15 * z12)*detJinv
		      Jinv(2,2) = (x12 * z15 - z12 * x15)*detJinv
		      Jinv(2,3) = (-x12 * y15 + x15 * y12)*detJinv
		      Jinv(3,1) = (y12 * z14 - y14 * z12)*detJinv
		      Jinv(3,2) = (-x12 * z14 + z12 * x14)*detJinv
		      Jinv(3,3) = (x12 * y14 - y12 * x14)*detJinv
      
		      r = Jinv(1,1)*(xp-x1)+Jinv(1,2)*(yp-y1)+Jinv(1,3)*(zp-z1)
		      s = Jinv(2,1)*(xp-x1)+Jinv(2,2)*(yp-y1)+Jinv(2,3)*(zp-z1)
		      t = Jinv(3,1)*(xp-x1)+Jinv(3,2)*(yp-y1)+Jinv(3,3)*(zp-z1)
		      
		      if ((r >= 0 .and. r <= 1) .and. (s >= 0 .and. s <= 1) .and. (t >= 0 .and. t <= 1)) then

                         Nf(1) = 1 - r - s - t + r*s + s*t + t*r - r*s*t   !
                         Nf(2) = r - r*s - t*r + r*s*t                     !
                         Nf(3) = r*s - r*s*t                               !
                         Nf(4) = s - r*s - s*t + r*s*t                     !   Funcoes de forma para o octa 8
                         Nf(5) = t - s*t - t*r + r*s*t                     !
                         Nf(6) = t*r - r*s*t                               !
                         Nf(7) = r*s*t                                     !
                         Nf(8) = s*t - r*s*t                               !

						do inoel = 1, nnoel0
							ap(id(1,ino)) = ap(id(1,ino)) + Nf(inoel)*u0(incid0(iel0,inoel),1)
							ap(id(2,ino)) = ap(id(2,ino)) + Nf(inoel)*u0(incid0(iel0,inoel),2)
							ap(id(3,ino)) = ap(id(3,ino)) + Nf(inoel)*u0(incid0(iel0,inoel),3)
					    enddo

					    exit

		      endif
		 enddo
	enddo

    do i=1, neq
        disp(i) = ap(i)
    enddo

      deallocate (incid0)
      deallocate (x0    )
      deallocate (y0    )
      deallocate (z0    )
      deallocate (u0    )
      deallocate (ap    )

      return
      end
