
subroutine cuthillMcKee()

      use    modvar,    only: numnp,nume,nummat,nnoel,ncprop,etype,ireordering,nnomax,nad,isolver_type,ntria,nbar
      use    modmesh,   only: x, y, z, incid, incidaux2, prop, mtype , jdiag , ip , mask , xls  , invp , invp2 ,list_tria3
      use    modtapes,  only: iin1, iin2, iout, ierror


      IMPLICIT NONE
      integer:: i, jnoel, iel, ino, k , i0, i1, i2, i3, i4, i5, i6, itwo, imod ,nume_aux ,numnp_aux , nume1, ino_old, ino_new
      real*8  , allocatable :: xr(:) , yr(:) , zr (:)
      integer , allocatable :: incidaux(:,:) , noviz (:) , ips(:)

          numnp_aux = numnp
          nume_aux = nume
          nume1 = 1

          allocate (incidaux2         (nume,nnoel))   ; incidaux2 (:,:) = 0
          allocate (xr(numnp))
          allocate (yr(numnp))
          allocate (zr(numnp))
          allocate (incidaux         (nnoel,nume))   ; incidaux (:,:) = 0
          incidaux2(:,:) = incid(:,:)
          nnomax = nnoel ! numero max. de nos por elemento
          itwo = 2
          i0 = 1
          i1 = i0 + nume*(nnomax+1)
          imod = mod(i1-1,itwo)
          i1 = i1 + imod
          i2 = i1 + numnp*3*itwo
          i3 = i2 + numnp + 1
          i4 = i3 + numnp

          call matrixTrans2 (incid,nume,nnoel,incidaux) ! Retorna atraves de incidAUX a transposta da matriz contida no parametro incid.
! PROFIL : determina o perfil da matriz de rigidez
! INPUT: numnp  ;nume  ;nnomax;  incidaux = incidencia dos elementos
          call profil  (jdiag, incidaux , nnomax , numnp, nume , nad , nume1)
! OUTPUT: jdiag = vetor apontador do perfil da matriz ;   nad  = numero de coeficientes do perfil da matriz
          call aponta  (incidaux,ip,numnp,nume,nnomax,i5,nnoel) ! ALTERADO POR LIAD

          i6 = i5
          i5 = i4 + i5

          allocate (noviz                    (i6))   ; noviz (:) = 0
          allocate (ips                      (i6))   ; ips (:) = 0

          call grafos  (incidaux,ip,ips,noviz,nume,nnomax,numnp,i3,i6,nume1,nnoel) ! ALTERADO POR LIAD
          call sort    (ip,ips,noviz,numnp,i3,i4,i5,i6) ! ???

          !INPUT: numnp, ip, ips
          call genrcm  (numnp, ip, ips, noviz, mask, xls, i6, numnp)   ! general reverse cuthill mckee
          ! OUTPUT: noviz

          call splot (incidaux , noviz , invp , nnomax, numnp, nume , i6,nume1,nnoel) ! ALTERADO POR LIAD ! grava o arquivo para plotter.
          call matrixtrans2 (incidaux,nnoel,nume,incid)! Retorna atraves de incid a transposta da matriz contida no parametro incidAUX.

          deallocate (noviz)
          deallocate (ips)

          xr(:) = 0
          yr(:) = 0
          zr(:) = 0

          do i = 1 , numnp
               xr(invp(i)) = x(i)
               yr(invp(i)) = y(i)
               zr(invp(i)) = z(i)
          enddo

          x(:) = xr(:)
          y(:) = yr(:)
          z(:) = zr(:)

          deallocate(xr)
          deallocate(yr)
          deallocate(zr)
          deallocate (incidaux)

      return
  end
