      SUBROUTINE INMESH
	
!     LEITURA DOS DADOS DA MALHA E DEFINICAO DE VARIOS ARRANJOS

      use    modvar,    only: numnp,nume,nummat,nnoel,ncprop,etype,ireordering,isolver_type
      use    modmesh,   only: x, y, z, incid, prop, mtype
      use    modtapes,  only: iin1, iin2, iout, ierror
      use modMPIvar

      implicit none

      integer :: numei0, numef0, numnpi0, numnpf0, iproc, i
       
      do iproc = 0, numprocs
              numnpi0 = (numnp/numprocs)*iproc+1
              numei0 = (nume/numprocs)*iproc+1
              if (iproc .lt. mod(numnp,numprocs)) then
                  numnpi0 = numnpi0 + iproc
              else
                  numnpi0 = numnpi0 + mod(numnp,numprocs)
              endif
              if (iproc .lt. mod(nume,numprocs)) then
                  numei0 = numei0 + iproc
              else
                  numei0 = numei0 + mod(nume,numprocs)
              endif
              procsnp(iproc) = numnpi0
              procsel(iproc) = numei0
      enddo

      CALL JOINTS    (X,Y,Z,NUMNP,IIN2,IOUT,IERROR)
      write(*,*)
      write(*,*) '-----------------------------'
      write(*,*) 'DEBUG PROCCESSOR: ', myID
      do i = numnpi,numnpf
        write(*,*) i, x(i), y(i), z(i)
      enddo
      write(*,*) '-----------------------------'
    call MPI_Barrier(MPI_Comm_World,ierr)

      CALL MATSET    (PROP,NUMMAT,NCPROP,IIN1,IOUT,IERROR,ETYPE,isolver_type)
      CALL INELEM    (INCID,NUME,NNOEL,IIN2,NUMMAT,IOUT,IERROR,MTYPE,etype)
      write(*,*)
      write(*,*) '-----------------------------'
      write(*,*) 'DEBUG PROCCESSOR: ', myID
      do i = numei,numef
        write(*,*) i, incid(i,1),incid(i,2),incid(i,3),incid(i,4),incid(i,5),incid(i,6),incid(i,7),incid(i,8)
      enddo
      write(*,*) '-----------------------------'
      
!      if (ireordering .eq. 1 ) then
!            call CuthillMcKee()
!      endif

      RETURN
      END


