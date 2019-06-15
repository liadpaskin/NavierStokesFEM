!    ---------------
      PROGRAM FEMSIMULATION
!    ---------------

      use modMPIvar
      use modvar

      integer :: inoi,inof,ieli,ielf
      integer, allocatable :: incid0(:,:)
      real*8, allocatable :: x0(:), y0(:), z0(:)

!
! **********************************************************************
! *                                    *
! *                           FEMSIMULATION                    *
! *                                                                    *
! *             FINITE ELEMENT SIMULATION By LIAD/KADU            *
! **********************************************************************
! COPPE-PEC/LAMCE/UFRJ/2011 JLDA/KADU/IC'S

    PRINT *,'hello world MPI'

    call mpi_init(ierr)
    call mpi_comm_rank(mpi_comm_world,myid,ierr)
    call mpi_comm_size(mpi_comm_world,numprocs,ierr)

      i_root = 0
      l_root = .false.
      if(myid .eq. i_root) l_root = .true.

      l_debug = .false.
      if(idebug .eq. 1) l_debug = .true.

      ts_mpi = mpi_wtime()

    call iomngr (0)
    call contrl
    CALL ALLOC (1)
    call inmesh

    call MPI_Barrier(MPI_Comm_World,ierr)

   ! call loads

   call partition

   ! call ALLOC(3)
   ! call solver

    if (myID == 0) then
      allocate (x0(numnp))
      allocate (y0(numnp))
      allocate (z0(numnp))
      allocate (incid0(nume,nnoel))
      do iproc = 1, numprocs-1
          ieli = procsel(iproc)
          ielf = procsel(iproc+1)-1
          inoi = procsnp(iproc)
          inof = procsnp(iproc+1)-1
          call MPI_Recv(incid0(1:ielf-ieli+1,:), (numef-numei+1)*NNOEL_MESH, MPI_REAL*8, iproc, 4, MPI_Comm_world, status)
          call MPI_Recv(x0(inoi:inof), numnpf-numnpi+1, MPI_REAL*8, iproc, 1, MPI_Comm_world, status,ierr)
          call MPI_Recv(y0(inoi:inof), numnpf-numnpi+1, MPI_REAL*8, iproc, 2, MPI_Comm_world, status,ierr)
          call MPI_Recv(z0(inoi:inof), numnpf-numnpi+1, MPI_REAL*8, iproc, 3, MPI_Comm_world, status,ierr)
      enddo
      call OUTTECPLT3D(x0,y0,z0,x0,up,id,incid0,ngl,numnp,nume,nnoel,neq,neqp,nd,kensight,T,TFUNC,'U')
    else
          call MPI_Send(incid, (numef-numei+1)*NNOEL_MESH, MPI_REAL*8, 0, 4, MPI_Comm_world, status)
          call MPI_Send(x, numnpf-numnpi+1, MPI_REAL*8, 0, 1, MPI_Comm_world, status,ierr)
          call MPI_Send(y, numnpf-numnpi+1, MPI_REAL*8, 0, 2, MPI_Comm_world, status,ierr)
          call MPI_Send(z, numnpf-numnpi+1, MPI_REAL*8, 0, 3, MPI_Comm_world, status,ierr)
    endif

    ! call output

   ! call iomngr (1)

    !call copyright

      te_mpi = mpi_wtime()
      write(6,9997) (te_mpi - ts_mpi)
      write(nprt,9997) (te_mpi - ts_mpi)
 9997 format(' Job Execution Time = ',e15.6)
    write (*,*)' FEMSIMULATION - END PROGRAM  '

    call mpi_finalize(ierr)

    stop
    end
