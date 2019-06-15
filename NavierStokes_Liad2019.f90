!    ---------------
      PROGRAM FEMSIMULATION
!    ---------------

      include 'mpif.h'

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

    !Serial IO input
    if l_root then
        call iomngr (0)
        call contrl
        CALL ALLOC (1)
        call inmesh
        call loads
    endif
    call MPI_Barrier(MPI_Comm comm)

    call ALLOC(3)
    call solver

    ! call output

    call iomngr (1)

    !call copyright

      te_mpi = mpi_wtime()
      write(6,9997) (te_mpi - ts_mpi)
      write(nprt,9997) (te_mpi - ts_mpi)
 9997 format(' Job Execution Time = ',e15.6)
    write (*,*)' FEMSIMULATION - END PROGRAM  '

    call mpi_finalize(ierr)

    stop
    end
