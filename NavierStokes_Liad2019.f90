!    ---------------
      PROGRAM FEMSIMULATION
!    ---------------

!
! **********************************************************************
! *                                    *
! *                           FEMSIMULATION                    *
! *                                                                    *
! *             FINITE ELEMENT SIMULATION By LIAD/KADU            *
! **********************************************************************
! COPPE-PEC/LAMCE/UFRJ/2011 JLDA/KADU/IC'S

    PRINT *,'hello world'

    call iomngr (0)

    call contrl

    CALL ALLOC (1)

    call inmesh

    call loads

    call ALLOC(3)

    call solver

    ! call output

    call iomngr (1)

    !call copyright

    write (*,*)' FEMSIMULATION - END PROGRAM  '
    stop
    end
