 
SUBROUTINE JOINTS (X,Y,Z,NUMNP,IIN2,IOUT,IERROR)
	
      use modMPIvar

      implicit none

      INTEGER :: NUMNP_MESH, I, INO

      REAL*8 , allocatable   ::   X0  (:)          , &
                                  Y0  (:)          , &
                                  Z0  (:)
          
      REAL*8 , INTENT(out)   ::   X  (NUMNP_proc)          , &
                                  Y  (NUMNP_proc)          , &
                                  Z  (NUMNP_proc)

	  INTEGER , INTENT(in)   ::   NUMNP               , &
	                              IIN2                , &
	                              IOUT                , &
	                              IERROR	

	  CHARACTER*20  :: teste
       integer :: inoi, inof, iproc

	  if (myID == 0) then

          allocate(X0(NUMNP_proc))
          allocate(Y0(NUMNP_proc))
          allocate(Z0(NUMNP_proc))

          READ (IIN2,'(A)') teste
          READ (IIN2,'(A)') teste
          READ (IIN2,'(A)') teste
          READ (IIN2,'(A)') teste
          READ (IIN2,'(A)') teste

          READ (IIN2,'(I8)') NUMNP_MESH
      
          IF (NUMNP.NE.NUMNP_MESH) THEN
             WRITE (*,100)
             WRITE (IERROR,100)
             read (*,*)
             STOP
          ENDIF

          WRITE (IOUT,200)

          DO I = numnpi, numnpf
              READ  (IIN2,300) INO, x(i),y(i),z(i)
          ENDDO
          do iproc = 1, numprocs-1
              inoi = procsnp(iproc)
              inof = procsnp(iproc+1)-1
              DO I = inoi, inof
                 READ  (IIN2,300) INO, x0(i-inoi+1),y0(i-inoi+1),z0(i-inoi+1)
              ENDDO
              call MPI_Send(x0(1:(inof-inoi+1)), inof-inoi+1, MPI_DOUBLE_PRECISION, iproc, 1, MPI_Comm_world,ierr)
              call MPI_Send(y0(1:(inof-inoi+1)), inof-inoi+1, MPI_DOUBLE_PRECISION, iproc, 2, MPI_Comm_world,ierr)
              call MPI_Send(z0(1:(inof-inoi+1)), inof-inoi+1, MPI_DOUBLE_PRECISION, iproc, 3, MPI_Comm_world,ierr)
          enddo

          deallocate(X0)
          deallocate(Y0)
          deallocate(Z0)

      else

      call MPI_Recv(x, numnpf-numnpi+1, MPI_DOUBLE_PRECISION, 0, 1, MPI_Comm_world, status,ierr)
      call MPI_Recv(y, numnpf-numnpi+1, MPI_DOUBLE_PRECISION, 0, 2, MPI_Comm_world, status,ierr)
      call MPI_Recv(z, numnpf-numnpi+1, MPI_DOUBLE_PRECISION, 0, 3, MPI_Comm_world, status,ierr)

      ENDIF
      
      RETURN

100  FORMAT (/,/,&
             '******************************************************************',/,&
             /,'ERROR(JOINTS.f90): NUMNP(jobname.cntr) .NE. NUMNP_MESH(jobname.mesh)',/,/,&
             '******************************************************************') 
200  FORMAT (' DADOS DOS PONTOS NODAIS',/)
300  FORMAT (I8,3E12.5)

     END
