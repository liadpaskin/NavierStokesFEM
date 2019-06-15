     
SUBROUTINE INELEM (INCID,NUME,NNOEL,IIN2,NUMMAT,IOUT,IERROR,MTYPE,ETYPE2)!,IOUT,LM,ND,ID,MTYPE,
! Variavel criada:
! NUMPART: No. total de partes.
! ETYPE: Tipo de elemento.
! NUMEPART: No. de Elementos da parte em questao.
! NUMMAT: No total de materiais (partes).
      
      use modMPIvar

       IMPLICIT NONE
     
       INTEGER             :: I,N,J,INO,NUME_MESH,NUMEPART,NUMPARTS,NNOEL_MESH,NPART, ios
       INTEGER,INTENT(in)  :: NUMMAT,NUME,IIN2,IOUT,NNOEL,IERROR
       INTEGER,INTENT(out) :: INCID(NUME_proc,NNOEL),MTYPE (NUME)
       CHARACTER*6         :: ETYPE
       CHARACTER*5         :: AUX2
       CHARACTER*13        :: AUX,ETYPE2

       integer :: ieli, ielf, iel, iproc
       integer, allocatable :: incid0(:,:)
      

      WRITE(IOUT,2040)

! ... LEITURA DAS CONECTIVIDADES, MATERIAL E CORPO REFERENTE AOS ELEMENTOS

      NUME_MESH = 0
      NUMPARTS = 0
      

    if (myID == 0) then

      DO

           READ (IIN2,'(A)', iostat=ios)  AUX
           if (ios  /= 0) exit  ! EOF
              IF(AUX.eq.'      ') GOTO 1
              READ (AUX,'(A4,I9)') AUX2,NPART
              READ (IIN2,*)
              READ (IIN2,'(A6)') ETYPE
              READ (IIN2,'(I8)') NUMEPART

              SELECT CASE (ETYPE)
                  CASE ('tria6')
                      NNOEL_MESH = 6
                  CASE ('quad4')
                      NNOEL_MESH = 4
                  CASE ('hexa8')
                      NNOEL_MESH = 8
                  CASE ('hexa27')
                      NNOEL_MESH = 27
                  CASE DEFAULT
                      WRITE (IOUT,2000)
               END SELECT

                  allocate(incid0(NUME_proc,NNOEL_MESH))
              
                  DO I = numei, numef
                         READ  (IIN2,*) N,(INCID(N,INO),INO=1,NNOEL_MESH)
                         NUME_MESH = NUME_MESH+1
                  ENDDO
                  do iproc = 1, numprocs-1
                      if (iproc .ge. mod(nume,numprocs)) then
                          deallocate(incid0)
                          allocate(incid0(NUME_proc-1,NNOEL_MESH))
                      endif
                      ieli = procsel(iproc)
                      ielf = procsel(iproc+1)-1
                      DO iel = ieli, ielf
                         READ  (IIN2,*) N,(INCID0(iel-ieli+1,INO),INO=1,NNOEL_MESH)
                         NUME_MESH = NUME_MESH+1
                      ENDDO
                      call MPI_Send(incid0, (ielf-ieli+1)*NNOEL_MESH, MPI_INTEGER, iproc, 1, MPI_Comm_world,ierr)
                  enddo

                  deallocate(incid0)

                    NUMPARTS = NUMPARTS + 1

                    IF (NUME_MESH.NE.NUME) THEN
                        WRITE (*,200)
                        WRITE (IERROR,200)
                        read (*,*)
                        STOP
                    ENDIF

                    IF(NUMMAT.NE.NUMPARTS) THEN
                       WRITE (*,300)
                       WRITE (IERROR,300)
                       read (*,*)
                       STOP
                    ENDIF

                ENDDO

              else

                    call MPI_Recv(incid, (numef-numei+1)*NNOEL, MPI_INTEGER, 0, 1, MPI_Comm_world, status,ierr)

              ENDIF
           


1 continue
    
      RETURN
  100 FORMAT (/,/,'***********************************************************************',/,/,'ERROR &
  (INELEM.f90): NNOEL(jobname.cntr) .NE. NNOEL_MESH(jobname.mesh)',/,/,'*************************************')
      
  200 FORMAT (/,/,'***********************************************************************',/,/,'ERROR &
  (INELEM.f90): NUME(jobname.cntr) .NE. NUME_MESH(jobname.mesh)',/,/,'*************************************')
  
  300 FORMAT (/,/,'***********************************************************************',/,/,'ERROR &
  (INELEM.f90): NUMMAT(jobname.cntr) .NE. NUMEPARTS(jobname.mesh)',/,/,'*************************************')
  
 1020 FORMAT (6I8)
 2040 FORMAT (                                                          &
          ' INFORMACOES DOS ELEMENTOS ',//,                             &
          ' ____________________________________________________ ',/)
 1030 FORMAT (6I8)
 
 2000 FORMAT ('TIPO DE ELEMENTO INVALIDO.')
      END
      
     
 !********
