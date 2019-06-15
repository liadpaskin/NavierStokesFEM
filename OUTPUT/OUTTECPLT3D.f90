
subroutine OUTTECPLT3D(x,y,z,disp,up,id,incid,ngl,numnp,nume,nnoel,neq,neqp,nd,istep,time,TFUNC,tipo)

      use modtapes,     only: filer, iblank, iplt

    implicit none

    Include "tecio.f90"

    REAL*8 solTime
    INTEGER*4 debug, i, j, k, dIsDouble, vIsDouble, zoneType, strandID, parentZn, isBlock,ieq,inoel,iel, zero
    INTEGER*4 iCellMax, jCellMax, kCellMax, nFConns, fNMode, shrConn, fileType,nnd(3,3,3)
    INTEGER*4 nNodes, nCells, nFaces, connectivityCount, index
    INTEGER*4 fileFormat !// 0 == PLT, 1 == SZPLT

    REAL*8, intent(in) :: disp(0:neq),up(0:neqp),x(numnp),y(numnp),z(numnp),TFUNC
    INTEGER, intent(in) :: ngl, numnp, neq, neqp, istep, id(ngl,numnp),nd, incid(nume,nnoel),time,nume,nnoel
    real*8, allocatable ::  ap(:,:)
    INTEGER*4, allocatable ::  connectivity(:)
      character*70 filename

    CHARACTER*1 NULCHAR
    character  tipo
    POINTER     (NullPtr,Null)
    INTEGER*4   Null(*), iblank0,iel0,ielx,iely,ielz

      filename = filer

!     arquivo Ensight para vetores
      write(*,*)'stop',istep
      write (filename(1:4),'(i4.4)') istep

      write(filename,'(a)') './OUT2/'//filer(1:iblank-1)//'.'//tipo//filename(1:4)

    DO I = 1,80
      IF (filename(I:I).EQ.' ') THEN
        IBLANK0=I
        GO TO 10
      ENDIF
    ENDDO
10    CONTINUE

    NULCHAR    = CHAR(0)
    NullPtr    = 0
    fileFormat = 0;
    debug     = 1;
    vIsDouble = 1;
    dIsDouble = 1;
    nNodes = numnp;
    zoneType  = 5;      !/* Brick */
    solTime   = time;
    strandID  = 0;      !/* StaticZone */
    parentZn  = 0;      !/* No Parent */
    isBlock   = 1;      !/* Block */
    iCellMax  = 0;
    jCellMax  = 0;
    kCellMax  = 0;
    nFConns   = 0;
    fNMode    = 0;
    shrConn   = 0;
    fileType  = 0;      !0=Full 1=Grid 2=Solution
    VIsDouble  = 1;
    zero = 0

    if (nnoel == 8) then
        nCells = nume;
    elseif (nnoel == 27) then
        nCells = nume*8
    endif
     !
     ! Open the file and write the tecplot datafile
     ! header information


      I = TECINI142('SIMPLE DATASET, '//NULCHAR, &
                    'x, y, z, u, v, w'//NULCHAR, &
                    filename(1:iblank0-1)//NULCHAR, &
                    '.'//NULCHAR, &
                     FileFormat, &
                     FileType, &
                     Debug, &
                     VIsDouble)

    if (nnoel == 8) then
         allocate(connectivity(nume*nnoel))
        i = 0
        do iel = 1, nume
            do inoel = 1, nnoel
                i = i + 1
                connectivity(i) = incid( iel , inoel )
            enddo
        enddo
        connectivityCount = nnoel * nume;
    elseif (nnoel == 27) then

        allocate(connectivity(8*nCells))
        nnd(1, 1, 1) = 1;
        nnd(3, 1, 1) = 2;
        nnd(3, 3, 1) = 3;
        nnd(1, 3, 1) = 4;
        nnd(1, 1, 3) = 5;
        nnd(3, 1, 3) = 6;
        nnd(3, 3, 3) = 7;
        nnd(1, 3, 3) = 8;
        nnd(1, 1, 2) = 9;
        nnd(3, 1, 2) = 10;
        nnd(3, 3, 2) = 11;
        nnd(1, 3, 2) = 12;
        nnd(2, 2, 1) = 13;
        nnd(2, 2, 2) = 14;
        nnd(2, 2, 3) = 15;
        nnd(2, 1, 1) = 16;
        nnd(3, 2, 1) = 17;
        nnd(2, 3, 1) = 18;
        nnd(1, 2, 1) = 19;
        nnd(2, 1, 2) = 20;
        nnd(3, 2, 2) = 21;
        nnd(2, 3, 2) = 22;
        nnd(1, 2, 2) = 23;
        nnd(2, 1, 3) = 24;
        nnd(3, 2, 3) = 25;
        nnd(2, 3, 3) = 26;
        nnd(1, 2, 3) = 27;

        i = 0
        do iel0 = 1, nume
            do ielx = 1, 2
                do iely = 1, 2
                    do ielz = 1, 2
                        i = i + 1
                        connectivity(i) = incid( iel0 , nnd(ielx, iely, ielz) )
                        i = i + 1
                        connectivity(i) = incid( iel0 , nnd(ielx+1, iely, ielz) )
                        i = i + 1
                        connectivity(i) = incid( iel0 , nnd(ielx+1, iely+1, ielz) )
                        i = i + 1
                        connectivity(i) = incid( iel0 , nnd(ielx, iely+1, ielz) )
                        i = i + 1
                        connectivity(i) = incid( iel0 , nnd(ielx, iely, ielz+1) )
                        i = i + 1
                        connectivity(i) = incid( iel0 , nnd(ielx+1, iely, ielz+1) )
                        i = i + 1
                        connectivity(i) = incid( iel0 , nnd(ielx+1, iely+1, ielz+1) )
                        i = i + 1
                        connectivity(i) = incid( iel0 , nnd(ielx, iely+1, ielz+1) )
                    enddo
                 enddo
            enddo
        enddo
        connectivityCount = 8 * nCells;
    endif

    !connectivityCount = 8 * nCells;
    !connectivity = (INTEGER4*)malloc(connectivityCount * sizeof(INTEGER4));

     !* Write the zone header information.

    i = TECZNE142('Simple Zone',&
                  zoneType,&
                  nNodes,&
                  nCells,&
                  nFaces,&
                  iCellMax,&
                  jCellMax,&
                  kCellMax,&
                  solTime,&
                  strandID,&
                  parentZn,&
                  isBlock,&
                  nFConns,&
                  fNMode,&
                  zero,      &        !/* TotalNumFaceNodes */
                  zero,      &        !/* NumConnectedBoundaryFaces */
                  zero,      &        !/* TotalNumBoundaryConnections */
                  NULL,   &        !/* PassiveVarList */
                  NULL,   &        !/* ValueLocation = Nodal */
                  NULL,   &        !/* SharVarFromZone */
                  &shrConn)
     !
     ! Write out the field data.
     !

if (neqp.gt.0) then

    allocate (ap(numnp,ngl)) ; ap = 0.d0

    do i = 1, numnp
       do j = 1, ngl
       ieq = id(j,i)
       if (ieq.ge.0) then
         ap(i,j)=disp(ieq)
       else
         ap(i,j)=up(-ieq)*TFUNC
       endif
     enddo
    enddo

              IF (ngl .eq. 3) THEN
                i = TECDAT142(nNodes, x, dIsDouble);
                i = TECDAT142(nNodes, y, dIsDouble);
                i = TECDAT142(nNodes, z, dIsDouble);
                i = TECDAT142(nNodes, ap(:,1), dIsDouble);
                i = TECDAT142(nNodes, ap(:,2), dIsDouble);
                i = TECDAT142(nNodes, ap(:,3), dIsDouble);
              ELSE
                write (*,*) 'Error in: Subroutine ensg_vec'
                read (*,*)
                STOP
              ENDIF

    deallocate (ap)

else

                write (*,*) 'Error in: Subroutine ensg_vec'
                read (*,*)
                STOP

ENDIF


    i = TECNOD142(connectivity);
    i = TECEND142();

end
