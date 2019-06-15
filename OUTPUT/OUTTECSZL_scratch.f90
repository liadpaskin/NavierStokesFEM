
subroutine OUTTECSZL(x,y,z,disp,up,id,incid,ngl,numnp,nume,nnoel,neq,neqp,nd,istep,time,TFUNC,ireordering,tipo)

      use modtapes,     only: filer, iblank, iplt

    implicit none

    Include "tecio.f90"

    REAL*8 solTime;
    INTEGER debug, i, j, k, dIsDouble, vIsDouble, zoneType, strandID, parentZn, isBlock,ieq,inoel,iel
    INTEGER iCellMax, jCellMax, kCellMax, nFConns, fNMode, shrConn, fileType
    INTEGER nNodes, nCells, nFaces, connectivityCount, index
    INTEGER fileFormat !// 0 == PLT, 1 == SZPLT

    REAL*8, intent(in) :: disp(0:neq),up(0:neqp),x(numnp),y(numnp),z(numnp),TFUNC
    INTEGER, intent(in) :: ngl, numnp, neq, neqp, istep, id(ngl,numnp),nd, incid(nume,nnoel),time,nume,ireordering,nnoel
    real*8, allocatable ::  ap(:,:)
    INTEGER, allocatable ::  connectivity(:)
      character*70 filename

    CHARACTER*1 NULCHAR
    character  tipo
    POINTER     (NullPtr,Null)
    INTEGER*4   Null(*), iblank0

    allocate(connectivity(nume*nnoel))

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
    nCells = nume;
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
    VIsDouble  = 0


    fileFormat = 1; !0: Tecplot binary (.plt)    1: Tecplot subzone (.szplt)
    fileType  = 0;  !0=Full 1=Grid 2=Solution
    dataType = 1    !1: Double-precision (8) 2: Single-precision (4)
                    !3: 32-bit integer (4)   4: 16-bit integer (2)   5: 8-bit unsigned int (1)

     !
     ! Open the file and write the tecplot datafile
     ! header information


      I = tecFileWriterOpen(filename(1:iblank0-1),&
                                'SIMPLE DATASET, ', &
                                'x, y, z, u, v, w', &
                                 FileFormat, &
                                 FileType, &
                                 0, &
                                 fh)

      I = tecDataSetAddAuxData(fh, Name, Value)
      I = tecVarAddAuxData    (fh, vID, Name, Value)
      I = tecZoneAddAuxData   (fh, zID, Name, Value    )

      I = tecZoneCreateFEd(fh,'Zone 01',zoneType, numnp, nume,dataType,'x','y','z','u','v','w', Value Locations,
Passive Variables, Connectivity Sharing Source Zone, Number of Face
Connections, Face Neighbor Mode, Zone Indexf)

I = tecZoneSetUnsteadyOptions (fh, zID, T, istep)
I = tecZoneSetParentZone File Handle, Zone Index, Parent Zone

I = tecZoneVarWriteDoubleValues(fh, zID, 1, 1, numnp,x)
I = tecZoneVarWriteDoubleValues(fh, zID, 2, 1, numnp,y)
I = tecZoneVarWriteDoubleValues(fh, zID, 3, 1, numnp,z)
I = tecZoneVarWriteDoubleValues(fh, zID, 4, 1, numnp,u)
I = tecZoneVarWriteDoubleValues(fh, zID, 5, 1, numnp,v)
I = tecZoneVarWriteDoubleValues(fh, zID, 6, 1, numnp,w)





      ('SIMPLE DATASET, '//NULCHAR, &
                    'x, y, z, u, v, w'//NULCHAR, &
                    filename(1:iblank0-1)//NULCHAR, &
                    '.'//NULCHAR, &
                     FileFormat, &
                     FileType, &

                     Debug, &
                     VIsDouble)



    i = 0
    do iel = 1, nume
        do inoel = 1, nnoel
            i = i + 1
            connectivity(i) = incid( iel , inoel )
        enddo
    enddo
    connectivityCount = nnoel * nume;

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
                  0,      &        !/* TotalNumFaceNodes */
                  0,      &        !/* NumConnectedBoundaryFaces */
                  0,      &        !/* TotalNumBoundaryConnections */
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

    IF (ireordering .eq. 0) then
                write (*,*) 'Error in: Subroutine ensg_vec'
                read (*,*)
                STOP
    else
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
    endif

    deallocate (ap)

else

                write (*,*) 'Error in: Subroutine ensg_vec'
                read (*,*)
                STOP

ENDIF


    i = TECNOD142(connectivity);
    i = TECEND142();

end
