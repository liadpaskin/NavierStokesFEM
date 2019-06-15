
subroutine OUTTECPLT2D(x,y,disp,up,id,incid,ngl,numnp,nume,nnoel,neq,neqp,nd,istep,time,TFUNC,tipo)

      use modtapes,     only: filer, iblank, iplt

    implicit none

    Include "tecio.f90"

    REAL*8 solTime;
    INTEGER*4 debug, i, j, k, dIsDouble, vIsDouble, zoneType, strandID, parentZn, isBlock,ieq,inoel,iel,zero
    INTEGER*4 iCellMax, jCellMax, kCellMax, nFConns, fNMode, shrConn, fileType
    INTEGER*4 nNodes, nCells, nFaces, connectivityCount, index
    INTEGER*4 fileFormat !// 0 == PLT, 1 == SZPLT

    REAL*8, intent(in) :: disp(0:neq),up(0:neqp),x(numnp),y(numnp),TFUNC
    INTEGER, intent(in) :: ngl, numnp, neq, neqp, istep, id(ngl,numnp),nd, incid(nume,nnoel),time,nume,nnoel
    real*8, allocatable ::  ap(:,:),x0(:),y0(:),u0(:),v0(:)
    INTEGER*4, allocatable ::  connectivity(:),connectivity0(:),map(:)
      character*70 filename

    CHARACTER*1 NULCHAR
    character  tipo
    POINTER     (NullPtr,Null)
    INTEGER*4   Null(*), iblank0, marc, cont,ino

    allocate(connectivity(nume*4))
    allocate(map(numnp))

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
    zoneType  = 3;      !/* Quadrilateral */
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
    zero = 0;

     !
     ! Open the file and write the tecplot datafile
     ! header information


      I = TECINI142('SIMPLE DATASET, '//NULCHAR, &
                    'x, y, u, v'//NULCHAR, &
                    filename(1:iblank0-1)//NULCHAR, &
                    '.'//NULCHAR, &
                     FileFormat, &
                     FileType, &
                     Debug, &
                     VIsDouble)

    i = 0
    do iel = 1, nume
        do inoel = 1, 3
            i = i + 1
            connectivity(i) = incid( iel , inoel )
        enddo
        i = i + 1
        if (nnoel == 6) then
            connectivity(i) = incid( iel , 3 )
        elseif (nnoel == 4) then
            connectivity(i) = incid( iel , 4 )
        endif
    enddo
    connectivityCount = nnoel * nume;

    allocate (ap(numnp,ngl)) ; ap = 0.d0

    cont = 0
    do i = 1, numnp
       if (nnoel == 6) then
           marc = 0
           do iel = 1, nume
                if (i == incid(iel,4) .or. i == incid(iel,5) .or. i == incid(iel,6)) then
                    marc = 1
                    exit
                endif
           enddo
           if (marc == 1) then
                cycle
           endif
           cont =  cont + 1
           map(cont) = i
       endif
       do j = 1, ngl
       ieq = id(j,i)
       if (ieq.ge.0) then
         ap(i,j)=disp(ieq)
       else
         ap(i,j)=up(-ieq)*TFUNC
       endif
     enddo
    enddo

if (nnoel == 6) then
    nNodes = cont
    allocate(x0(cont))
    allocate(y0(cont))
    allocate(u0(cont))
    allocate(v0(cont))
    allocate(connectivity0(nume*4))
    connectivity0(:)=connectivity(:)
    do i = 1, cont
        ino = map(i)
        x0(i) = x(ino)
        y0(i) = y(ino)
        u0(i) = ap(ino,1)
        v0(i) = ap(ino,2)
        do j = 1,(nume*4)
            if (connectivity(j) == ino) then
                connectivity0(j) = i
            endif
        enddo
    enddo
endif

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

if (nnoel == 6) then
                i = TECDAT142(nNodes, x0, dIsDouble);
                i = TECDAT142(nNodes, y0, dIsDouble);
                i = TECDAT142(nNodes, u0, dIsDouble);
                i = TECDAT142(nNodes, v0, dIsDouble);
                i = TECNOD142(connectivity0);
else
                i = TECDAT142(nNodes, x, dIsDouble);
                i = TECDAT142(nNodes, y, dIsDouble);
                i = TECDAT142(nNodes, ap(:,1), dIsDouble);
                i = TECDAT142(nNodes, ap(:,2), dIsDouble);
                i = TECNOD142(connectivity);
endif

    deallocate (ap)

    i = TECEND142();

end
