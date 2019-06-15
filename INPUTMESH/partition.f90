      SUBROUTINE partition
	
!     LEITURA DOS DADOS DA MALHA E DEFINICAO DE VARIOS ARRANJOS

      use    modvar,    only: numnp,nume,nummat,nnoel,ncprop,etype,ireordering,isolver_type
      use    modmesh,   only: x, y, z, incid, prop, mtype
      use    modtapes,  only: iin1, iin2, iout, ierror
      use modMPIvar

      implicit none

      integer :: numei0, numef0, numnpi0, numnpf0, iproc, i

      allocate (eptr(nume_proc)); eptr(:) = 0
      allocate (eind(nume_proc*nnoel)); eind(:) = 0
      icount = 0
      eptr(1) = 1
      do iel = numei, numef
        do inoel = 1, nnoel
            icount = icount+1
            ino = incid(iel,inoel)
            eptr(icount+1) = eptr(icount)+nnoel
            eind(icount)=ino
      enddo
      
      elmwgt = null
      wgtflag = 0
      numflag = 1
      ncon = 1
      ncommonnodes = 4
      nparts = numprocs
      allocate (tpwgts(ncon*nparts)); tpwgts(:) = (1.0d0/nparts)
      allocate (ubvec(ncon)); ubvec(:) = (1.05d0)
      allocate (part(numnp_proc); part(:) = (0)
      options = 0

      call ParMETIS_V3_PartMeshKway ( &
            procsel, eptr, eind, elmwgt, wgtflag, numflag, &
            ncon, ncommonnodes, nparts, tpwgts, ubvec, &
            options, edgecut, part, MPI_COMM_WORLD)

      deallocate(tpwgts)
      deallocate(ubvec)

      RETURN
      END


