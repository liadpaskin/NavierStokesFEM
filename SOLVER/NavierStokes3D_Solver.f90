subroutine NavierStokes3D_Solver

    use    modvar
    use    modsolver
    use    modmesh        
    use    modloads
    use    modtapes
    use    modtim 
    
    IMPLICIT NONE
     
    INTEGER              :: inel , ind ,kensight ,i, j, n,ivec, ima, iel, Phi0seted, ForceFact, n0, passoutput
    REAL*8               :: tempo_cpu , tempo_cpu1,tempo_cpu2,tempo_cpu3 ,xtempo,t,func, Normed_Ku, Normed_Ku_Conv, &
    Normed_Ku_Conv_Ant, ddot, test
    REAL*8               :: Normed_Ku_Ant, Normed_Fp, Normed_Ku_Convec, Normed_Ku_Difus, areatot
    real*8, allocatable  :: K_x_u(:), stiff_Mass(:,:), stiff_virtMass(:,:), stiff_KN(:,:), Acel(:), Pressure(:),&
     dAcel(:), dF(:), Lumped_Mass(:)
    
    allocate (stiff_KN  (nd*nd,nume));   stiff_KN  (:,:) = 0.d0
    allocate (stiff_Mass(nd*nd,nume));   stiff_Mass(:,:) = 0.d0
    allocate (Lumped_Mass(0:neq));        Lumped_Mass(:)  = 0.d0
    allocate (stiff_virtMass(nd*nd,nume));   stiff_virtMass(:,:) = 0.d0
    allocate (K_x_u          (0:neq)  );   K_x_u (:) = 0.d0
    allocate (Acel           (0:neq)  );   Acel  (:) = 0.d0
    allocate (dAcel           (0:neq)  );   dAcel  (:) = 0.d0
    allocate (dF              (0:neq)  );   dF     (:) = 0.d0
    allocate (Pressure       (nume) );   Pressure(:) = 0.d0
    
    forcefact = 1
    CFLmax = 0.d0
    T = 0.d0
    Normed_Ku   = 0.d0
    Normed_Ku_Conv  = 0.d0
    Normed_Ku_Conv_Ant  = 0.d0
    N0 = 1
    kensight = 0  
    passoutput = 0

    !call tec_vec_read('U',u,up,id,ngl,numnp,neq,neqp,80)     ! Le condicao inicial
    call initialInterpolate ('U',u,up,id,ngl,numnp,neq,neqp,129)
    !call tec_vec_read('U',uPred,up,id,ngl,numnp,neq,neqp,81) ! Le condicao inicial
    call initialInterpolate ('U',uPred,up,id,ngl,numnp,neq,neqp,130)
    do i = 1, neq                                                !
        Acel(i) = (u(i) - uPred(i)) / dT                         !
    enddo
    TFUNC = FUNC(T)
    call OUTTECPLT3D(x,y,z,u,up,id,incid,ngl,numnp,nume,nnoel,neq,neqp,nd,kensight,T,TFUNC,'U')
    call PlotLine3D(u, up, id, neq, neqp, ngl, kensight)

    do inel=1, nume                                 !
       do ind =1, nd                                !
          if (lm(inel,ind).lt.0) then               !
              lm_aux(inel,ind)=0                    !
        	       else                             !  Prepara a matriz lm em formato para
              lm_aux(inel,ind)=lm(inel,ind)         !    ser lida pelo GMres e seus utilitarios
           endif                                    !
       enddo                                        !
    enddo                                           !
    ielblk(1) = nume                                !
    
        if (etype .eq. 'NavierStokes3D')  then
                  ! Prepara matriz de massa para o campo de velocidades
            if (nnoel == 8) then
                CALL hexa8_NavierStokes3D_Mass (stiff_Mass,x,y,z,incid,mtype,prop,numnp,nume,nummat,nd,nnoel,ncprop)
                      ! Prepara matriz de difusao
                CALL hexa8_NavierStokes3D_Difus (lm,x,y,z,incid,mtype,prop,numnp,nume,nummat,nwk,nd,neq+1,neqp,&
                    nnoel,ncprop,stiff_difus_Lambda, stiff_difus_xMu,areatot)
                      ! Prepara matriz convectiva
                CALL hexa8_NavierStokes3D_Convec (stiff,fp,up,u,lm,x,y,z,incid,mtype,prop,numnp,nume,nummat,nwk,&
                nd,neq+1,neqp,nnoel,ncprop,stiff_Convec,stiff_difus_Lambda, stiff_difus_xMu,areatot)
            elseif (nnoel == 27) then
                CALL hexa27_NavierStokes3D_Mass (stiff_Mass,x,y,z,incid,mtype,prop,numnp,nume,nummat,nd,nnoel,ncprop)
                      ! Prepara matriz de difusao
                CALL hexa27_NavierStokes3D_Difus (lm,x,y,z,incid,mtype,prop,numnp,nume,nummat,nwk,nd,neq+1,neqp,&
                    nnoel,ncprop,stiff_difus_Lambda, stiff_difus_xMu,areatot)
                      ! Prepara matriz convectiva
                CALL hexa27_NavierStokes3D_Convec (stiff,fp,up,u,lm,x,y,z,incid,mtype,prop,numnp,nume,nummat,nwk,&
                nd,neq+1,neqp,nnoel,ncprop,stiff_Convec,stiff_difus_Lambda, stiff_difus_xMu,areatot)
            endif
        elseif (etype .eq. 'NavierStokes2D')  then
            if (nnoel == 6) then
                CALL tria6_NavierStokes2D_Mass (stiff_Mass,x,y,incid,mtype,prop,numnp,nume,nummat,nd,nnoel,ncprop)
                CALL tria6_NavierStokes2D_Difus (lm,x,y,incid,mtype,prop,numnp,nume,nummat,nwk,nd,neq+1,neqp,&
                    nnoel,ncprop,stiff_difus_Lambda, stiff_difus_xMu,areatot)
                CALL tria6_NavierStokes2D_Convec (stiff,fp,up,u,lm,x,y,incid,mtype,prop,numnp,nume,nummat,nwk,&
                nd,neq+1,neqp,nnoel,ncprop,stiff_Convec,stiff_difus_Lambda, stiff_difus_xMu,areatot)
             elseif (nnoel == 4) then
                CALL quad4_NavierStokes2D_Mass (stiff_Mass,x,y,incid,mtype,prop,numnp,nume,nummat,nd,nnoel,ncprop)
                CALL quad4_NavierStokes2D_Difus (lm,x,y,incid,mtype,prop,numnp,nume,nummat,nwk,nd,neq+1,neqp,&
                    nnoel,ncprop,stiff_difus_Lambda, stiff_difus_xMu,areatot)
                CALL quad4_NavierStokes2D_Convec (stiff,fp,up,u,lm,x,y,incid,mtype,prop,numnp,nume,nummat,nwk,&
                nd,neq+1,neqp,nnoel,ncprop,stiff_Convec,stiff_difus_Lambda, stiff_difus_xMu,areatot)
             endif
        endif
                  
    OPEN  (UNIT=22,FILE='./OUT1/factored.OUT', form = 'formatted')  
    write (22,'(a)') 'Tempo, CFL maximo'
    CLOSE (unit=22)
    
    DO N = N0, NSTEPS
                
            T = T + DT
            TFUNC = FUNC(T)
            !if (TFUNC .ne. FUNC(T)) then
            !    call acelbound(t,stiff_mass)
            !endif

            if (isolver_type .eq. 4) then
                call LinearImplicit (Acel, stiff_Mass, stiff_KN, N, T, ngl, nnoel, numnp, neq, neqp, &
                nd, nume, u, stiff, lm, up, u_n, uPred, ALHS, CLHS, nALHS, ForceFact, K_x_u,passoutput)
            else
                write(*,*) 'SOLVER NOT IMPLEMENTED'
                read (*,*)
                STOP
            endif

             !Chama SUBROTINA DE SAIDA
            if (passoutput .ne. 1) then
                Call NavierStokes3D_OUT (kensight, K_x_u, t, N, Normed_Ku, Normed_Ku_Ant, Normed_Ku_Conv,&
                 Normed_Ku_Conv_Ant, Normed_Ku_Difus, Normed_Ku_Convec, Normed_Fp, Pressure)
            endif
                                                                                                                  !
    ENDDO  ! fim loop tempo   
    
    CALL SYSTEM_CLOCK(it2, count_rate = clockcount2) 
    write (IOUT, *) 'Tempo de CPU (s) : ', (it2/clockcount1 - it1/clockcount2)
    

    deallocate (stiff_KN )
    deallocate (stiff_Mass)
    deallocate (Lumped_Mass)
    deallocate (stiff_virtMass)
    deallocate (K_x_u)
    deallocate (Acel )
    deallocate (dAcel)
    deallocate (dF )
    deallocate (Pressure)
    deallocate(Alhs)
    deallocate(Clhs)

RETURN
END

subroutine NavierStokes3D_OUT (kensight, K_x_u, t, N, Normed_Ku, Normed_Ku_Ant, Normed_Ku_Conv, &
Normed_Ku_Conv_Ant, Normed_Ku_Difus, Normed_Ku_Convec, Normed_Fp, Pressure)

    use    modvar
    use    modsolver
    use    modmesh
    use    modloads
    use    modtapes
    use    modtim

IMPLICIT NONE

INTEGER      :: kensight, N
REAL*8       :: Normed_Ku, Normed_Ku_Ant, Normed_Ku_Conv, Normed_Ku_Conv_Ant, Normed_Fp, ddot, &
Normed_Ku_Difus, Normed_Ku_Convec, t, K_x_u(neq), Pressure(nume)
CHARACTER*70 :: filename1, filename2, filename3, filename4

            !FILENAME1 = './OUT1/PICARDcntr.OUT'
            FILENAME2 = './OUT1/RESIDUALcntr.OUT'
            FILENAME3 = './OUT2/time.out'
            FILENAME4 = './OUT1/CPUtime.out'
            if (N .eq. 1) then
                  !OPEN  (UNIT=22,FILE=FILENAME1, form = 'formatted')
                  OPEN  (UNIT=23,FILE=FILENAME2, form = 'formatted')
                  OPEN  (UNIT=24,FILE=FILENAME3, form = 'formatted')
                  OPEN  (UNIT=25,FILE=FILENAME4, form = 'formatted')
                  write(23,'(a)') 'T, ABS(Ku - Ku_Ant)/Ku, ABS(Ku_Conv - Ku_Conv_Ant)/Ku_Conv, Ku_difus/Ku, &
                  Ku_convec/Ku '
                  write(25,'(a)') 'T, CPU time '
            else
                 !OPEN (UNIT=22,FILE=FILENAME1, position='append', Status= 'unknown')
                 OPEN (UNIT=23,FILE=FILENAME2, position='append', Status= 'unknown')
                 OPEN (UNIT=24,FILE=FILENAME3, position='append', Status= 'unknown')
                 OPEN (UNIT=25,FILE=FILENAME4, position='append', Status= 'unknown')
            endif

            IF (MOD(N-1,nflag) .EQ. 0)  THEN                                 !

                    K_x_u(:) = 0.d0                                                                                 !
                    CALL matvec ( u, K_x_u, ve1, ve2, neq, ielblk, nsize, nd, nelblk, nume, stiff, lm_aux )         ! Calcula K(u_n-1)*u_n-1                                                                                                               !
                    Normed_Ku = ddot(neq, K_x_u, 1, K_x_u, 1)
                    K_x_u(:) = 0.d0                                                                          !
                    CALL matvec ( u, K_x_u, ve1, ve2, neq, ielblk, nsize, nd, nelblk, nume, stiff_Convec, lm_aux )  ! Calcula K(u_n-2)*u_n-1                                                                                                                    !
                    Normed_Ku_Conv = ddot(neq, K_x_u, 1, K_x_u, 1)                                                  !                                                                                                                                      !
                    K_x_u(:) = 0.d0
                    CALL addmatvec ( u, K_x_u, ve1, ve2, neq, ielblk, nsize, nd, nelblk, nume,  stiff_difus_Lambda,&
                         stiff_difus_xMu, lm_aux )                              !
                    Normed_Ku_Difus  = ddot(neq, K_x_u, 1, K_x_u, 1)                                                                                     !                                                                                          !  OUTPUT
                    K_x_u(:) = 0.d0                                                                                                                                      !
                    CALL matvec ( u, K_x_u, ve1, ve2, neq, ielblk, nsize, nd, nelblk, nume, stiff_convec, lm_aux )                                                       !   em tela
                    Normed_Ku_Convec        = ddot(neq, K_x_u, 1, K_x_u, 1)
                    Normed_Fp = ddot(neq, fp, 1, fp, 1)

                  WRITE(*,*)                                                 !
                  WRITE(*,*)'PASSO DE TEMPO'                                 !  Output formato
                  WRITE(*,*)N,NSTEPS                                         !     Ensight
                  WRITE(*,*)                                                 !
                  WRITE(*,*)'TEMPO DA AN�LISE'                               !
                  WRITE(*,*)T,TIMEF                                          !
                  WRITE(*,*)                                                 !
                  kensight = kensight +1                                     !

                    if (etype .eq. 'NavierStokes3D')  then
                        if (nnoel == 8) then
                            call Calc_Pressure_Hexa8(Pressure, TFUNC,N,T)
                        elseif (nnoel == 27) then
                            call Calc_Pressure_Hexa27(Pressure, TFUNC,N,T)
                        endif
                        call OUTTECPLT3D(x,y,z,u,up,id,incid,ngl,numnp,nume,nnoel,neq,neqp,nd,kensight,T,TFUNC,'U')
                        call PlotLine3D(u, up, id, neq, neqp, ngl, kensight)
                    elseif (etype .eq. 'NavierStokes2D')  then
                        if (nnoel == 6) then
                            call Calc_Pressure_Tria6(Pressure, TFUNC,N,T)
                        elseif (nnoel == 4) then
                            call Calc_Pressure_Quad4(Pressure, TFUNC,N,T)
                        endif
                        call OUTTECPLT2D(x,y,u,up,id,incid,ngl,numnp,nume,nnoel,neq,neqp,nd,kensight,T,TFUNC,'U')
                        call PlotLine2D(u, up, id, neq, neqp, ngl, kensight)
                    endif
                  !call ensg_vec ('U',u,up,id,ngl,numnp,neq,neqp,kensight)    !
                  ! call ensg_vec ('F',forca,up,id,ngl,numnp,neq,neqp,kensight)    !
                  ! call ensg_tensor2_elem ('T',sig11,sig12,sig22,nume,etype,kensight)
                  !call ensg_scalar_elem ('Pre',Pressure,nume,etype,kensight)
                  !call ensg_scalar_elem ('Pec',elPe    ,nume,etype,kensight)

                  write(24,*) T
                                                                                                                                 !
                  write(*  ,*)   ' '
                  write(23  ,'(5e12.4E3)') T, ABS(Normed_Ku - Normed_Ku_Ant)/ Normed_Ku, &
                  ABS(Normed_Ku_Conv - Normed_Ku_Conv_Ant)/ Normed_Ku_Conv, Normed_Ku_Difus  &
                   / Normed_Ku  , Normed_Ku_Convec / Normed_Ku
                  write(*   ,'(a,1e12.4E3,a)') ' [(Ku-Ku_ant)/Ku]_Conv = ' , &
                  100 * ABS(Normed_Ku_Conv - Normed_Ku_Conv_Ant)/ Normed_Ku_Conv, ' % '
                  write(*   ,'(a,1e12.4E3,a)') ' (Ku-Ku_ant)/Ku = ' , &
                  100 * ABS(Normed_Ku - Normed_Ku_Ant)/ Normed_Ku, ' % '
                  write(*   ,'(a,1e12.4E3,a)') '  Ku_Difus  /Ku = ' , &
                  100 *               Normed_Ku_Difus / Normed_Ku, ' % '
                  write(*   ,'(a,1e12.4E3,a)') '  Ku_Convec /Ku = ' , &
                  100 *              Normed_Ku_Convec / Normed_Ku, ' % '

                  CALL SYSTEM_CLOCK(it2, count_rate = clockcount2)
                  write (25, *) T, (it2/clockcount1 - it1/clockcount2)

                  !CLOSE (UNIT=22)
                  CLOSE (UNIT=23)
                  CLOSE (UNIT=24)
                  CLOSE (UNIT=25)
            endif

            return
    end
