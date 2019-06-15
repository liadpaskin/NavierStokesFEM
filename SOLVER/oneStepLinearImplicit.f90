
subroutine LinearImplicit (Acel, stiff_Mass, stiff_KN, N, T, ngl, nnoel, numnp, neq, neqp, nd, nume, &
u, stiff, lm, up, u_n, uPred, ALHS, CLHS, nALHS, ForceFact, K_x_u,passoutput)

    use modtim, only: dT, Alpha, Nflag, CFLmax_ant, CFLmax, tfunc
    use modloads, only: f, fbody
    use modsolver, only: stiff_difus_lambda, stiff_difus_xmu, stiff_convec, fp, lm_aux , ielblk, ve1, ve2, fpAcel
    use modmesh, only: x,y,z,incid,mtype, prop
    use modvar, only: nwk, nummat, ncprop, niter, etol, etype, nelblk, nsize
    
    IMPLICIT NONE
     
    !      Global     
    INTEGER, INTENT(in)  :: ngl, nnoel, numnp, neq, nd, nume, lm(nume,nd), neqp
    INTEGER, INTENT(inout)  :: nALHS
    REAL*8 , INTENT(in)  :: up(0:neqp)
    REAL*8 , INTENT(out) :: u(0:neq), u_n(0:neq), uPred (0:neq)
    REAL*8 , INTENT(inout) :: ALHS(nALHS), CLHS(nALHS), stiff(nd*nd,nume),Acel(0:neq), stiff_KN(nd*nd,nume),K_x_u (neq)
    
    ! Local
    INTEGER              :: iel, ind, ieq, kk, jj, ii, jnd, N, passoutput, ForceFact, iter, iOutNL
    REAL*8               :: ddot, areatot
    REAL*8               :: Normed_Ku_Conv, Normed_Ku_Conv_Ant
    REAL*8, intent(In)   :: T, stiff_Mass(nd*nd,nume)
    character*70 :: filename
    
    real*8, allocatable  :: sk (:,:), fe(:), ue(:)
    allocate (sk (nd,nd))
    allocate (fe(nd))
    allocate (ue(nd))

    iOutNL = 72

    CFLmax = 1.0d0
    CFLmax_ant = 0.0d0
    !forcefact = 1
    
    FILENAME = './OUT1/LinearImplicitAlgortm.OUT'
    if (N .eq. 1) then
          OPEN  (UNIT=iOutNL,FILE=FILENAME, form = 'formatted',ACTION='write')
          write(iOutNL,*) 'niter, etol:  '  , niter, ' ; ' , etol
          write(iOutNL,*) 'T, niter realizadas, discrepancia(ku_convec) '
    else
         OPEN (UNIT=iOutNL,FILE=FILENAME, position='append', Status= 'unknown',ACTION='write')
    endif
    
    if ((N .eq. 1) .or. (forcefact .eq. 1)) then
        do ind = 1, nd*nd                                                                                                             !
            do iel = 1, nume                                                                                                          !
                stiff_KN (ind, iel) =  (Alpha* dT* (Stiff_difus_Lambda(ind, iel)+Stiff_difus_xMu(ind, iel)))  &
                + stiff_Mass (ind, iel) !  Matriz de rigidez final = Matriz de massa + Matriz de rigidez
            enddo                                                                                                                     !
        enddo
    endif
    
    do ind = 1, neq                                            
        uPred (ind) =  u(ind)  + (1-Alpha) *dT *Acel(ind)      
    enddo          
            
    do ind = 1, neq                                            
        u (ind) =  uPred(ind)      
    enddo     
    
    Normed_Ku_Conv = 0.0d0
    
    do iter = 1, niter
        
        if (etype .eq. 'NavierStokes3D')  then
             if (nnoel == 8) then
                CALL hexa8_NavierStokes3D_Convec (stiff,fp,up,u,lm,x,y,z,incid,mtype,prop,numnp,nume,nummat,nwk,&
                nd,neq+1,neqp,nnoel,ncprop,stiff_Convec,stiff_difus_Lambda, stiff_difus_xMu,areatot)
             elseif (nnoel == 27) then
                CALL hexa27_NavierStokes3D_Convec (stiff,fp,up,u,lm,x,y,z,incid,mtype,prop,numnp,nume,nummat,nwk,&
                nd,neq+1,neqp,nnoel,ncprop,stiff_Convec,stiff_difus_Lambda, stiff_difus_xMu,areatot)
             endif
        elseif (etype .eq. 'NavierStokes2D')  then
             if (nnoel == 6) then
                CALL tria6_NavierStokes2D_Convec (stiff,fp,up,u,lm,x,y,incid,mtype,prop,numnp,nume,nummat,nwk,&
                nd,neq+1,neqp,nnoel,ncprop,stiff_Convec,stiff_difus_Lambda, stiff_difus_xMu,areatot)
             elseif (nnoel == 4) then
                CALL quad4_NavierStokes2D_Convec (stiff,fp,up,u,lm,x,y,incid,mtype,prop,numnp,nume,nummat,nwk,&
                nd,neq+1,neqp,nnoel,ncprop,stiff_Convec,stiff_difus_Lambda, stiff_difus_xMu,areatot)
             endif
        endif 
            
        do ind = 1, neq
            u_n(ind) = f(ind)- fp(ind)*TFUNC - fbody(ind) - fpAcel(ind)
        enddo
    
        do iel = 1, nume                                      
               fe(:) = 0.d0                                   
               ue(:) = 0.d0                                   
               do ind =1, nd                                  
                   ieq = lm(iel,ind)                          
                   if(ieq.lt.0) then                          
                       ue(ind) = 0.d0                         
                   else                                       
                       ue(ind) = u(ieq)                      
                   endif                                      
               enddo                                          
               kk = 0                                         
               do jj = 1 ,nd                                  
                    do ii = 1 ,nd                             
                       kk = kk + 1                            
                       sk(ii,jj) = stiff_convec(kk,iel)       
                    enddo                                     
               enddo                                          
               do ind = 1, nd                                 
                    do jnd = 1, nd                            
                       fe(ind) = fe(ind) + sk(ind,jnd)*ue(jnd)
                    enddo                                     
               enddo                                          
               do ind = 1, nd                                 
                    ieq = lm(iel,ind)                         
                    if (ieq.gt.0) then     
                       u_n(ieq) = u_n(ieq) - fe(ind)          
                    endif                                     
               enddo                                          
        enddo
        
        do ind = 1, neq
            u_n(ind) = u_n(ind) * alpha * dT
        enddo
        
        do iel = 1, nume                                      
               fe(:) = 0.d0                                   
               ue(:) = 0.d0                                   
               do ind =1, nd                                  
                   ieq = lm(iel,ind)                          
                   if(ieq.lt.0) then                          
                       ue(ind) = 0.d0                         
                   else                                       
                       ue(ind) = uPred(ieq)                      
                   endif                                      
               enddo                                          
               kk = 0                                         
               do jj = 1 ,nd                                  
                    do ii = 1 ,nd                             
                       kk = kk + 1                            
                       sk(ii,jj) = stiff_mass(kk,iel)    !   
                    enddo                                     
               enddo                                          
               do ind = 1, nd                                 
                    do jnd = 1, nd                            
                       fe(ind) = fe(ind) + sk(ind,jnd)*ue(jnd)
                    enddo                                     
               enddo                                          
               do ind = 1, nd                                 
                    ieq = lm(iel,ind)                         
                    if (ieq.gt.0) then                        
                       u_n(ieq) = u_n(ieq) + fe(ind)          
                    endif                                     
               enddo                                          
        enddo

        ! BEGIN Solver direto Crout factoration
        call CroutFact_solver(nume, nd, neq, ngl, nnoel, numnp, Lm, stiff_KN, u_n, u, N, T, Nflag, ALHS, CLHS, &
        nALHS, ForceFact)
        ! END   Solver direto Crout factoration
        
        forcefact = 2
        
        K_x_u(:) = 0.d0
        Normed_Ku_Conv_Ant = Normed_Ku_Conv
        CALL matvec ( u, K_x_u, ve1, ve2, neq, ielblk, nsize, nd, nelblk, nume, stiff_Convec, lm_aux )  ! Calcula K(u_n-2)*u_n-1                                                                                                                    !
        Normed_Ku_Conv = ddot(neq, K_x_u, 1, K_x_u, 1)
        
        if (abs(Normed_Ku_Conv_Ant - Normed_Ku_Conv)/Normed_Ku_Conv .lt. etol) then
            exit
        endif
         
    enddo    
    
    do ind = 1, neq
        Acel(ind) = (u(ind) - uPred(ind)) / (Alpha*dT)
    enddo
    
    if (etype .eq. 'NavierStokes3D')  then
             if (nnoel == 8) then
                call hexa8_CFLcalc (N, T, passoutput, dT)
             elseif (nnoel == 27) then
                call hexa27_CFLcalc (N, T, passoutput, dT)
             endif
    elseif (etype .eq. 'NavierStokes2D') then
             if (nnoel == 6) then
                call tria6_CFLcalc (N, T, passoutput, dT)
             elseif (nnoel == 4) then
                call quad4_CFLcalc (N, T, passoutput, dT)
             endif
    endif
    
    write(iOutNL,*) T, iter-1, abs(Normed_Ku_Conv_Ant - Normed_Ku_Conv)/Normed_Ku_Conv
    CLOSE (UNIT=iOutNL)

    deallocate (sk )
    deallocate (fe)
    deallocate (ue)

    return

end
