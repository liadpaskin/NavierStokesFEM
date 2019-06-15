subroutine alloc (index)
  
      use modvar  
      use modmesh
      use modloads
      use modsolver
      use modtapes
      use modtim
      
      implicit none
      integer :: index,M_INTEGER,M_REAL,teste1,teste2, iostatus
      REAL*8  :: XVALUE
     
     !*****************************************************************	   
      if (index.eq.1) then
     !*****************************************************************	 
          
          M_REAL = 0
          
          allocate (found  (ngl+1,NUMNP  ));                found(:,:) = .false. 
          allocate (Fno    (NUMNP,NGL+1,2));                Fno(:,:,:) = 0.d0 
          allocate (id        (ngl,numnp  ));                   id(:,:) = 0
          allocate (incid     (nume,nnoel));                incid(:,:) = 0 
          allocate (mtype           (nume));                  mtype(:) = 0 
          allocate (lm           (nume,nd));                   lm(:,:) = 0
          allocate (x              (numnp));                      x(:) = 0.d0 
          allocate (y              (numnp));                      y(:) = 0.d0 
          allocate (z              (numnp));                      z(:) = 0.d0  
          allocate (xi             (numnp));                      x(:) = 0.d0 
          allocate (yi             (numnp));                      y(:) = 0.d0 
          allocate (zi             (numnp));                      z(:) = 0.d0
          allocate (prop   (nummat,ncprop));                 prop(:,:) = 0.d0
          allocate (jdiag          (numnp));                  jdiag(:) = 0 
          allocate (ip           (numnp+1));                     ip(:) = 0 
          allocate (mask           (0:numnp));                   mask(:) = 0
          allocate (xls            (numnp));                    xls(:) = 0
          allocate (invp           (numnp));                   invp(:) = 0
          allocate (lm_tetra     (nume,nd));             lm_tetra(:,:) = 0
          allocate (lm_bar        (nume,6));               lm_bar(:,:) = 0   
          
          M_INTEGER    = (NGL * NUMNP) + (nume * nnoel) + nume + (nume * nd)   
          M_REAL       = (NUMNP * 3)   + (nummat * ncprop)    
          
          CALL MEMORY (M_INTEGER,1,XVALUE)
          VALUE_MEM = XVALUE 
          CALL MEMORY (M_REAL,2,XVALUE) 
          VALUE_MEM = VALUE_MEM + XVALUE 
          
          WRITE(iout,*)            ' *** First Partial: Memory Allocated *** (MEGABYTES) = ',VALUE_MEM                
          WRITE(*,*)            ' *** First Partial: Memory Allocated *** (MEGABYTES) = ',VALUE_MEM                
      
 
     !*****************************************************************	   
      elseif (index.eq.2) then
     !*****************************************************************	  

           M_REAL    = 0.d0
           M_INTEGER = 0.d0
           
                allocate(stiff_Convec       (nd*nd,nume)   );    stiff_Convec       (:,:)  = 0.d0
                allocate(stiff_Difus_Lambda (nd*nd,nume)   );    stiff_Difus_Lambda (:,:)  = 0.d0
                allocate(stiff_Difus_xMu    (nd*nd,nume)   );    stiff_Difus_xMu    (:,:)  = 0.d0
                allocate(stiff              (nd*nd,nume)   );    stiff              (:,:)  = 0.d0
                allocate(Forca(numnp,ngl))                  ;    Forca(:,:)                = 0.d0
                allocate(sig11(nume))                       ;    sig11(:)                  = 0.d0
                allocate(sig12(nume))                       ;    sig12(:)                  = 0.d0
                allocate(sig22(nume))                       ;    sig22(:)                  = 0.d0
        if (etype .eq. 'NavierStokes3D')  then
                allocate(sig13(nume))                       ;    sig13(:)                  = 0.d0
                allocate(sig23(nume))                       ;    sig23(:)                  = 0.d0
                allocate(sig33(nume))                       ;    sig33(:)                  = 0.d0
        endif
                allocate(Fint(nd,nume))                     ;    Fint(:,:)                 = 0.d0
                allocate(elPe(nume))                        ;    elPe(:)                   = 0.d0
                allocate (fp              (0:neq));                 fp   (:)  = 0.d0
                allocate (u               (0:neq));                 u    (:)  = 0.d0
                M_REAL    = M_REAL    + (nd*nd*nume)*4 + (neqp+neq) + (nume*3) + (nume*nd)
           !elseif (isolver_type .eq. 3) then
               ia = nd*nd                                                                      ! 
               ib = nume                                                                       ! 
               kmax1 = kmax + 1                                                                !  
               allocate( ielblk(1))                        ;    ielblk(:)                 = 0  !
               allocate( lm_aux( nume , nd ))              ;    lm_aux(:,:)  = 0.d0            !  ESTRUTURAS DE DADOS GMRES NECESSARIAS
               allocate   ( v                       ( 0 : neq ))  ;     v(:)         = 0.d0    !   PARA OPERACOES SECUNDARIAS
               allocate   ( diag                    ( 0 : neq ))  ;     diag(:)      = 0.d0    ! 
               allocate   ( ve1           ( nume , nd ))  ;     ve1(:,:)     = 0.d0    !
               allocate   ( ve2           ( nume , nd ))  ;     ve2(:,:)     = 0.d0    !
               M_REAL    = M_REAL    + (neq+1)*2 + (nd*nume)*2                                 !
               M_INTEGER = M_INTEGER + (1) + (nume*nd)  
           !endif
           allocate(fpAcel              (0:neq));   fpAcel   (:)  = 0.d0
           call getNALHS(ngl,nnoel,nume,neq,lm,nALHS)
           allocate(Alhs(nALHS)); Alhs(:) = 0.d0
           allocate(Clhs(nALHS)); Clhs(:) = 0.d0
                
           M_REAL = M_REAL + (neq*3) + 3
           
           IF (TIMEF.gt.0.D0) THEN ! Transiente
               allocate (du                (0:neq));                 du     (:)  = 0.d0   
               allocate (upred             (0:neq));               upred    (:)  = 0.d0
               allocate (u_n               (0:neq));                 u_n    (:)  = 0.d0
               allocate (v_n               (0:neq));                 v_n    (:)  = 0.d0
               allocate (fint_pred         (0:neq));     fint_pred    (:)  = 0.d0
               M_REAL = M_REAL + (neq*4) + 4
               !if ( (  etype.eq.'elasticidade3d').or.(etype.eq.'viga3d        ').or.(etype.eq.'tria3_plate        ').or.(etype.eq.'tria3_shell  ') .or.(etype.eq.'ept  ').or.(etype.eq.'epd ').or.(etype.eq.'epdmohrc ') .or. (etype.eq."truss3d") .or.(etype.eq."tdpflex") )  then  
               !   allocate (vpred             (0:neq));               vpred    (:)  = 0.d0
               !   allocate (a_n               (0:neq));                 a_n    (:)  = 0.d0
               !   M_REAL = M_REAL + (neq*2) + 2
               !endif
           endif
           
           CALL MEMORY (M_INTEGER,1,XVALUE)
           VALUE_MEM = XVALUE 
           CALL MEMORY (M_REAL,2,XVALUE) 
           VALUE_MEM = VALUE_MEM + XVALUE 
           
           WRITE(iout,*)            ' *** Second Partial: Memory Allocated *** (MEGABYTES) = ',VALUE_MEM                
           WRITE(*,*)            ' *** Second Partial: Memory Allocated *** (MEGABYTES) = ',VALUE_MEM                
            

     endif
       
return
end
	
	
