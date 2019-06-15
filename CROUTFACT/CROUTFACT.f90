function coldot(a,b,n) 
 
   !  program to compute the dot product of vectors stored column-wise 
 
      implicit real*8 (a-h,o-z) 

      dimension a(n),b(n) 
      common /consts/ zero,pt1667,pt25,pt5,one,two,three,four,five,six 
      
      coldot = 0.d0 
      do 100 i=1,n 
            coldot = coldot + a(i)*b(i) 
  100 continue 
      return 
end
    
subroutine addlhs(alhs,clhs,eleffm,idiag,lm,nee,ldiag,lsym,nALHS)  ! Transfere os termos da matriz de rigidez do elemento para a matriz de rigidez global

             ! IN:
         !  nee:     nen*ned = nnoel*ngl = nd
         !  Lm:      entra como um pointer para onde come�a a matriz lm do elemento da vez
         !  eleffm:  Matriz de rigidez do elemento
         !  idiag:   ???
         !  (...)
         !  ldiag = .true.   -> matriz diagonal                
         !  ldiag  = .false. -> matriz n�o diagonal
             !lsym = .true.    -> matriz sim�trica
             !lsym = .false.   -> matriz n�o sim�trica
           
             ! IN / OUT:
         !  alhs:    ???
         !  clhs:    ???
             
      implicit real*8 (a-h,o-z) 
 
      integer*8 :: nALHS

      logical ldiag,lsym 
      dimension idiag(*),lm(*) 
      
      real*8 :: alhs(nALHS), clhs(nALHS), eleffm(nee,nee)
 
      if (ldiag) then ! Matriz diagonal
         do 100 j=1,nee 
            k = lm(j) 
            if (k.gt.0) then 
               l = idiag(k)
               alhs(l) = alhs(l) + eleffm(j,j) 
            endif 
 100     continue 
 
      else ! Matriz n�o diagonal
         do 400 j=1,nee 
             k = lm(j) 
             if (k.gt.0) then 
                do 200 i=1,j ! Banda superior da matriz local
                      m = lm(i) 
                      if (m.gt.0) then 
                         if (k.ge.m) then ! Banda superior da matriz global
                            l = idiag(k) - k + m 
                         else  ! Banda inferior da matriz global
                            l = idiag(m) - m + k 
                         endif 
                         alhs(l) = alhs(l) + eleffm(i,j) 
                      endif 
 200           continue 
             
                if (.not. lsym) then ! Matriz nao simetrica
                   do 300 i = j,nee 
                      m = lm(i) 
                      if (m .gt. 0) then 
                         if (k .ge. m) then 
                           l = idiag(k) - k + m 
                         else 
                           l = idiag(m) - m + k 
                         endif 
                      clhs(l) = clhs(l) + eleffm(i,j) 
                      endif 
 300               continue 
                endif 
             endif 
 400    continue 
 
      endif 
 
      return 
end 

subroutine diag(idiag,neq,n) 

       ! program to compute diagonal addresses of left-hand-side matrix 
       
       Implicit none
       ! Globais
      INTEGER, intent(In)     :: neq
      INTEGER*8, Intent(Out)    :: n
      INTEGER*8, intent(Inout) :: idiag(neq)
       ! Locais
      INTEGER :: ieq, test

      n = 1 
      idiag(1) = 1 
      if (neq.eq.1) return 
      
      do 100 ieq=2,neq 
           !if (ieq.eq.897661) pause
           !test = huge(idiag(1))
           idiag(ieq) = idiag(ieq) + idiag(ieq-1) + 1 
  100 continue 
      n = idiag(neq) 

      return 
end 
    
subroutine colht(idiag,lm,ned,nen,numel,neq) ! routine to compute column heights in global left-hand-side matrix 
 
         ! IN:
     !  lm (ned,nen,numel) -> retorna o numero da equa��o de determinado grau de liberdade
     !  nen   = nnoel (Numero de nos por elemento)
     !  ned   = ngl (Numero de graus de liberdade do n�)
     !  numel = nume (Numero de elementos na malha) 
     !  neq   = neq (Numero de equa��es no problema = numnp * ngl)
         ! OUT:
     !  idiag
     
       ! Globais
     INTEGER, intent (In)  :: lm(ned,nen,numel) , ned, nen, numel, neq
     INTEGER*8, intent (Out) :: idiag(neq)
 
     do 500 k=1,numel 
         
         min = neq 
         do 200 j=1,nen 
             do 100 i=1,ned 
                  num = lm(i,j,k) 
                  if (num.gt.0) min = min0(min,num) 
             100 continue 
         200 continue 
     
         do 400 j=1,nen 
              do 300 i=1,ned 
                  num = lm(i,j,k) 
                  if (num.gt.0) then 
                     m = num - min 
                     if (m.gt.idiag(num)) idiag(num) = m 
                  endif 
              300 continue 
         400 continue 
 
 500 continue 
 
      return 
end 
    
subroutine factns(a,c,idiag,neq) 
              ! 
              !.... program to perform crout factorization: a = l * d * u 
              ! 
              !        a(i):  coefficient matrix stored in compacted column form; 
              !               after factorization contains d and u 
              ! 
              !        c(i):  non-symmetric lower triangular coefficient matrix stored in 
              !                compacted row form; after factorization contains l 
              ! 
              ! 
      implicit real*8 (a-h,o-z) 
       ! deactivate above card(s) for single-precision operation 
       
      dimension a(*),c(*),idiag(*) 
      common /consts/ zero,pt1667,pt25,pt5,one,two,three,four,five,six 
      
      jj = 0 
      do 300 j=1,neq 
           jjlast = jj 
           jj     = idiag(j) 
           jcolht = jj - jjlast 
           if (jcolht.gt.2) then 
               ! for column j and i.le.j-1, replace a(i,j) with d(i,i)*u(i,j) 
              istart = j - jcolht + 2 
              jm1    = j - 1 
              ij     = jjlast + 2 
              ii     = idiag(istart-1) 
              do 100 i=istart,jm1 
                  iilast = ii 
                  ii     = idiag(i) 
                  icolht = ii - iilast 
                  length = min0(icolht-1,i - istart + 1) 
                  if (length.gt.0)  then 
                     a(ij) = a(ij) - coldot(a(ij-length),c(ii-length),length) 
                     c(ij) = c(ij) - coldot(c(ij-length),a(ii-length),length) 
                  endif 
                  ij = ij + 1 
  100         continue 
           endif 
           
           if (jcolht.ge.2) then ! for column j and i.le.j-1, replace a(i,j) with u(i,j); replace a(j,j) with d(j,j). 
              jtemp = j - jj 
              do 200 ij=jjlast+1,jj-1 
                  ii = idiag(jtemp + ij) 
                  
                 ! warning: the following calculations are skipped if a(ii) equals zero 
                  if (a(ii).ne.0.d0) then
                      c(ij) = c(ij)/a(ii) 
                      a(jj) = a(jj) - c(ij)*a(ij) 
                      a(ij) = a(ij)/a(ii) 
                  endif 
  200         continue 
           endif 
  300  continue 
       return 
end        
           
subroutine backns(a,c,b,idiag,neq) ! program to perform forward reduction and back substitution 
           
           implicit real*8 (a-h,o-z) 
           
           dimension a(*),c(*),b(*),idiag(*) 
           common /consts/ zero,pt1667,pt25,pt5,one,two,three,four,five,six 
             
             ! forward reduction 
           jj = 0 
           do 100 j=1,neq 
              jjlast = jj 
              jj     = idiag(j) 
              jcolht = jj - jjlast 
              if (jcolht.gt.1) then 
                  b(j) = b(j) - coldot(c(jjlast+1),b(j-jcolht+1),jcolht-1) 
              endif 
100        continue 
           
             ! diagonal scaling 
           do 200 j=1,neq 
               ajj = a(idiag(j)) 
                 ! warning: diagonal scaling is not performed if ajj equals zero 
               if (ajj.ne.zero) b(j) = b(j)/ajj 
200        continue 
           
             ! back substitution 
           if (neq.eq.1) return 
           jjnext = idiag(neq) 
           do 400 j=neq,2,-1 
               jj     = jjnext 
               jjnext = idiag(j-1) 
               jcolht = jj - jjnext 
               if (jcolht.gt.1) then 
                  bj = b(j) 
                  istart = j - jcolht + 1 
                  jtemp  = jjnext - istart + 1 
                  do 300 i=istart,j-1 
                      b(i) = b(i) - a(jtemp+i)*bj 
  300        continue 
          endif 
  400 continue 
      return 
    end       
    
subroutine CroutFact_solver(nume, nd, neq, ngl, nnoel, numnp, Lm, stiff, u_n, u, N, T, Nflag, ALHS, CLHS, nALHS, ForceFact)

    use modtim, only: CFLmax, CFLmax_ant

    implicit none
    
    !      Global  
    INTEGER, INTENT(in)  :: nume, nd, neq, ngl, numnp, nnoel, Lm(nume,nd)
    REAL*8 , INTENT(in)  :: stiff (nd*nd,nume), u_n(0:neq), T
    REAL*8 , INTENT(inout)  :: u(0:neq), ALHS(nALHS), CLHS(nALHS)
    
    !      Local
    INTEGER*8              :: nALHS
    INTEGER              :: iel, jnd, iCol, iLin, ieq, innoel, ingl, N, Nflag, ForceFact
    INTEGER*8, allocatable :: lm_Aux(:,:), lm_Parab(:,:,:), idiag(:)
    REAL*8 , allocatable :: Brhs(:), eleFFM(:,:)
    
    allocate(Lm_aux(nd,Nume))          ; Lm_aux  (:,  :) = 0
    allocate(Lm_Parab(ngl,nnoel,Nume)) ; Lm_Parab(:,:,:) = 0
    
    allocate(idiag(neq))        ; idiag(:)      = 0
    allocate(Brhs(neq))         ; Brhs(:)       = 0.d0
    allocate(eleFFM(nd,nd))     ; eleFFM(:,:)   = 0.d0
    
    do iel=1,nume                         !
       do jnd=1, nd                       !
           lm_aux(jnd,iel) = lm(iel,jnd)  ! Transforma lm(nume,nd) -> lm(nd,nume)
       enddo                              !                         = lm (ned,nen,numel)
    enddo                                 !
    do iel=1,nume                                                           !
       do innoel=1, nnoel                                                   !
           do ingl = 1, ngl                                                 !
                lm_parab(ingl,innoel,iel) = lm_aux((innoel-1)*ngl+ingl,iel) ! Transforma lm(nd,nume) -> lm (ned,nen,numel)
           enddo                                                            !       
       enddo                                                                !
    enddo                                                                   !
    
    call colht(idiag,lm_Parab,ngl,nnoel,nume,neq) ! Forma o idiag
    call diag(idiag,neq,nALHS) 
    
    
    
    IF (((CFLmax .gt. CFLmax_ant) .or. (ForceFact .eq. 1) .or. (N.eq.1) ) .and. (ForceFact .ne. 2)) then    !((MOD(N-1,1) .EQ. 0) .or. ((N.le.0).and.(MOD(N-1,5) .EQ. 0)))  THEN             
        OPEN (UNIT=22,FILE='./OUT1/factored.OUT', position='append', Status= 'unknown') 
        write(22,'(2e12.4E3,a,1i4)') T, CFLmax
        CLOSE(unit=22)
        Alhs(:) = 0.d0
        Clhs(:) = 0.d0
        do iel = 1, nume
            do iCol = 1, nd                                                                      !
                do iLin = 1, nd                                                                  !
                     eleffm(iLin,iCol) = stiff ((iCol-1)*nd+iLin, iel)                           ! Transforma matrizes locais para serem
                enddo                                                                            !  lidas por addlhs
            enddo                                                                                !
            call addlhs_Liad(Alhs,Clhs,eleffm,idiag,lm_Parab(1,1,iel),nd,.false.,.false.,nALHS)  !
        Enddo      
        !write (*,*) 'Factoring...'                 ! 
        call factns(alhs,clhs,idiag,neq)            !
    ENDIF                                           !
                                                    !
    do ieq = 1, neq                                 
        brhs(ieq) = u_n(ieq)        
    enddo  
   ! write (*,*) 'Backing...'                    !  Fatora a matriz de rigidez
    call backns(alhs,clhs,brhs,idiag,neq)        !    e soluciona o sistema
    do ieq = 1, neq                              !                              
        u(ieq) = brhs(ieq)                       !                              
    enddo                                        !
    
    call CroutFactOUT (idiag, neq)
    
return                                           
end
    
subroutine CroutFactOUT (idiag, neq)

    use modvar, only: ireordering

    implicit none
    
    INTEGER*8               :: idiag(neq)
    INTEGER               :: neq
    CHARACTER*70 :: FILENAME
    
    FILENAME = './OUT1/BAND.OUT'
    OPEN (UNIT=22,FILE=FILENAME) 
    write (22,2201)  idiag(neq)
    CLOSE (UNIT = 22)    
    
 2201 FORMAT ('Tamanho da meia banda:   ___ ', i10/)

    end
    
    subroutine addlhs_LIAD(alhs,clhs,eleffm,idiag,lm,nee,ldiag,lsym,nALHS)  ! Transfere os termos da matriz de rigidez do elemento para a matriz de rigidez global

             ! IN:
         !  nee:     nen*ned = nnoel*ngl = nd
         !  Lm:      entra como um pointer para onde come�a a matriz lm do elemento da vez
         !  eleffm:  Matriz de rigidez do elemento
         !  idiag:   ???
         !  (...)
         !  ldiag = .true.   -> matriz diagonal                
         !  ldiag  = .false. -> matriz n�o diagonal
             !lsym = .true.    -> matriz sim�trica
             !lsym = .false.   -> matriz n�o sim�trica
           
             ! IN / OUT:
         !  alhs:    ???
         !  clhs:    ???
             
      implicit real*8 (a-h,o-z) 
 
      logical ldiag,lsym 
      dimension idiag(*),lm(*) 
      
      integer*8 :: nALHS
      real*8 :: alhs(nALHS), clhs(nALHS), eleffm(nee,nee)
 
      if (ldiag) then ! Matriz diagonal
         do 100 j=1,nee 
            k = lm(j) 
            if (k.gt.0) then 
               l = idiag(k)
               alhs(l) = alhs(l) + eleffm(j,j) 
            endif 
 100     continue 
 
      else ! Matriz n�o diagonal
         do 400 j=1,nee 
             k = lm(j) 
             if (k.gt.0) then 
                do 200 i=1,nee ! Banda superior da matriz local
                      m = lm(i) 
                      if (m.gt.0) then 
                         if (k.ge.m) then ! Banda superior da matriz global
                            l = idiag(k) - k + m 
                            alhs(l) = alhs(l) + eleffm(i,j) 
                         else  ! Banda inferior da matriz global
                            l = idiag(m) + k - m 
                            clhs(l) = clhs(l) + eleffm(i,j) 
                         endif 
                      endif 
200             continue 
                
             endif 
 400    continue 
 
      endif 
 
      return 
end
