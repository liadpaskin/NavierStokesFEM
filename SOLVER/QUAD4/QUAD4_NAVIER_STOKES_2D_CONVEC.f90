subroutine quad4_NavierStokes2D_Convec (stiff,fp,up,u,lm,x,y,incid,mtype,prop,numnp,nume,&
nummat,nwk,nd,neq1,neqp,nnoel,ncprop,stiff_Convec,stiff_difus_Lambda, stiff_difus_xMu, areatot)

use modtapes , only: iout,ierror
use modmesh  , only: ID
use modvar   , only: ngl, grav
use modloads , only: fbody
    
IMPLICIT NONE
     
!      Global     
INTEGER, INTENT(in)  :: nwk, nume, nd, neq1, numnp, nummat, nnoel, ncprop, neqp
INTEGER, INTENT(in)  :: mtype(nume), lm(nume,nd), incid(nume,nnoel)
REAL*8 , INTENT(in)  :: x(NUMNP), y(NUMNP), prop(nummat,ncprop), up(0:neqp), u(0:neq1-1), &
stiff_Difus_Lambda(nd*nd,nume),stiff_Difus_xMu(nd*nd,nume)
REAL*8 , INTENT(out) :: stiff(nd*nd,nume),stiff_Convec(nd*nd,nume), fp(0:neq1-1), areatot

!      Local
INTEGER              :: i, iel, lmaux(nd), ii, jj, kk, ind, jnd, ieq, ima, inoel, ino, iNgl, col
REAL*8               :: x12, x14, y12, y14, DetJ, DetJ_inv, area
REAL*8               :: rpg(4), spg(4), Wpg, r, s, t, xRho, xMu, xLambda
REAL*8               :: dN1dx, dN2dx, dN3dx, dN4dx, dN5dx, dN6dx, dN1dy, dN2dy, dN3dy, dN4dy, dN5dy, dN6dy
REAL*8               :: N1, N2, N3, N4, N5, N6, u1, u2

real*8, allocatable  :: sk (:,:), fe(:), ue(:) , ske_Convec(:), vecforcante(:), u_ant(:,:)
allocate (sk (nd,nd))
allocate (fe(nd))
allocate (ue(nd))
allocate (ske_Convec(nd*nd))
allocate (vecforcante(nd))
allocate (u_ant(nnoel,ngl))

    

! ske(  1) : ske( 13) : ske( 25) : ske( 37) : ske( 49) : ske( 61) : ske( 73) : ske( 85) : ske( 97) : ske(109) : --- : ske(133)
! ske(  2) : ske( 14) : ske( 26) : ske( 38) : ske( 50) : ske( 62) : ske( 74) : ske( 86) : ske( 98) : ske(110) : --- : ske(134)
! ske(  3) : ske( 15) : ske( 27) : ske( 39) : ske( 51) : ske( 63) : ske( 75) : ske( 87) : ske( 99) : ske(111) : --- : ske(135)
! ske(  4) : ske( 16) : ske( 28) : ske( 40) : ske( 52) : ske( 64) : ske( 76) : ske( 88) : ske(100) : ske(112) : --- : ske(136)
! ske(  5) : ske( 17) : ske( 29) : ske( 41) : ske( 53) : ske( 65) : ske( 77) : ske( 89) : ske(101) : ske(113) : --- : ske(137)
! ske(  6) : ske( 18) : ske( 30) : ske( 42) : ske( 54) : ske( 66) : ske( 78) : ske( 90) : ske(102) : ske(114) : --- : ske(138)
! ske(  7) : ske( 19) : ske( 31) : ske( 43) : ske( 55) : ske( 67) : ske( 79) : ske( 91) : ske(103) : ske(115) : --- : ske(139)
! ske(  8) : ske( 20) : ske( 32) : ske( 44) : ske( 56) : ske( 68) : ske( 80) : ske( 92) : ske(104) : ske(116) : --- : ske(140)
! ske(  9) : ske( 21) : ske( 33) : ske( 45) : ske( 57) : ske( 69) : ske( 81) : ske( 93) : ske(105) : ske(117) : --- : ske(141)
! ske( 10) : ske( 22) : ske( 34) : ske( 46) : ske( 58) : ske( 70) : ske( 82) : ske( 94) : ske(106) : ske(118) : --- : ske(142)
! ............................................................................................................................
! ske( 12) : ske( 24) : ske( 36) : ske( 48) : ske( 60) : ske( 72) : ske( 84) : ske( 96) : ske(108) : ske(120) : --- : ske(576)


!   Pontos de Gauss:


          rpg(1) =   0.2113248654
          rpg(2) =   0.7886751346
          rpg(3) =   0.7886751346
          rpg(4) =   0.2113248654

          spg(1) =   0.2113248654
          spg(2) =   0.2113248654
          spg(3) =   0.7886751346
          spg(4) =   0.7886751346

          Wpg    = 0.2500d0 ! Fator de ponderação dos pontos de Gauss (1/4)

fbody(:)   = 0.d0
fp(:)      = 0.d0
stiff(:,:) = 0.d0
areatot    = 0.d0
    
do iel = 1, nume ! loop nos elementos
    
       ue(:)    = 0.d0      
       fe(:)    = 0.d0      
       lmaux(:) = 0      
       
       xRho    = prop(mtype(iel),4)       ! Densidade

       x12 = x(incid(iel,2)) - x(incid(iel,1))  !
       x14 = x(incid(iel,4)) - x(incid(iel,1))  ! Vetores que definem o quadrilatero
       y12 = y(incid(iel,2)) - y(incid(iel,1))  !
       y14 = y(incid(iel,4)) - y(incid(iel,1))  ! Vetores que definem o quadrilatero
       
       DetJ = x12*y14-x14*y12

       DetJ_inv = 1.d0/DetJ       ! Inverso do determinante da tranformação
       area = DetJ        ! Area do elemento
       areatot = areatot + area   ! Contador: Area total da malha
           
        if(area.le.0.d0) then         !
                write(ierror,100)iel  !  Controle, confere que a area do elemento e positiva
            stop                      !     (nos do elemento ordenados corretamente)
        endif                         !                       !
           
       ske_Convec      (:) = 0.d0 
       vecForcante     (:) = 0.d0
       do i =1, 4  ! loop pontos gauss - FULL INTEGRATION
             r = rpg(i)
             s = spg(i)
                           
             N1 = 1-s-r+r*s                         !
             N2 = r-r*s                         !
             N3 = r*s                         ! Funções de forma
             N4 = s-r*s                                !  avaliadas nos pontos de Gauss

            dN1dx = (y14  *(-1+s)-y12 *(-1+r) )*detJ_Inv
            dN1dy = (-x14 *(-1+s)+x12 *(-1+r) )*detJ_Inv
            dN2dx = (y14  *(1-s) +y12 *r      )*detJ_Inv
            dN2dy = (-x14 *(1-s) -x12 *r      )*detJ_Inv
            dN3dx = (y14  *s     -y12 *r      )*detJ_Inv
            dN3dy = (-x14 *s     +x12 *r      )*detJ_Inv
            dN4dx = (-y14 *s     -y12 *(1-r)  )*detJ_Inv
            dN4dy = (x14  *s     +x12 *(1-r)  )*detJ_Inv
                                                                                                                  !
             
             !!!!!!!!!!!!!!!!! MATRIZ CONVECTIVA
             
             do inoel = 1, nnoel                              !
                 ino  = INCID(iel,inoel)                  !
                 do ingl = 1, ngl                           !
                     ieq = ID(ingl,ino)                   !
       	             if (ieq.gt.0) then                   !  Armazena velocidades do passo anterior
       	                u_ant(inoel,ingl) = u(ieq)        !
                     else                                 !
                        u_ant(inoel,ingl) =  up(-ieq)     !
                     endif                                !
                 enddo                                    !
             enddo       
       
             u1 = (u_ant(1,1) * N1) + (u_ant(2,1) * N2) + (u_ant(3,1) * N3) + (u_ant(4,1) * N4) ! Calcula velocidades
             u2 = (u_ant(1,2) * N1) + (u_ant(2,2) * N2) + (u_ant(3,2) * N3) + (u_ant(4,2) * N4)   !  do passo anterior

            ske_Convec (1) =  ske_Convec (1) + N1 * u1 * dN1dx + N1 * u2 * dN1dy        !
            ske_Convec (3) =  ske_Convec (3) + N2 * u1 * dN1dx + N2 * u2 * dN1dy        !      MATRIZ CONVECTIVA
            ske_Convec (5) =  ske_Convec (5) + N3 * u1 * dN1dx + N3 * u2 * dN1dy        !
            ske_Convec (7) =  ske_Convec (7) + N4 * u1 * dN1dx + N4 * u2 * dN1dy        !    Non Sym
            ske_Convec (10) =  ske_Convec (10) + N1 * u1 * dN1dx + N1 * u2 * dN1dy      !
            ske_Convec (12) =  ske_Convec (12) + N2 * u1 * dN1dx + N2 * u2 * dN1dy      ! FULL INTEGRATION
            ske_Convec (14) =  ske_Convec (14) + N3 * u1 * dN1dx + N3 * u2 * dN1dy      !
            ske_Convec (16) =  ske_Convec (16) + N4 * u1 * dN1dx + N4 * u2 * dN1dy      !
            ske_Convec (17) =  ske_Convec (17) + N1 * u1 * dN2dx + N1 * u2 * dN2dy      !
            ske_Convec (19) =  ske_Convec (19) + N2 * u1 * dN2dx + N2 * u2 * dN2dy      !
            ske_Convec (21) =  ske_Convec (21) + N3 * u1 * dN2dx + N3 * u2 * dN2dy      !
            ske_Convec (23) =  ske_Convec (23) + N4 * u1 * dN2dx + N4 * u2 * dN2dy      !
            ske_Convec (26) =  ske_Convec (26) + N1 * u1 * dN2dx + N1 * u2 * dN2dy      !
            ske_Convec (28) =  ske_Convec (28) + N2 * u1 * dN2dx + N2 * u2 * dN2dy      !
            ske_Convec (30) =  ske_Convec (30) + N3 * u1 * dN2dx + N3 * u2 * dN2dy      !
            ske_Convec (32) =  ske_Convec (32) + N4 * u1 * dN2dx + N4 * u2 * dN2dy      !
            ske_Convec (33) =  ske_Convec (33) + N1 * u1 * dN3dx + N1 * u2 * dN3dy      !
            ske_Convec (35) =  ske_Convec (35) + N2 * u1 * dN3dx + N2 * u2 * dN3dy      !
            ske_Convec (37) =  ske_Convec (37) + N3 * u1 * dN3dx + N3 * u2 * dN3dy      !
            ske_Convec (39) =  ske_Convec (39) + N4 * u1 * dN3dx + N4 * u2 * dN3dy      !
            ske_Convec (42) =  ske_Convec (42) + N1 * u1 * dN3dx + N1 * u2 * dN3dy      !
            ske_Convec (44) =  ske_Convec (44) + N2 * u1 * dN3dx + N2 * u2 * dN3dy      !
            ske_Convec (46) =  ske_Convec (46) + N3 * u1 * dN3dx + N3 * u2 * dN3dy      !      MATRIZ CONVECTIVA
            ske_Convec (48) =  ske_Convec (48) + N4 * u1 * dN3dx + N4 * u2 * dN3dy      !
            ske_Convec (49) =  ske_Convec (49) + N1 * u1 * dN4dx + N1 * u2 * dN4dy      !    Non Sym
            ske_Convec (51) =  ske_Convec (51) + N2 * u1 * dN4dx + N2 * u2 * dN4dy      !
            ske_Convec (53) =  ske_Convec (53) + N3 * u1 * dN4dx + N3 * u2 * dN4dy      ! FULL INTEGRATION
            ske_Convec (55) =  ske_Convec (55) + N4 * u1 * dN4dx + N4 * u2 * dN4dy      !
            ske_Convec (58) =  ske_Convec (58) + N1 * u1 * dN4dx + N1 * u2 * dN4dy
            ske_Convec (60) =  ske_Convec (60) + N2 * u1 * dN4dx + N2 * u2 * dN4dy
            ske_Convec (62) =  ske_Convec (62) + N3 * u1 * dN4dx + N3 * u2 * dN4dy
            ske_Convec (64) =  ske_Convec (64) + N4 * u1 * dN4dx + N4 * u2 * dN4dy
       
             !!!!!!!!!!!!!!!! Vetor forcante de corpo (GRAVIDADE)
             
             !vecForcante(1)  = vecForcante(1)  + N1 * grav  (1)
             vecForcante(2)  = vecForcante(2)  + N1 * grav  !(2)
             !vecForcante(3)  = vecForcante(3)  + N2 * grav  (1)
             vecForcante(4)  = vecForcante(4)  + N2 * grav  !(2)
             !vecForcante(5)  = vecForcante(5)  + N3 * grav  (1)
             vecForcante(6)  = vecForcante(6)  + N3 * grav  !(2)
             !vecForcante(7)  = vecForcante(7)  + N4 * grav  (1)
             vecForcante(8)  = vecForcante(8)  + N4 * grav  !(2)
             
       enddo !Fim do loop nos pontos de Gauss FULL INTEGRATION
           
       !!!!! Multiplica constantes
       do ima = 1, nd * nd
           stiff_convec       (ima, iel) = ske_convec       (ima) * xRho     * area *  Wpg       !       no caso de Full Integration
       enddo
       do ind = 1, nd
           vecForcante (ind) = vecForcante(ind) * xRho * area * Wpg
       enddo
       
       !!!!! Monta Matriz difusiva para multiplicar o termo for�ante
       do ima = 1, nd*nd                                                     
           stiff       (ima, iel) =  stiff_difus_Lambda (ima, iel) + stiff_difus_xMu (ima, iel) !+ stiff_convec (ima, iel)   !!!CONVECCAO NAO ESTA SOMANDO NO TERMO FOR�ANTE
       enddo                                                                                      
       
       !!!!! Begin Avalia f = K*[up]                           !
       do ind = 1, nd                                          ! 
           ieq = lm (iel,ind)                                  !  Monta termo for�ante com  
           if (ieq.ge.0) lmaux(ind) = lm (iel,ind)             !    matriz de rigidez e "u" prescrito                                  
       enddo                                                   !                                    
       do ind =1, nd                                           !
           ieq = lm(iel,ind)                                   !
           if(ieq.lt.0) then                                   !
               ue(ind) = up(-ieq)                              !
          else                                                 !
               ue(ind) = 0.d0                                  !
           endif                                               !
       enddo                                                   !
       kk = 0                                                  !
       do jj = 1 ,nd                                           !   
            do ii = 1 ,nd                                      !   Monta termo for�ante com
               kk = kk + 1                                     !     matriz de rigidez e "u" prescrito
               sk(ii,jj) = stiff(kk,iel)                       !
            enddo                                              !
       enddo                                                   !
       do ind = 1, nd                                          !
            do jnd = 1, nd                                     !
               fe(ind) = fe(ind) + sk(ind,jnd)*ue(jnd)         !
            enddo                                              !
       enddo                                                   !
       do ind = 1, nd                                          !
            ieq = lm(iel,ind)                                  !
            if (ieq.gt.0) then                                 !
               fp(ieq) = fp(ieq) + fe(ind)                     !  Monta termo for�ante com
               fbody(ieq)    = fbody(ieq) + vecForcante(ind)   !
            endif                                              !    matriz de rigidez e "u" prescrito
       enddo                                                   !
       !!!!! end Avalia fp = K*[up]  
       
       do ima = 1, nd*nd         !  Matriz de rigidez final
           stiff (ima, iel) =  stiff (ima, iel)  + stiff_convec (ima, iel)      !   soma das parcelas convectiva e difusivas  !!! CONVECCAO NAO ESTA SOMANDO NO TERMO FORCANTE
       enddo    
       
enddo ! end loop nos elementos

!write (iout,200) areatot
!write (*,200)    areatot

deallocate (sk )
deallocate (fe)
deallocate (ue)
deallocate (ske_Convec)



deallocate (vecforcante)
deallocate (u_ant)

100 format ('***(TRIED2D.F90) area nao positivo p/ o elemento (',i8,')')
!200 format (//,' *** Total Volume of the Mesh= ', f10.5,' *** '//)

return
end subroutine
