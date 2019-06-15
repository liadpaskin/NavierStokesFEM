subroutine quad4_NavierStokes2D_Difus (lm,x,y,incid,mtype,prop,numnp,nume,nummat,nwk, &
nd,neq1,neqp,nnoel,ncprop,stiff_difus_Lambda, stiff_difus_xMu, areatot)

use modtapes , only: iout,ierror
use modmesh  , only: ID
use modvar   , only: ngl
    
IMPLICIT NONE
     
!      Global     
INTEGER, INTENT(in)  :: nwk, nume, nd, neq1, numnp, nummat, nnoel, ncprop, neqp
INTEGER, INTENT(in)  :: mtype(nume), lm(nume,nd), incid(nume,nnoel)
REAL*8 , INTENT(in)  :: x(NUMNP), y(NUMNP), prop(nummat,ncprop)
REAL*8 , INTENT(out) :: stiff_Difus_Lambda(nd*nd,nume),stiff_Difus_xMu(nd*nd,nume), areatot

!      Local
INTEGER              :: i, iel, lmaux(nd), ii, jj, kk, ind, jnd, ieq, ima, inoel, ino, iNgl, col
REAL*8               :: ske_difus_xMu(nd*nd),ske_difus_Lambda(nd*nd)
REAL*8               :: x12, x14, y12, y14, DetJ, DetJ_inv, area
REAL*8               :: rpg(4), spg(4), Wpg, r, s, t, xMu, xLambda
REAL*8               :: dN1dx, dN2dx, dN3dx, dN4dx, dN5dx, dN6dx, dN1dy, dN2dy, dN3dy, dN4dy, dN5dy, dN6dy
REAL*8               :: N1, N2, N3, N4, N5, N6

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

          rpg(1) =   0.2113248654
          rpg(2) =   0.7886751346
          rpg(3) =   0.7886751346
          rpg(4) =   0.2113248654

          spg(1) =   0.2113248654
          spg(2) =   0.2113248654
          spg(3) =   0.7886751346
          spg(4) =   0.7886751346

          Wpg    = 0.2500d0 ! Fator de ponderação dos pontos de Gauss (1/4)

do iel = 1, nume ! loop nos elementos
      
       lmaux(:) = 0      
       
       xMu     = prop(mtype(iel),5)       ! Viscosidade din�mica
       xLambda = prop(mtype(iel),6)       ! Coeficiente de incompressibilidade
           
       x12 = x(incid(iel,2)) - x(incid(iel,1))  !
       x14 = x(incid(iel,4)) - x(incid(iel,1))  ! Vetores que definem o quadrilatero
       y12 = y(incid(iel,2)) - y(incid(iel,1))  !
       y14 = y(incid(iel,4)) - y(incid(iel,1))  ! Vetores que definem o quadrilatero

       DetJ = x12*y14-x14*y12
       
       DetJ_inv = 1.d0/DetJ       ! Inverso do determinante da tranformação
       area = DetJ        ! Area do elemento
       areatot = areatot + area   ! Contador: Area total da malha

        if(area.le.0.d0) then         !
                write(ierror,100)iel  !  Controle, confere que o area do elemento e positivo
            stop                      !     (nos do elemento ordenados corretamente)
        endif                         !
           
       ske_difus_xMu   (:) = 0.d0 
       ske_difus_Lambda(:) = 0.d0 
       do i =1, 4  ! loop pontos gauss - FULL INTEGRATION
             r = rpg(i)
             s = spg(i)

            dN1dx = (y14  *(-1+s)-y12 *(-1+r) )*detJ_Inv
            dN1dy = (-x14 *(-1+s)+x12 *(-1+r) )*detJ_Inv
            dN2dx = (y14  *(1-s) +y12 *r      )*detJ_Inv
            dN2dy = (-x14 *(1-s) -x12 *r      )*detJ_Inv
            dN3dx = (y14  *s     -y12 *r      )*detJ_Inv
            dN3dy = (-x14 *s     +x12 *r      )*detJ_Inv
            dN4dx = (-y14 *s     -y12 *(1-r)  )*detJ_Inv
            dN4dy = (x14  *s     +x12 *(1-r)  )*detJ_Inv
                                                !
             !!!!!!!!!!!!!!!! Matriz DIFUSIVA (Parcela Mu)


            ske_difus_xMu (1) =  ske_difus_xMu (1) + dN1dx * dN1dx + dN1dy * dN1dy
            ske_difus_xMu (10) =  ske_difus_xMu (10) + dN1dx * dN1dx + dN1dy * dN1dy
            ske_difus_xMu (17) =  ske_difus_xMu (17) + dN2dx * dN1dx + dN2dy * dN1dy
            ske_difus_xMu (19) =  ske_difus_xMu (19) + dN2dx * dN2dx + dN2dy * dN2dy
            ske_difus_xMu (26) =  ske_difus_xMu (26) + dN2dx * dN1dx + dN2dy * dN1dy
            ske_difus_xMu (28) =  ske_difus_xMu (28) + dN2dx * dN2dx + dN2dy * dN2dy
            ske_difus_xMu (33) =  ske_difus_xMu (33) + dN3dx * dN1dx + dN3dy * dN1dy
            ske_difus_xMu (35) =  ske_difus_xMu (35) + dN3dx * dN2dx + dN3dy * dN2dy
            ske_difus_xMu (37) =  ske_difus_xMu (37) + dN3dx * dN3dx + dN3dy * dN3dy
            ske_difus_xMu (42) =  ske_difus_xMu (42) + dN3dx * dN1dx + dN3dy * dN1dy
            ske_difus_xMu (44) =  ske_difus_xMu (44) + dN3dx * dN2dx + dN3dy * dN2dy
            ske_difus_xMu (46) =  ske_difus_xMu (46) + dN3dx * dN3dx + dN3dy * dN3dy
            ske_difus_xMu (49) =  ske_difus_xMu (49) + dN4dx * dN1dx + dN4dy * dN1dy
            ske_difus_xMu (51) =  ske_difus_xMu (51) + dN4dx * dN2dx + dN4dy * dN2dy
            ske_difus_xMu (53) =  ske_difus_xMu (53) + dN4dx * dN3dx + dN4dy * dN3dy
            ske_difus_xMu (55) =  ske_difus_xMu (55) + dN4dx * dN4dx + dN4dy * dN4dy
            ske_difus_xMu (58) =  ske_difus_xMu (58) + dN4dx * dN1dx + dN4dy * dN1dy
            ske_difus_xMu (60) =  ske_difus_xMu (60) + dN4dx * dN2dx + dN4dy * dN2dy
            ske_difus_xMu (62) =  ske_difus_xMu (62) + dN4dx * dN3dx + dN4dy * dN3dy
            ske_difus_xMu (64) =  ske_difus_xMu (64) + dN4dx * dN4dx + dN4dy * dN4dy
             
      enddo !Fim do loop nos pontos de Gauss FULL INTEGRATION

           ! Begin: Matriz DIFUSIVA (Parcela LAMBDA) integra��o reduzida

             r = 0.5d0      ! PONTO DE GAUSS PARA INTEGRAÇÂO REDUZIDA
             s = 0.5d0      !

            dN1dx = (y14  *(-1+s)-y12 *(-1+r) )*detJ_Inv
            dN1dy = (-x14 *(-1+s)+x12 *(-1+r) )*detJ_Inv
            dN2dx = (y14  *(1-s) +y12 *r      )*detJ_Inv
            dN2dy = (-x14 *(1-s) -x12 *r      )*detJ_Inv
            dN3dx = (y14  *s     -y12 *r      )*detJ_Inv
            dN3dy = (-x14 *s     +x12 *r      )*detJ_Inv
            dN4dx = (-y14 *s     -y12 *(1-r)  )*detJ_Inv
            dN4dy = (x14  *s     +x12 *(1-r)  )*detJ_Inv

       
            ske_difus_Lambda (1) =  ske_difus_Lambda (1) + dN1dx * dN1dx       !
            ske_difus_Lambda (9) =  ske_difus_Lambda (9) + dN1dy * dN1dx       !
            ske_difus_Lambda (10) =  ske_difus_Lambda (10) + dN1dy * dN1dy     ! MATRIZ DIFUSIVA
            ske_difus_Lambda (17) =  ske_difus_Lambda (17) + dN2dx * dN1dx     !
            ske_difus_Lambda (18) =  ske_difus_Lambda (18) + dN2dx * dN1dy     !    PARCELA LAMBDA
            ske_difus_Lambda (19) =  ske_difus_Lambda (19) + dN2dx * dN2dx     !
            ske_difus_Lambda (25) =  ske_difus_Lambda (25) + dN2dy * dN1dx     ! REDUCED INTEGRATION
            ske_difus_Lambda (26) =  ske_difus_Lambda (26) + dN2dy * dN1dy     !
            ske_difus_Lambda (27) =  ske_difus_Lambda (27) + dN2dy * dN2dx     !
            ske_difus_Lambda (28) =  ske_difus_Lambda (28) + dN2dy * dN2dy     !
            ske_difus_Lambda (33) =  ske_difus_Lambda (33) + dN3dx * dN1dx     !
            ske_difus_Lambda (34) =  ske_difus_Lambda (34) + dN3dx * dN1dy     !
            ske_difus_Lambda (35) =  ske_difus_Lambda (35) + dN3dx * dN2dx     !
            ske_difus_Lambda (36) =  ske_difus_Lambda (36) + dN3dx * dN2dy     !
            ske_difus_Lambda (37) =  ske_difus_Lambda (37) + dN3dx * dN3dx     !
            ske_difus_Lambda (41) =  ske_difus_Lambda (41) + dN3dy * dN1dx     !
            ske_difus_Lambda (42) =  ske_difus_Lambda (42) + dN3dy * dN1dy     !
            ske_difus_Lambda (43) =  ske_difus_Lambda (43) + dN3dy * dN2dx     !
            ske_difus_Lambda (44) =  ske_difus_Lambda (44) + dN3dy * dN2dy     !
            ske_difus_Lambda (45) =  ske_difus_Lambda (45) + dN3dy * dN3dx     !
            ske_difus_Lambda (46) =  ske_difus_Lambda (46) + dN3dy * dN3dy     !
            ske_difus_Lambda (49) =  ske_difus_Lambda (49) + dN4dx * dN1dx     !
            ske_difus_Lambda (50) =  ske_difus_Lambda (50) + dN4dx * dN1dy     !
            ske_difus_Lambda (51) =  ske_difus_Lambda (51) + dN4dx * dN2dx     ! MATRIZ DIFUSIVA
            ske_difus_Lambda (52) =  ske_difus_Lambda (52) + dN4dx * dN2dy     !
            ske_difus_Lambda (53) =  ske_difus_Lambda (53) + dN4dx * dN3dx     !    PARCELA LAMBDA
            ske_difus_Lambda (54) =  ske_difus_Lambda (54) + dN4dx * dN3dy     !
            ske_difus_Lambda (55) =  ske_difus_Lambda (55) + dN4dx * dN4dx     ! REDUCED INTEGRATION
            ske_difus_Lambda (57) =  ske_difus_Lambda (57) + dN4dy * dN1dx     !
            ske_difus_Lambda (58) =  ske_difus_Lambda (58) + dN4dy * dN1dy     !
            ske_difus_Lambda (59) =  ske_difus_Lambda (59) + dN4dy * dN2dx     !
            ske_difus_Lambda (60) =  ske_difus_Lambda (60) + dN4dy * dN2dy     !
            ske_difus_Lambda (61) =  ske_difus_Lambda (61) + dN4dy * dN3dx     !
            ske_difus_Lambda (62) =  ske_difus_Lambda (62) + dN4dy * dN3dy
            ske_difus_Lambda (63) =  ske_difus_Lambda (63) + dN4dy * dN4dx
            ske_difus_Lambda (64) =  ske_difus_Lambda (64) + dN4dy * dN4dy
                                                                                     !
       
       !   End: Matriz DIFUSIVA (Parcela LAMBDA) integra��o reduzida

           !!!!! Espelha simetria das matrizes difusivas
           do col = 1, nd                                                                                             !
              do ima = col+1, nd                                                                                      !
                  ske_difus_Lambda (ima+((col-1)*nd)) = ske_difus_Lambda ((ima+((col-1)*nd)+((nd-1)*(ima-col))))      ! Algoritmo para espelhar as matrizes
                  ske_difus_xMu    (ima+((col-1)*nd)) = ske_difus_xMu    ((ima+((col-1)*nd)+((nd-1)*(ima-col))))      !    sim�tricas
              enddo                                                                                                   !
           enddo                                                                                                      !
       
           !!!!! Multiplica constantes
           do ima = 1, nd * nd
               stiff_difus_Lambda (ima, iel) = ske_difus_Lambda (ima) * xLambda  * area              !  Multiplica constantes de cada matriz
               stiff_difus_xMu    (ima, iel) = ske_difus_xMu    (ima) * xMu      * area *  Wpg       !     inclusive a pondera��o de Gauss
           enddo
       
enddo ! end loop nos elementos

write (iout,200) areatot
write (*,200)    areatot

return

100 format ('***(TRIED2D.F90) area nao positivo p/ o elemento (',i8,')')
200 format (//,' *** Total area of the Mesh= ', f10.5,' *** '//)

end subroutine
