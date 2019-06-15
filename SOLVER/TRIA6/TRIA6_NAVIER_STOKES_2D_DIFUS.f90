subroutine tria6_NavierStokes2D_Difus (lm,x,y,incid,mtype,prop,numnp,nume,nummat,nwk, &
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
REAL*8               :: x23, x13, y23, y13, x12, y12, DetJ, DetJ_inv, area
REAL*8               :: rpg(3), spg(3), Wpg, r, s, t, xMu, xLambda
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

rpg(1) = 0.666666666666667d0 ! Coordenadas ponto 1 (2/3; 1/6)
spg(1) = 0.166666666666667d0
rpg(2) = 0.166666666666667d0 ! Coordenadas ponto 2 (1/6; 2/3)
spg(2) = 0.666666666666667d0
rpg(3) = 0.166666666666667d0 ! Coordenadas ponto 3 (1/6; 1/6)
spg(3) = 0.166666666666667d0
Wpg    = 0.333333333333333d0 ! Fator de ponderação dos pontos de Gauss (1/3)

do iel = 1, nume ! loop nos elementos
      
       lmaux(:) = 0      
       
       xMu     = prop(mtype(iel),5)       ! Viscosidade din�mica
       xLambda = prop(mtype(iel),6)       ! Coeficiente de incompressibilidade
           
       x13 = x(incid(iel,1)) - x(incid(iel,3))  !
       x23 = x(incid(iel,2)) - x(incid(iel,3))  ! Vetores que definem o elemento triangular
       y13 = y(incid(iel,1)) - y(incid(iel,3))  !
       y23 = y(incid(iel,2)) - y(incid(iel,3))  !
       
       DetJ = + (x(incid(iel,1))*y(incid(iel,2))) + (x(incid(iel,3))*y(incid(iel,1))) + (x(incid(iel,2))*y(incid(iel,3))) &      ! Determinante da transformação em
              - (x(incid(iel,3))*y(incid(iel,2))) - (x(incid(iel,1))*y(incid(iel,3))) - (x(incid(iel,2))*y(incid(iel,1)))        !  coordenadas naturais do triangulo

       DetJ_inv = 1.d0/DetJ       ! Inverso do determinante da tranformação
       area = 0.5d0 * DetJ        ! Area do elemento
       areatot = areatot + area   ! Contador: Area total da malha

        if(area.le.0.d0) then         !
                write(ierror,100)iel  !  Controle, confere que o area do elemento e positivo
            stop                      !     (nos do elemento ordenados corretamente)
        endif                         !
           
       ske_difus_xMu   (:) = 0.d0 
       ske_difus_Lambda(:) = 0.d0 
       do i =1, 3  ! loop pontos gauss - FULL INTEGRATION
             r = rpg(i)
             s = spg(i)
             t = 1 - r - s

             dN1dx = (    (y23)*( 4*r - 1      )                               ) *  DetJ_inv    !
             dN2dx = (   -(y13)*( 4*s - 1      )                               ) *  DetJ_inv    !
             dN3dx = (    (y23)*(-3 + 4*r + 4*s) -   (y13) * (-3 + 4*r + 4*s)  ) *  DetJ_inv    !
             dN4dx = (  4*(y23)*( s            ) - 4*(y13) * ( r            )  ) *  DetJ_inv    !
             dN5dx = ( -4*(y23)*( s            ) -   (y13) * ( 4 - 4*r - 8*s)  ) *  DetJ_inv    !
             dN6dx = (    (y23)*(-8*r + 4 - 4*s) + 4*(y13) * ( r            )  ) *  DetJ_inv    !    Primeiras derivadas da função de forma
             dN1dy = (  - (x23)*( 4*r - 1      )                               ) *  DetJ_inv    !         avaliadas nos pontos de Gauss
             dN2dy = (    (x13)*( 4*s - 1      )                               ) *  DetJ_inv    !
             dN3dy = (  - (x23)*(-3 + 4*r + 4*s) +   (x13) * (-3 + 4*r + 4*s)  ) *  DetJ_inv    !
             dN4dy = ( -4*(x23)*( s            ) + 4*(x13) * ( r            )  ) *  DetJ_inv    !
             dN5dy = (  4*(x23)*( s            ) +   (x13) * ( 4 - 4*r - 8*s)  ) *  DetJ_inv    !
             dN6dy = (  - (x23)*(-8*r + 4 - 4*s) - 4*(x13) * ( r            )  ) *  DetJ_inv    !    avaliadas nos pontos de Gauss
                                                !
             !!!!!!!!!!!!!!!! Matriz DIFUSIVA (Parcela Mu)

ske_difus_xMu (1) =  ske_difus_xMu (1) + dN1dx * dN1dx + dN1dy * dN1dy
ske_difus_xMu (14) =  ske_difus_xMu (14) + dN1dx * dN1dx + dN1dy * dN1dy
ske_difus_xMu (25) =  ske_difus_xMu (25) + dN2dx * dN1dx + dN2dy * dN1dy
ske_difus_xMu (27) =  ske_difus_xMu (27) + dN2dx * dN2dx + dN2dy * dN2dy
ske_difus_xMu (38) =  ske_difus_xMu (38) + dN2dx * dN1dx + dN2dy * dN1dy
ske_difus_xMu (40) =  ske_difus_xMu (40) + dN2dx * dN2dx + dN2dy * dN2dy
ske_difus_xMu (49) =  ske_difus_xMu (49) + dN3dx * dN1dx + dN3dy * dN1dy
ske_difus_xMu (51) =  ske_difus_xMu (51) + dN3dx * dN2dx + dN3dy * dN2dy
ske_difus_xMu (53) =  ske_difus_xMu (53) + dN3dx * dN3dx + dN3dy * dN3dy
ske_difus_xMu (62) =  ske_difus_xMu (62) + dN3dx * dN1dx + dN3dy * dN1dy
ske_difus_xMu (64) =  ske_difus_xMu (64) + dN3dx * dN2dx + dN3dy * dN2dy
ske_difus_xMu (66) =  ske_difus_xMu (66) + dN3dx * dN3dx + dN3dy * dN3dy
ske_difus_xMu (73) =  ske_difus_xMu (73) + dN4dx * dN1dx + dN4dy * dN1dy
ske_difus_xMu (75) =  ske_difus_xMu (75) + dN4dx * dN2dx + dN4dy * dN2dy
ske_difus_xMu (77) =  ske_difus_xMu (77) + dN4dx * dN3dx + dN4dy * dN3dy
ske_difus_xMu (79) =  ske_difus_xMu (79) + dN4dx * dN4dx + dN4dy * dN4dy
ske_difus_xMu (86) =  ske_difus_xMu (86) + dN4dx * dN1dx + dN4dy * dN1dy
ske_difus_xMu (88) =  ske_difus_xMu (88) + dN4dx * dN2dx + dN4dy * dN2dy
ske_difus_xMu (90) =  ske_difus_xMu (90) + dN4dx * dN3dx + dN4dy * dN3dy
ske_difus_xMu (92) =  ske_difus_xMu (92) + dN4dx * dN4dx + dN4dy * dN4dy
ske_difus_xMu (97) =  ske_difus_xMu (97) + dN5dx * dN1dx + dN5dy * dN1dy
ske_difus_xMu (99) =  ske_difus_xMu (99) + dN5dx * dN2dx + dN5dy * dN2dy
ske_difus_xMu (101) =  ske_difus_xMu (101) + dN5dx * dN3dx + dN5dy * dN3dy
ske_difus_xMu (103) =  ske_difus_xMu (103) + dN5dx * dN4dx + dN5dy * dN4dy
ske_difus_xMu (105) =  ske_difus_xMu (105) + dN5dx * dN5dx + dN5dy * dN5dy
ske_difus_xMu (110) =  ske_difus_xMu (110) + dN5dx * dN1dx + dN5dy * dN1dy
ske_difus_xMu (112) =  ske_difus_xMu (112) + dN5dx * dN2dx + dN5dy * dN2dy
ske_difus_xMu (114) =  ske_difus_xMu (114) + dN5dx * dN3dx + dN5dy * dN3dy
ske_difus_xMu (116) =  ske_difus_xMu (116) + dN5dx * dN4dx + dN5dy * dN4dy
ske_difus_xMu (118) =  ske_difus_xMu (118) + dN5dx * dN5dx + dN5dy * dN5dy
ske_difus_xMu (121) =  ske_difus_xMu (121) + dN6dx * dN1dx + dN6dy * dN1dy
ske_difus_xMu (123) =  ske_difus_xMu (123) + dN6dx * dN2dx + dN6dy * dN2dy
ske_difus_xMu (125) =  ske_difus_xMu (125) + dN6dx * dN3dx + dN6dy * dN3dy
ske_difus_xMu (127) =  ske_difus_xMu (127) + dN6dx * dN4dx + dN6dy * dN4dy
ske_difus_xMu (129) =  ske_difus_xMu (129) + dN6dx * dN5dx + dN6dy * dN5dy
ske_difus_xMu (131) =  ske_difus_xMu (131) + dN6dx * dN6dx + dN6dy * dN6dy
ske_difus_xMu (134) =  ske_difus_xMu (134) + dN6dx * dN1dx + dN6dy * dN1dy
ske_difus_xMu (136) =  ske_difus_xMu (136) + dN6dx * dN2dx + dN6dy * dN2dy
ske_difus_xMu (138) =  ske_difus_xMu (138) + dN6dx * dN3dx + dN6dy * dN3dy
ske_difus_xMu (140) =  ske_difus_xMu (140) + dN6dx * dN4dx + dN6dy * dN4dy
ske_difus_xMu (142) =  ske_difus_xMu (142) + dN6dx * dN5dx + dN6dy * dN5dy
ske_difus_xMu (144) =  ske_difus_xMu (144) + dN6dx * dN6dx + dN6dy * dN6dy
             
      enddo !Fim do loop nos pontos de Gauss FULL INTEGRATION

           ! Begin: Matriz DIFUSIVA (Parcela LAMBDA) integra��o reduzida

             r = 0.333333333333333d0      ! PONTO DE GAUSS PARA INTEGRAÇÂO REDUZIDA
             s = 0.333333333333333d0      !

             dN1dx = (    (y23)*( 4*r - 1      )                               ) *  DetJ_inv    !
             dN2dx = (   -(y13)*( 4*s - 1      )                               ) *  DetJ_inv    !
             dN3dx = (    (y23)*(-3 + 4*r + 4*s) -   (y13) * (-3 + 4*r + 4*s)  ) *  DetJ_inv    !
             dN4dx = (  4*(y23)*( s            ) - 4*(y13) * ( r            )  ) *  DetJ_inv    !
             dN5dx = ( -4*(y23)*( s            ) -   (y13) * ( 4 - 4*r - 8*s)  ) *  DetJ_inv    !
             dN6dx = (    (y23)*(-8*r + 4 - 4*s) + 4*(y13) * ( r            )  ) *  DetJ_inv    !    Primeiras derivadas da função de forma
             dN1dy = (  - (x23)*( 4*r - 1      )                               ) *  DetJ_inv    !         avaliadas nos pontos de Gauss
             dN2dy = (    (x13)*( 4*s - 1      )                               ) *  DetJ_inv    !
             dN3dy = (  - (x23)*(-3 + 4*r + 4*s) +   (x13) * (-3 + 4*r + 4*s)  ) *  DetJ_inv    !
             dN4dy = ( -4*(x23)*( s            ) + 4*(x13) * ( r            )  ) *  DetJ_inv    !
             dN5dy = (  4*(x23)*( s            ) +   (x13) * ( 4 - 4*r - 8*s)  ) *  DetJ_inv    !
             dN6dy = (  - (x23)*(-8*r + 4 - 4*s) - 4*(x13) * ( r            )  ) *  DetJ_inv    !    avaliadas nos pontos de Gauss

       
            ske_difus_Lambda (1) =  ske_difus_Lambda (1) + dN1dx * dN1dx             !
            ske_difus_Lambda (13) =  ske_difus_Lambda (13) + dN1dy * dN1dx           !
            ske_difus_Lambda (14) =  ske_difus_Lambda (14) + dN1dy * dN1dy           ! MATRIZ DIFUSIVA
            ske_difus_Lambda (25) =  ske_difus_Lambda (25) + dN2dx * dN1dx           !
            ske_difus_Lambda (26) =  ske_difus_Lambda (26) + dN2dx * dN1dy           !    PARCELA LAMBDA
            ske_difus_Lambda (27) =  ske_difus_Lambda (27) + dN2dx * dN2dx           !
            ske_difus_Lambda (37) =  ske_difus_Lambda (37) + dN2dy * dN1dx           ! REDUCED INTEGRATION
            ske_difus_Lambda (38) =  ske_difus_Lambda (38) + dN2dy * dN1dy           !
            ske_difus_Lambda (39) =  ske_difus_Lambda (39) + dN2dy * dN2dx           !
            ske_difus_Lambda (40) =  ske_difus_Lambda (40) + dN2dy * dN2dy           !
            ske_difus_Lambda (49) =  ske_difus_Lambda (49) + dN3dx * dN1dx           !
            ske_difus_Lambda (50) =  ske_difus_Lambda (50) + dN3dx * dN1dy           !
            ske_difus_Lambda (51) =  ske_difus_Lambda (51) + dN3dx * dN2dx           !
            ske_difus_Lambda (52) =  ske_difus_Lambda (52) + dN3dx * dN2dy           !
            ske_difus_Lambda (53) =  ske_difus_Lambda (53) + dN3dx * dN3dx           !
            ske_difus_Lambda (61) =  ske_difus_Lambda (61) + dN3dy * dN1dx           !
            ske_difus_Lambda (62) =  ske_difus_Lambda (62) + dN3dy * dN1dy           !
            ske_difus_Lambda (63) =  ske_difus_Lambda (63) + dN3dy * dN2dx           !
            ske_difus_Lambda (64) =  ske_difus_Lambda (64) + dN3dy * dN2dy           !
            ske_difus_Lambda (65) =  ske_difus_Lambda (65) + dN3dy * dN3dx           !
            ske_difus_Lambda (66) =  ske_difus_Lambda (66) + dN3dy * dN3dy           !
            ske_difus_Lambda (73) =  ske_difus_Lambda (73) + dN4dx * dN1dx           !
            ske_difus_Lambda (74) =  ske_difus_Lambda (74) + dN4dx * dN1dy           !
            ske_difus_Lambda (75) =  ske_difus_Lambda (75) + dN4dx * dN2dx           ! MATRIZ DIFUSIVA
            ske_difus_Lambda (76) =  ske_difus_Lambda (76) + dN4dx * dN2dy           !
            ske_difus_Lambda (77) =  ske_difus_Lambda (77) + dN4dx * dN3dx           !    PARCELA LAMBDA
            ske_difus_Lambda (78) =  ske_difus_Lambda (78) + dN4dx * dN3dy           !
            ske_difus_Lambda (79) =  ske_difus_Lambda (79) + dN4dx * dN4dx           ! REDUCED INTEGRATION
            ske_difus_Lambda (85) =  ske_difus_Lambda (85) + dN4dy * dN1dx           !
            ske_difus_Lambda (86) =  ske_difus_Lambda (86) + dN4dy * dN1dy           !
            ske_difus_Lambda (87) =  ske_difus_Lambda (87) + dN4dy * dN2dx           !
            ske_difus_Lambda (88) =  ske_difus_Lambda (88) + dN4dy * dN2dy           !
            ske_difus_Lambda (89) =  ske_difus_Lambda (89) + dN4dy * dN3dx           !
            ske_difus_Lambda (90) =  ske_difus_Lambda (90) + dN4dy * dN3dy
            ske_difus_Lambda (91) =  ske_difus_Lambda (91) + dN4dy * dN4dx           !
            ske_difus_Lambda (92) =  ske_difus_Lambda (92) + dN4dy * dN4dy           !
            ske_difus_Lambda (97) =  ske_difus_Lambda (97) + dN5dx * dN1dx           !
            ske_difus_Lambda (98) =  ske_difus_Lambda (98) + dN5dx * dN1dy           !
            ske_difus_Lambda (99) =  ske_difus_Lambda (99) + dN5dx * dN2dx           ! MATRIZ DIFUSIVA
            ske_difus_Lambda (100) =  ske_difus_Lambda (100) + dN5dx * dN2dy         !
            ske_difus_Lambda (101) =  ske_difus_Lambda (101) + dN5dx * dN3dx         !    PARCELA LAMBDA
            ske_difus_Lambda (102) =  ske_difus_Lambda (102) + dN5dx * dN3dy         !
            ske_difus_Lambda (103) =  ske_difus_Lambda (103) + dN5dx * dN4dx         ! REDUCED INTEGRATION
            ske_difus_Lambda (104) =  ske_difus_Lambda (104) + dN5dx * dN4dy         !
            ske_difus_Lambda (105) =  ske_difus_Lambda (105) + dN5dx * dN5dx         !
            ske_difus_Lambda (109) =  ske_difus_Lambda (109) + dN5dy * dN1dx         !
            ske_difus_Lambda (110) =  ske_difus_Lambda (110) + dN5dy * dN1dy         !
            ske_difus_Lambda (111) =  ske_difus_Lambda (111) + dN5dy * dN2dx         !
            ske_difus_Lambda (112) =  ske_difus_Lambda (112) + dN5dy * dN2dy         !
            ske_difus_Lambda (113) =  ske_difus_Lambda (113) + dN5dy * dN3dx         !
            ske_difus_Lambda (114) =  ske_difus_Lambda (114) + dN5dy * dN3dy         !
            ske_difus_Lambda (115) =  ske_difus_Lambda (115) + dN5dy * dN4dx         !
            ske_difus_Lambda (116) =  ske_difus_Lambda (116) + dN5dy * dN4dy         !
            ske_difus_Lambda (117) =  ske_difus_Lambda (117) + dN5dy * dN5dx         !
            ske_difus_Lambda (118) =  ske_difus_Lambda (118) + dN5dy * dN5dy         !
            ske_difus_Lambda (121) =  ske_difus_Lambda (121) + dN6dx * dN1dx         !
            ske_difus_Lambda (122) =  ske_difus_Lambda (122) + dN6dx * dN1dy         !
            ske_difus_Lambda (123) =  ske_difus_Lambda (123) + dN6dx * dN2dx         !
            ske_difus_Lambda (124) =  ske_difus_Lambda (124) + dN6dx * dN2dy         !
            ske_difus_Lambda (125) =  ske_difus_Lambda (125) + dN6dx * dN3dx         ! MATRIZ DIFUSIVA
            ske_difus_Lambda (126) =  ske_difus_Lambda (126) + dN6dx * dN3dy         !
            ske_difus_Lambda (127) =  ske_difus_Lambda (127) + dN6dx * dN4dx         !    PARCELA LAMBDA
            ske_difus_Lambda (128) =  ske_difus_Lambda (128) + dN6dx * dN4dy         !
            ske_difus_Lambda (129) =  ske_difus_Lambda (129) + dN6dx * dN5dx         ! REDUCED INTEGRATION
            ske_difus_Lambda (130) =  ske_difus_Lambda (130) + dN6dx * dN5dy         !
            ske_difus_Lambda (131) =  ske_difus_Lambda (131) + dN6dx * dN6dx         !
            ske_difus_Lambda (133) =  ske_difus_Lambda (133) + dN6dy * dN1dx         !
            ske_difus_Lambda (134) =  ske_difus_Lambda (134) + dN6dy * dN1dy         !
            ske_difus_Lambda (135) =  ske_difus_Lambda (135) + dN6dy * dN2dx         !
            ske_difus_Lambda (136) =  ske_difus_Lambda (136) + dN6dy * dN2dy         !
            ske_difus_Lambda (137) =  ske_difus_Lambda (137) + dN6dy * dN3dx         !
            ske_difus_Lambda (138) =  ske_difus_Lambda (138) + dN6dy * dN3dy         !
            ske_difus_Lambda (139) =  ske_difus_Lambda (139) + dN6dy * dN4dx         !
            ske_difus_Lambda (140) =  ske_difus_Lambda (140) + dN6dy * dN4dy         !
            ske_difus_Lambda (141) =  ske_difus_Lambda (141) + dN6dy * dN5dx         !
            ske_difus_Lambda (142) =  ske_difus_Lambda (142) + dN6dy * dN5dy         !
            ske_difus_Lambda (143) =  ske_difus_Lambda (143) + dN6dy * dN6dx         !
            ske_difus_Lambda (144) =  ske_difus_Lambda (144) + dN6dy * dN6dy         !
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
