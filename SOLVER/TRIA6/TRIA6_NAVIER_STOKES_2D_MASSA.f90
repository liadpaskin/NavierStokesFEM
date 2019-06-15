subroutine tria6_NavierStokes2D_Mass (stiff_Mass,x,y,incid,mtype,prop,numnp,nume,nummat,nd,nnoel,ncprop)

!use modtapes , only: iout,ierror
!use modmesh  , only: ID
!use modvar   , only: ngl
    
IMPLICIT NONE
     
!      Global     
INTEGER, INTENT(in)  :: nume, nd, numnp, nummat, nnoel, ncprop
INTEGER, INTENT(in)  :: mtype(nume), incid(nume,nnoel)
REAL*8 , INTENT(in)  :: x(NUMNP), y(NUMNP), prop(nummat,ncprop)
REAL*8 , INTENT(out) :: stiff_Mass(nd*nd,nume)

!      Local
INTEGER              :: i, iel, lmaux(nd), ii, jj, kk, ind, jnd, ieq, ima, inoel, ino, iNgl, col
REAL*8               :: ske_Mass (nd*nd)
REAL*8               :: x23, x13, y23, y13, x12, y12, DetJ, DetJ_inv, area
REAL*8               :: rpg(3), spg(3), Wpg, r, s, t, xRho
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


!   Pontos de Gauss:

rpg(1) = 0.666666666666667d0 ! Coordenadas ponto 1 (2/3; 1/6)
spg(1) = 0.166666666666667d0
rpg(2) = 0.166666666666667d0 ! Coordenadas ponto 2 (1/6; 2/3)
spg(2) = 0.666666666666667d0
rpg(3) = 0.166666666666667d0 ! Coordenadas ponto 3 (1/6; 1/6)
spg(3) = 0.166666666666667d0
Wpg    = 0.333333333333333d0 ! Fator de ponderação dos pontos de Gauss (1/3)
    
do iel = 1, nume ! loop nos elementos        
       
       xRho    = prop(mtype(iel),4)       ! Densidade
           
       x13 = x(incid(iel,1)) - x(incid(iel,3))  !
       x23 = x(incid(iel,2)) - x(incid(iel,3))  ! Vetores que definem o elemento triangular
       y13 = y(incid(iel,1)) - y(incid(iel,3))  !
       y23 = y(incid(iel,2)) - y(incid(iel,3))  !
       
       DetJ = + (x(incid(iel,1))*y(incid(iel,2))) + (x(incid(iel,3))*y(incid(iel,1))) + (x(incid(iel,2))*y(incid(iel,3))) &      ! Determinante da transformação em
              - (x(incid(iel,3))*y(incid(iel,2))) - (x(incid(iel,1))*y(incid(iel,3))) - (x(incid(iel,2))*y(incid(iel,1)))        !  coordenadas naturais do triangulo

       DetJ_inv = 1.d0/DetJ       ! Inverso do determinante da tranformação
       area = 0.5d0 * DetJ        ! Area do elemento

       ske_Mass (:) = 0.d0 
       do i =1, 3  ! loop pontos gauss - FULL INTEGRATION
             r = rpg(i)
             s = spg(i)
             t = 1 - r - s

             N1 = r * (2*r - 1)                         !
             N2 = s * (2*s - 1)                         !
             N3 = t * (2*t - 1)                         ! Funções de forma
             N4 = 4*s*r                                 !  avaliadas nos
             N5 = 4*s*t                                 ! pontos de Gauss
             N6 = 4*r*t
             
   !   BEGIN: MATRIZ DE MASSA
       
            ske_Mass (1) =  ske_Mass (1) + N1 * N1      !
            ske_Mass (14) =  ske_Mass (14) + N1 * N1    !  MATRIZ MASSA
            ske_Mass (25) =  ske_Mass (25) + N1 * N2    !
            ske_Mass (27) =  ske_Mass (27) + N2 * N2    !
            ske_Mass (38) =  ske_Mass (38) + N1 * N2    !
            ske_Mass (40) =  ske_Mass (40) + N2 * N2    !  MATRIZ MASSA
            ske_Mass (49) =  ske_Mass (49) + N1 * N3    !
            ske_Mass (51) =  ske_Mass (51) + N2 * N3    !
            ske_Mass (53) =  ske_Mass (53) + N3 * N3    !
            ske_Mass (62) =  ske_Mass (62) + N1 * N3    !
            ske_Mass (64) =  ske_Mass (64) + N2 * N3    !
            ske_Mass (66) =  ske_Mass (66) + N3 * N3    !
            ske_Mass (73) =  ske_Mass (73) + N1 * N4    !
            ske_Mass (75) =  ske_Mass (75) + N2 * N4    !
            ske_Mass (77) =  ske_Mass (77) + N3 * N4    !
            ske_Mass (79) =  ske_Mass (79) + N4 * N4    !
            ske_Mass (86) =  ske_Mass (86) + N1 * N4    !
            ske_Mass (88) =  ske_Mass (88) + N2 * N4    !
            ske_Mass (90) =  ske_Mass (90) + N3 * N4    !
            ske_Mass (92) =  ske_Mass (92) + N4 * N4    !
            ske_Mass (97) =  ske_Mass (97) + N1 * N5    !
            ske_Mass (99) =  ske_Mass (99) + N2 * N5    !
            ske_Mass (101) =  ske_Mass (101) + N3 * N5  !
            ske_Mass (103) =  ske_Mass (103) + N4 * N5  !
            ske_Mass (105) =  ske_Mass (105) + N5 * N5  !
            ske_Mass (110) =  ske_Mass (110) + N1 * N5  !
            ske_Mass (112) =  ske_Mass (112) + N2 * N5  !
            ske_Mass (114) =  ske_Mass (114) + N3 * N5  !  MATRIZ MASSA
            ske_Mass (116) =  ske_Mass (116) + N4 * N5  !
            ske_Mass (118) =  ske_Mass (118) + N5 * N5  !
            ske_Mass (121) =  ske_Mass (121) + N1 * N6  !
            ske_Mass (123) =  ske_Mass (123) + N2 * N6  !  MATRIZ MASSA
            ske_Mass (125) =  ske_Mass (125) + N3 * N6  !
            ske_Mass (127) =  ske_Mass (127) + N4 * N6  !
            ske_Mass (129) =  ske_Mass (129) + N5 * N6  !
            ske_Mass (131) =  ske_Mass (131) + N6 * N6  !
            ske_Mass (134) =  ske_Mass (134) + N1 * N6  !
            ske_Mass (136) =  ske_Mass (136) + N2 * N6  !
            ske_Mass (138) =  ske_Mass (138) + N3 * N6  !
            ske_Mass (140) =  ske_Mass (140) + N4 * N6  !
            ske_Mass (142) =  ske_Mass (142) + N5 * N6  !
            ske_Mass (144) =  ske_Mass (144) + N6 * N6  !

       enddo !Fim do loop nos pontos de Gauss FULL INTEGRATION
       
   !   End: Matriz de MASSA
   
       !!!!! Espelha simetria das matrizes difusivas
       do col = 1, nd                                                                            !
          do ima = col+1, nd                                                                     !
              ske_Mass (ima+((col-1)*nd)) = ske_Mass ((ima+((col-1)*nd)+((nd-1)*(ima-col))))     ! Algoritmo para espelhar as matrizes sim�tricas
          enddo                                                                                  !
       enddo                                                                                     !
   
       !!!!! Multiplica constantes
       do ima = 1, nd*nd
           stiff_Mass (ima, iel) = ske_Mass (ima) * xRho * detJ * Wpg       !  Full Integration
       enddo  
       
enddo ! end loop nos elementos

return

end subroutine
