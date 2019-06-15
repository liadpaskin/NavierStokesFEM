subroutine quad4_NavierStokes2D_Mass (stiff_Mass,x,y,incid,mtype,prop,numnp,nume,nummat,nd,nnoel,ncprop)

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
REAL*8               :: x12, x14, y12, y14, DetJ, DetJ_inv, area
REAL*8               :: rpg(4), spg(4), Wpg, r, s, t, xRho
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
       
       xRho    = prop(mtype(iel),4)       ! Densidade
           
       x12 = x(incid(iel,2)) - x(incid(iel,1))  !
       x14 = x(incid(iel,4)) - x(incid(iel,1))  ! Vetores que definem o quadrilatero
       y12 = y(incid(iel,2)) - y(incid(iel,1))  !
       y14 = y(incid(iel,4)) - y(incid(iel,1))  ! Vetores que definem o quadrilatero
       

       DetJ = x12*y14-x14*y12

       DetJ_inv = 1.d0/DetJ       ! Inverso do determinante da tranformação
       area = DetJ        ! Area do elemento

       ske_Mass (:) = 0.d0 
       do i =1, 4  ! loop pontos gauss - FULL INTEGRATION
             r = rpg(i)
             s = spg(i)

             N1 = 1-s-r+r*s                         !
             N2 = r-r*s                         !
             N3 = r*s                         ! Funções de forma
             N4 = s-r*s                                !  avaliadas nos pontos de Gauss
             
   !   BEGIN: MATRIZ DE MASSA
            ske_Mass (1) =  ske_Mass (1) + N1 * N1        !
            ske_Mass (10) =  ske_Mass (10) + N1 * N1      !
            ske_Mass (17) =  ske_Mass (17) + N1 * N2      !
            ske_Mass (19) =  ske_Mass (19) + N2 * N2      !
            ske_Mass (26) =  ske_Mass (26) + N1 * N2      !
            ske_Mass (28) =  ske_Mass (28) + N2 * N2      !
            ske_Mass (33) =  ske_Mass (33) + N1 * N3      !
            ske_Mass (35) =  ske_Mass (35) + N2 * N3      !  MATRIZ MASSA
            ske_Mass (37) =  ske_Mass (37) + N3 * N3      !
            ske_Mass (42) =  ske_Mass (42) + N1 * N3      !
            ske_Mass (44) =  ske_Mass (44) + N2 * N3      !
            ske_Mass (46) =  ske_Mass (46) + N3 * N3      !  MATRIZ MASSA
            ske_Mass (49) =  ske_Mass (49) + N1 * N4      !
            ske_Mass (51) =  ske_Mass (51) + N2 * N4      !
            ske_Mass (53) =  ske_Mass (53) + N3 * N4      !
            ske_Mass (55) =  ske_Mass (55) + N4 * N4      !
            ske_Mass (58) =  ske_Mass (58) + N1 * N4      !
            ske_Mass (60) =  ske_Mass (60) + N2 * N4      !
            ske_Mass (62) =  ske_Mass (62) + N3 * N4      !
            ske_Mass (64) =  ske_Mass (64) + N4 * N4      !

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
