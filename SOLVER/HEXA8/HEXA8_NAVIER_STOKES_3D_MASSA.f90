subroutine hexa8_NavierStokes3D_Mass (stiff_Mass,x,y,z,incid,mtype,prop,numnp,nume,nummat,nd,nnoel,ncprop)

!use modtapes , only: iout,ierror
!use modmesh  , only: ID
!use modvar   , only: ngl
    
IMPLICIT NONE
     
!      Global     
INTEGER, INTENT(in)  :: nume, nd, numnp, nummat, nnoel, ncprop
INTEGER, INTENT(in)  :: mtype(nume), incid(nume,nnoel)
REAL*8 , INTENT(in)  :: x(NUMNP), y(NUMNP), z(NUMNP), prop(nummat,ncprop)
REAL*8 , INTENT(out) :: stiff_Mass(nd*nd,nume)

!      Local
INTEGER              :: i, iel, lmaux(nd), ii, jj, kk, ind, jnd, ieq, ima, inoel, ino, iNgl, col
REAL*8               :: ske_Mass (nd*nd)
REAL*8               :: rpg(8), spg(8), tpg(8), Wpg, r, s, t, xRho
REAL*8               :: x15, y15, z15,x14, y14, z14, x12, y12, z12, DetJ, DetJ_inv, vol
REAL*8               :: N1, N2, N3, N4, N5, N6, N7, N8

    
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

          Wpg = 0.12500d0

          rpg(1) =   0.2113248654
          rpg(2) =   0.7886751346
          rpg(3) =   0.7886751346
          rpg(4) =   0.2113248654
          rpg(5) =   0.2113248654
          rpg(6) =   0.7886751346
          rpg(7) =   0.7886751346
          rpg(8) =   0.2113248654

          spg(1) =   0.2113248654
          spg(2) =   0.2113248654
          spg(3) =   0.7886751346
          spg(4) =   0.7886751346
          spg(5) =   0.2113248654
          spg(6) =   0.2113248654
          spg(7) =   0.7886751346
          spg(8) =   0.7886751346

          tpg(1) =   0.2113248654
          tpg(2) =   0.2113248654
          tpg(3) =   0.2113248654
          tpg(4) =   0.2113248654
          tpg(5) =   0.7886751346
          tpg(6) =   0.7886751346
          tpg(7) =   0.7886751346
          tpg(8) =   0.7886751346
    
do iel = 1, nume ! loop nos elementos        
       
       xRho    = prop(mtype(iel),4)       ! Densidade
           
       x12 = x(incid(iel,2)) - x(incid(iel,1))  !
       x14 = x(incid(iel,4)) - x(incid(iel,1))  ! Vetores que definem o hexaedro
       x15 = x(incid(iel,5)) - x(incid(iel,1))  !
       y12 = y(incid(iel,2)) - y(incid(iel,1))  !
       y14 = y(incid(iel,4)) - y(incid(iel,1))  ! Vetores que definem o hexaedro
       y15 = y(incid(iel,5)) - y(incid(iel,1))  !
       z12 = z(incid(iel,2)) - z(incid(iel,1))  !
       z14 = z(incid(iel,4)) - z(incid(iel,1))  ! Vetores que definem o hexaedro
       z15 = z(incid(iel,5)) - z(incid(iel,1))  !
       
       DetJ = (y14*z15-z14*y15)*x12+(x15*z14-x14*z15)*y12+(-x15*y14+x14*y15)*z12     ! Determinante da transformacao em coordenadas naturais
       
       DetJ_inv = 1.d0/DetJ       ! Inverso do determinante da tranformacao
       !vol = DetJ                ! Volume do elemento
           
       ske_Mass (:) = 0.d0 
       do i =1, 8  ! loop pontos gauss - FULL INTEGRATION
             r = rpg(i)
             s = spg(i)
             t = tpg(i)

             N1 = 1 - r - s - t + r*s + s*t + t*r - r*s*t   !
             N2 = r - r*s - t*r + r*s*t                     !
             N3 = r*s - r*s*t                               !
             N4 = s - r*s - s*t + r*s*t                     !   Funcoes de forma para o octa 8
             N5 = t - s*t - t*r + r*s*t                     !
             N6 = t*r - r*s*t                               !
             N7 = r*s*t                                     !
             N8 = s*t - r*s*t                               !
             
   !   BEGIN: MATRIZ DE MASSA
       
            ske_Mass (1) =  ske_Mass (1) + N1 * N1          !
            ske_Mass (26) =  ske_Mass (26) + N1 * N1        !
            ske_Mass (51) =  ske_Mass (51) + N1 * N1        !
            ske_Mass (73) =  ske_Mass (73) + N1 * N2        !
            ske_Mass (76) =  ske_Mass (76) + N2 * N2        !  MATRIZ MASSA
            ske_Mass (98) =  ske_Mass (98) + N1 * N2        !
            ske_Mass (101) =  ske_Mass (101) + N2 * N2      !
            ske_Mass (123) =  ske_Mass (123) + N1 * N2      !
            ske_Mass (126) =  ske_Mass (126) + N2 * N2      !  MATRIZ MASSA
            ske_Mass (145) =  ske_Mass (145) + N1 * N3      !
            ske_Mass (148) =  ske_Mass (148) + N2 * N3      !
            ske_Mass (151) =  ske_Mass (151) + N3 * N3      !
            ske_Mass (170) =  ske_Mass (170) + N1 * N3      !
            ske_Mass (173) =  ske_Mass (173) + N2 * N3      !
            ske_Mass (176) =  ske_Mass (176) + N3 * N3      !
            ske_Mass (195) =  ske_Mass (195) + N1 * N3      !
            ske_Mass (198) =  ske_Mass (198) + N2 * N3      !
            ske_Mass (201) =  ske_Mass (201) + N3 * N3      !
            ske_Mass (217) =  ske_Mass (217) + N1 * N4      !
            ske_Mass (220) =  ske_Mass (220) + N2 * N4      !
            ske_Mass (223) =  ske_Mass (223) + N3 * N4      !
            ske_Mass (226) =  ske_Mass (226) + N4 * N4      !
            ske_Mass (242) =  ske_Mass (242) + N1 * N4      !
            ske_Mass (245) =  ske_Mass (245) + N2 * N4      !
            ske_Mass (248) =  ske_Mass (248) + N3 * N4      !
            ske_Mass (251) =  ske_Mass (251) + N4 * N4      !  MATRIZ MASSA
            ske_Mass (267) =  ske_Mass (267) + N1 * N4      !
            ske_Mass (270) =  ske_Mass (270) + N2 * N4      !
            ske_Mass (273) =  ske_Mass (273) + N3 * N4      !
            ske_Mass (276) =  ske_Mass (276) + N4 * N4      !  MATRIZ MASSA
            ske_Mass (289) =  ske_Mass (289) + N1 * N5      !
            ske_Mass (292) =  ske_Mass (292) + N2 * N5      !
            ske_Mass (295) =  ske_Mass (295) + N3 * N5      !
            ske_Mass (298) =  ske_Mass (298) + N4 * N5      !
            ske_Mass (301) =  ske_Mass (301) + N5 * N5      !
            ske_Mass (314) =  ske_Mass (314) + N1 * N5      !
            ske_Mass (317) =  ske_Mass (317) + N2 * N5      !
            ske_Mass (320) =  ske_Mass (320) + N3 * N5      !
            ske_Mass (323) =  ske_Mass (323) + N4 * N5      !
            ske_Mass (326) =  ske_Mass (326) + N5 * N5      !
            ske_Mass (339) =  ske_Mass (339) + N1 * N5      !
            ske_Mass (342) =  ske_Mass (342) + N2 * N5      !
            ske_Mass (345) =  ske_Mass (345) + N3 * N5      !
            ske_Mass (348) =  ske_Mass (348) + N4 * N5      !
            ske_Mass (351) =  ske_Mass (351) + N5 * N5      !
            ske_Mass (361) =  ske_Mass (361) + N1 * N6      !
            ske_Mass (364) =  ske_Mass (364) + N2 * N6      !
            ske_Mass (367) =  ske_Mass (367) + N3 * N6      !
            ske_Mass (370) =  ske_Mass (370) + N4 * N6      !
            ske_Mass (373) =  ske_Mass (373) + N5 * N6      !
            ske_Mass (376) =  ske_Mass (376) + N6 * N6      !
            ske_Mass (386) =  ske_Mass (386) + N1 * N6      !  MATRIZ MASSA
            ske_Mass (389) =  ske_Mass (389) + N2 * N6      !
            ske_Mass (392) =  ske_Mass (392) + N3 * N6      !
            ske_Mass (395) =  ske_Mass (395) + N4 * N6      !
            ske_Mass (398) =  ske_Mass (398) + N5 * N6      !  MATRIZ MASSA
            ske_Mass (401) =  ske_Mass (401) + N6 * N6      !
            ske_Mass (411) =  ske_Mass (411) + N1 * N6      !
            ske_Mass (414) =  ske_Mass (414) + N2 * N6      !
            ske_Mass (417) =  ske_Mass (417) + N3 * N6      !
            ske_Mass (420) =  ske_Mass (420) + N4 * N6      !
            ske_Mass (423) =  ske_Mass (423) + N5 * N6      !
            ske_Mass (426) =  ske_Mass (426) + N6 * N6      !
            ske_Mass (433) =  ske_Mass (433) + N1 * N7      !
            ske_Mass (436) =  ske_Mass (436) + N2 * N7      !
            ske_Mass (439) =  ske_Mass (439) + N3 * N7      !
            ske_Mass (442) =  ske_Mass (442) + N4 * N7      !
            ske_Mass (445) =  ske_Mass (445) + N5 * N7      !
            ske_Mass (448) =  ske_Mass (448) + N6 * N7      !
            ske_Mass (451) =  ske_Mass (451) + N7 * N7      !
            ske_Mass (458) =  ske_Mass (458) + N1 * N7      !
            ske_Mass (461) =  ske_Mass (461) + N2 * N7      !
            ske_Mass (464) =  ske_Mass (464) + N3 * N7      !
            ske_Mass (467) =  ske_Mass (467) + N4 * N7      !
            ske_Mass (470) =  ske_Mass (470) + N5 * N7      !
            ske_Mass (473) =  ske_Mass (473) + N6 * N7      !
            ske_Mass (476) =  ske_Mass (476) + N7 * N7      !
            ske_Mass (483) =  ske_Mass (483) + N1 * N7      !
            ske_Mass (486) =  ske_Mass (486) + N2 * N7      !
            ske_Mass (489) =  ske_Mass (489) + N3 * N7      !  MATRIZ MASSA
            ske_Mass (492) =  ske_Mass (492) + N4 * N7      !
            ske_Mass (495) =  ske_Mass (495) + N5 * N7      !
            ske_Mass (498) =  ske_Mass (498) + N6 * N7      ! MATRIZ MASSA
            ske_Mass (501) =  ske_Mass (501) + N7 * N7      !
            ske_Mass (505) =  ske_Mass (505) + N1 * N8      !
            ske_Mass (508) =  ske_Mass (508) + N2 * N8      !
            ske_Mass (511) =  ske_Mass (511) + N3 * N8      !
            ske_Mass (514) =  ske_Mass (514) + N4 * N8      !
            ske_Mass (517) =  ske_Mass (517) + N5 * N8      !
            ske_Mass (520) =  ske_Mass (520) + N6 * N8      !
            ske_Mass (523) =  ske_Mass (523) + N7 * N8      !
            ske_Mass (526) =  ske_Mass (526) + N8 * N8      !
            ske_Mass (530) =  ske_Mass (530) + N1 * N8      !
            ske_Mass (533) =  ske_Mass (533) + N2 * N8      !
            ske_Mass (536) =  ske_Mass (536) + N3 * N8      !
            ske_Mass (539) =  ske_Mass (539) + N4 * N8      !
            ske_Mass (542) =  ske_Mass (542) + N5 * N8      !
            ske_Mass (545) =  ske_Mass (545) + N6 * N8      !
            ske_Mass (548) =  ske_Mass (548) + N7 * N8      !
            ske_Mass (551) =  ske_Mass (551) + N8 * N8      !
            ske_Mass (555) =  ske_Mass (555) + N1 * N8      !
            ske_Mass (558) =  ske_Mass (558) + N2 * N8      !
            ske_Mass (561) =  ske_Mass (561) + N3 * N8      !
            ske_Mass (564) =  ske_Mass (564) + N4 * N8      !
            ske_Mass (567) =  ske_Mass (567) + N5 * N8      !
            ske_Mass (570) =  ske_Mass (570) + N6 * N8      !
            ske_Mass (573) =  ske_Mass (573) + N7 * N8      !
            ske_Mass (576) =  ske_Mass (576) + N8 * N8      !

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
