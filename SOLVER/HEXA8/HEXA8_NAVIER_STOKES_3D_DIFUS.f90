subroutine hexa8_NavierStokes3D_Difus (lm,x,y,z,incid,mtype,prop,numnp,nume,nummat,nwk, &
nd,neq1,neqp,nnoel,ncprop,stiff_difus_Lambda, stiff_difus_xMu, voltot)

use modtapes , only: iout,ierror
use modmesh  , only: ID
use modvar   , only: ngl
    
IMPLICIT NONE
     
!      Global     
INTEGER, INTENT(in)  :: nwk, nume, nd, neq1, numnp, nummat, nnoel, ncprop, neqp
INTEGER, INTENT(in)  :: mtype(nume), lm(nume,nd), incid(nume,nnoel)
REAL*8 , INTENT(in)  :: x(NUMNP), y(NUMNP), z(NUMNP), prop(nummat,ncprop)
REAL*8 , INTENT(out) :: stiff_Difus_Lambda(nd*nd,nume),stiff_Difus_xMu(nd*nd,nume), voltot

!      Local
INTEGER              :: i, iel, lmaux(nd), ii, jj, kk, ind, jnd, ieq, ima, inoel, ino, iNgl, col
REAL*8               :: ske_difus_xMu(nd*nd),ske_difus_Lambda(nd*nd)
REAL*8               :: rpg(8), spg(8), tpg(8), Wpg, r, s, t, xMu, xLambda
REAL*8               :: x15, y15, z15,x14, y14, z14, x12, y12, z12, DetJ, DetJ_inv, vol
REAL*8               :: dN1dx, dN2dx, dN3dx, dN4dx, dN5dx, dN6dx, dN7dx, dN8dx, &
                        dN1dy, dN2dy, dN3dy, dN4dy, dN5dy, dN6dy, dN7dy, dN8dy, &
                        dN1dz, dN2dz, dN3dz, dN4dz, dN5dz, dN6dz, dN7dz, dN8dz
REAL*8               :: N1, N2, N3, N4, N5, N6, N7, N8
REAL*8               :: yz1512, xz1514, xy1514, yz1415, xz1215, xy1215, yz1214, xz1214, xy1214

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

          Wpg = 0.12500d0

!   Pontos de Gauss:
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
      
       lmaux(:) = 0      
       
       xMu     = prop(mtype(iel),5)       ! Viscosidade din�mica
       xLambda = prop(mtype(iel),6)       ! Coeficiente de incompressibilidade
           
       x12 = x(incid(iel,2)) - x(incid(iel,1))  !
       x14 = x(incid(iel,4)) - x(incid(iel,1))  ! Vetores que definem o hexaedro
       x15 = x(incid(iel,5)) - x(incid(iel,1))  !
       y12 = y(incid(iel,2)) - y(incid(iel,1))  !
       y14 = y(incid(iel,4)) - y(incid(iel,1))  ! Vetores que definem o hexaedro
       y15 = y(incid(iel,5)) - y(incid(iel,1))  !
       z12 = z(incid(iel,2)) - z(incid(iel,1))  !
       z14 = z(incid(iel,4)) - z(incid(iel,1))  ! Vetores que definem o hexaedro
       z15 = z(incid(iel,5)) - z(incid(iel,1))  !
       
       yz1512 = (y15*z12-y12*z15)      ! Caracteristicas
       xz1514 = (x15*z14-x14*z15)      !    Geometricas
       xy1514 = (-x15*y14+x14*y15)     !   Para Jacobiano
       yz1415 = (y14*z15-z14*y15)      !
       xz1215 = (x12*z15-z12*x15)      !
       xy1215 = (-x12*y15+y12*x15)     !
       yz1214 = (y12*z14-y14*z12)      !
       xz1214 = (-x12*z14+z12*x14)     !
       xy1214 = (x12*y14-y12*x14)      !

       DetJ = (y14*z15-z14*y15)*x12+(x15*z14-x14*z15)*y12+(-x15*y14+x14*y15)*z12     ! Determinante da transformacao em coordenadas naturais
       
       DetJ_inv = 1.d0/DetJ       ! Inverso do determinante da tranformacao
       vol = DetJ                ! Volume do elemento
       voltot = voltot + vol   ! Contador: volume total da malha
           
        if(vol.le.0.d0) then         !
                write(ierror,100)iel  !  Controle, confere que o volume do elemento e positivo
            stop                      !     (nos do elemento ordenados corretamente)
        endif                         !
           
       ske_difus_xMu   (:) = 0.d0 
       ske_difus_Lambda(:) = 0.d0 
       do i =1, 8  ! loop pontos gauss - FULL INTEGRATION
             r = rpg(i)
             s = spg(i)
             t = tpg(i)

            dN1dx = (yz1415*(-1+s+t-s*t)+yz1512*(-1+r+t-t*r)+yz1214*(-1+s+r-r*s))*DetJ_inv   !
            dN1dy = (xz1514*(-1+s+t-s*t)+xz1215*(-1+r+t-t*r)+xz1214*(-1+s+r-r*s))*DetJ_inv  !
            dN1dz = (xy1514*(-1+s+t-s*t)+xy1215*(-1+r+t-t*r)+xy1214*(-1+s+r-r*s))*DetJ_inv  !
                                                                                                                     !
            dN2dx = (yz1415*(1-s-t+s*t)+yz1512*(-r+t*r)+yz1214*(-r+r*s)         )*DetJ_inv  !
            dN2dy = (xz1514*(1-s-t+s*t)+xz1215*(-r+t*r)+xz1214*(-r+r*s)         )*DetJ_inv   !    Primeiras derivadas da funcao de forma
            dN2dz = (xy1514*(1-s-t+s*t)+xy1215*(-r+t*r)+xy1214*(-r+r*s)         )*DetJ_inv  !         avaliadas nos pontos de Gauss
                                                                                                                     !
            dN3dx = (yz1415*(s-s*t)+yz1512*(r-t*r)-yz1214*r*s                   )*DetJ_inv   !
            dN3dy = (xz1514*(s-s*t)+xz1215*(r-t*r)-xz1214*r*s                   )*DetJ_inv   !
            dN3dz = (xy1514*(s-s*t)+xy1215*(r-t*r)-xy1214*r*s                   )*DetJ_inv  !
                                                                                                                    !
            dN4dx = (yz1415*(-s+s*t)+yz1512*(1-r-t+t*r)+yz1214*(-s+r*s)         )*DetJ_inv   !
            dN4dy = (xz1514*(-s+s*t)+xz1215*(1-r-t+t*r)+xz1214*(-s+r*s)         )*DetJ_inv   !
            dN4dz = (xy1514*(-s+s*t)+xy1215*(1-r-t+t*r)+xy1214*(-s+r*s)         )*DetJ_inv  !
                                                                                                                     !
            dN5dx = (yz1415*(-t+s*t)+yz1512*(-t+t*r)+yz1214*(1-s-r+r*s)         )*DetJ_inv   !
            dN5dy = (xz1514*(-t+s*t)+xz1215*(-t+t*r)+xz1214*(1-s-r+r*s)         )*DetJ_inv   !    Primeiras derivadas da funcao de forma
            dN5dz = (xy1514*(-t+s*t)+xy1215*(-t+t*r)+xy1214*(1-s-r+r*s)         )*DetJ_inv  !         avaliadas nos pontos de Gauss
                                                                                                                     !
            dN6dx = (yz1415*(t-s*t)-yz1512*t*r+yz1214*(r-r*s)                   )*DetJ_inv   !
            dN6dy = (xz1514*(t-s*t)-xz1215*t*r+xz1214*(r-r*s)                   )*DetJ_inv   !
            dN6dz = (xy1514*(t-s*t)-xy1215*t*r+xy1214*(r-r*s)                   )*DetJ_inv  !
                                                                                                                     !
            dN7dx = (yz1415*s*t+yz1512*t*r+yz1214*r*s                           )*DetJ_inv   !
            dN7dy = (xz1514*s*t+xz1215*t*r+xz1214*r*s                           )*DetJ_inv   !
            dN7dz = (xy1514*s*t+xy1215*t*r+xy1214*r*s                           )*DetJ_inv  !
                                                                                                                     !
            dN8dx = -(yz1415*s*t+yz1512*(t-t*r)+yz1214*(s-r*s)                  )*DetJ_inv   !
            dN8dy = -(xz1514*s*t+xz1215*(t-t*r)+xz1214*(s-r*s)                  )*DetJ_inv   !    Primeiras derivadas da funcao de forma
            dN8dz = -(xy1514*s*t+xy1215*(t-t*r)+xy1214*(s-r*s)                  )*DetJ_inv  !         avaliadas nos pontos de Gauss
                                                                                                                     !                                            !
             !!!!!!!!!!!!!!!! Matriz DIFUSIVA (Parcela Mu)


ske_difus_xMu (1) =  ske_difus_xMu (1) + dN1dx * dN1dx + dN1dy * dN1dy + dN1dz * dN1dz
ske_difus_xMu (26) =  ske_difus_xMu (26) + dN1dx * dN1dx + dN1dy * dN1dy + dN1dz * dN1dz
ske_difus_xMu (51) =  ske_difus_xMu (51) + dN1dx * dN1dx + dN1dy * dN1dy + dN1dz * dN1dz
ske_difus_xMu (73) =  ske_difus_xMu (73) + dN2dx * dN1dx + dN2dy * dN1dy + dN2dz * dN1dz
ske_difus_xMu (76) =  ske_difus_xMu (76) + dN2dx * dN2dx + dN2dy * dN2dy + dN2dz * dN2dz
ske_difus_xMu (98) =  ske_difus_xMu (98) + dN2dx * dN1dx + dN2dy * dN1dy + dN2dz * dN1dz
ske_difus_xMu (101) =  ske_difus_xMu (101) + dN2dx * dN2dx + dN2dy * dN2dy + dN2dz * dN2dz
ske_difus_xMu (123) =  ske_difus_xMu (123) + dN2dx * dN1dx + dN2dy * dN1dy + dN2dz * dN1dz
ske_difus_xMu (126) =  ske_difus_xMu (126) + dN2dx * dN2dx + dN2dy * dN2dy + dN2dz * dN2dz
ske_difus_xMu (145) =  ske_difus_xMu (145) + dN3dx * dN1dx + dN3dy * dN1dy + dN3dz * dN1dz
ske_difus_xMu (148) =  ske_difus_xMu (148) + dN3dx * dN2dx + dN3dy * dN2dy + dN3dz * dN2dz
ske_difus_xMu (151) =  ske_difus_xMu (151) + dN3dx * dN3dx + dN3dy * dN3dy + dN3dz * dN3dz
ske_difus_xMu (170) =  ske_difus_xMu (170) + dN3dx * dN1dx + dN3dy * dN1dy + dN3dz * dN1dz
ske_difus_xMu (173) =  ske_difus_xMu (173) + dN3dx * dN2dx + dN3dy * dN2dy + dN3dz * dN2dz
ske_difus_xMu (176) =  ske_difus_xMu (176) + dN3dx * dN3dx + dN3dy * dN3dy + dN3dz * dN3dz
ske_difus_xMu (195) =  ske_difus_xMu (195) + dN3dx * dN1dx + dN3dy * dN1dy + dN3dz * dN1dz
ske_difus_xMu (198) =  ske_difus_xMu (198) + dN3dx * dN2dx + dN3dy * dN2dy + dN3dz * dN2dz
ske_difus_xMu (201) =  ske_difus_xMu (201) + dN3dx * dN3dx + dN3dy * dN3dy + dN3dz * dN3dz
ske_difus_xMu (217) =  ske_difus_xMu (217) + dN4dx * dN1dx + dN4dy * dN1dy + dN4dz * dN1dz
ske_difus_xMu (220) =  ske_difus_xMu (220) + dN4dx * dN2dx + dN4dy * dN2dy + dN4dz * dN2dz
ske_difus_xMu (223) =  ske_difus_xMu (223) + dN4dx * dN3dx + dN4dy * dN3dy + dN4dz * dN3dz
ske_difus_xMu (226) =  ske_difus_xMu (226) + dN4dx * dN4dx + dN4dy * dN4dy + dN4dz * dN4dz
ske_difus_xMu (242) =  ske_difus_xMu (242) + dN4dx * dN1dx + dN4dy * dN1dy + dN4dz * dN1dz
ske_difus_xMu (245) =  ske_difus_xMu (245) + dN4dx * dN2dx + dN4dy * dN2dy + dN4dz * dN2dz
ske_difus_xMu (248) =  ske_difus_xMu (248) + dN4dx * dN3dx + dN4dy * dN3dy + dN4dz * dN3dz
ske_difus_xMu (251) =  ske_difus_xMu (251) + dN4dx * dN4dx + dN4dy * dN4dy + dN4dz * dN4dz
ske_difus_xMu (267) =  ske_difus_xMu (267) + dN4dx * dN1dx + dN4dy * dN1dy + dN4dz * dN1dz
ske_difus_xMu (270) =  ske_difus_xMu (270) + dN4dx * dN2dx + dN4dy * dN2dy + dN4dz * dN2dz
ske_difus_xMu (273) =  ske_difus_xMu (273) + dN4dx * dN3dx + dN4dy * dN3dy + dN4dz * dN3dz
ske_difus_xMu (276) =  ske_difus_xMu (276) + dN4dx * dN4dx + dN4dy * dN4dy + dN4dz * dN4dz
ske_difus_xMu (289) =  ske_difus_xMu (289) + dN5dx * dN1dx + dN5dy * dN1dy + dN5dz * dN1dz
ske_difus_xMu (292) =  ske_difus_xMu (292) + dN5dx * dN2dx + dN5dy * dN2dy + dN5dz * dN2dz
ske_difus_xMu (295) =  ske_difus_xMu (295) + dN5dx * dN3dx + dN5dy * dN3dy + dN5dz * dN3dz
ske_difus_xMu (298) =  ske_difus_xMu (298) + dN5dx * dN4dx + dN5dy * dN4dy + dN5dz * dN4dz
ske_difus_xMu (301) =  ske_difus_xMu (301) + dN5dx * dN5dx + dN5dy * dN5dy + dN5dz * dN5dz
ske_difus_xMu (314) =  ske_difus_xMu (314) + dN5dx * dN1dx + dN5dy * dN1dy + dN5dz * dN1dz
ske_difus_xMu (317) =  ske_difus_xMu (317) + dN5dx * dN2dx + dN5dy * dN2dy + dN5dz * dN2dz
ske_difus_xMu (320) =  ske_difus_xMu (320) + dN5dx * dN3dx + dN5dy * dN3dy + dN5dz * dN3dz
ske_difus_xMu (323) =  ske_difus_xMu (323) + dN5dx * dN4dx + dN5dy * dN4dy + dN5dz * dN4dz
ske_difus_xMu (326) =  ske_difus_xMu (326) + dN5dx * dN5dx + dN5dy * dN5dy + dN5dz * dN5dz
ske_difus_xMu (339) =  ske_difus_xMu (339) + dN5dx * dN1dx + dN5dy * dN1dy + dN5dz * dN1dz
ske_difus_xMu (342) =  ske_difus_xMu (342) + dN5dx * dN2dx + dN5dy * dN2dy + dN5dz * dN2dz
ske_difus_xMu (345) =  ske_difus_xMu (345) + dN5dx * dN3dx + dN5dy * dN3dy + dN5dz * dN3dz
ske_difus_xMu (348) =  ske_difus_xMu (348) + dN5dx * dN4dx + dN5dy * dN4dy + dN5dz * dN4dz
ske_difus_xMu (351) =  ske_difus_xMu (351) + dN5dx * dN5dx + dN5dy * dN5dy + dN5dz * dN5dz
ske_difus_xMu (361) =  ske_difus_xMu (361) + dN6dx * dN1dx + dN6dy * dN1dy + dN6dz * dN1dz
ske_difus_xMu (364) =  ske_difus_xMu (364) + dN6dx * dN2dx + dN6dy * dN2dy + dN6dz * dN2dz
ske_difus_xMu (367) =  ske_difus_xMu (367) + dN6dx * dN3dx + dN6dy * dN3dy + dN6dz * dN3dz
ske_difus_xMu (370) =  ske_difus_xMu (370) + dN6dx * dN4dx + dN6dy * dN4dy + dN6dz * dN4dz
ske_difus_xMu (373) =  ske_difus_xMu (373) + dN6dx * dN5dx + dN6dy * dN5dy + dN6dz * dN5dz
ske_difus_xMu (376) =  ske_difus_xMu (376) + dN6dx * dN6dx + dN6dy * dN6dy + dN6dz * dN6dz
ske_difus_xMu (386) =  ske_difus_xMu (386) + dN6dx * dN1dx + dN6dy * dN1dy + dN6dz * dN1dz
ske_difus_xMu (389) =  ske_difus_xMu (389) + dN6dx * dN2dx + dN6dy * dN2dy + dN6dz * dN2dz
ske_difus_xMu (392) =  ske_difus_xMu (392) + dN6dx * dN3dx + dN6dy * dN3dy + dN6dz * dN3dz
ske_difus_xMu (395) =  ske_difus_xMu (395) + dN6dx * dN4dx + dN6dy * dN4dy + dN6dz * dN4dz
ske_difus_xMu (398) =  ske_difus_xMu (398) + dN6dx * dN5dx + dN6dy * dN5dy + dN6dz * dN5dz
ske_difus_xMu (401) =  ske_difus_xMu (401) + dN6dx * dN6dx + dN6dy * dN6dy + dN6dz * dN6dz
ske_difus_xMu (411) =  ske_difus_xMu (411) + dN6dx * dN1dx + dN6dy * dN1dy + dN6dz * dN1dz
ske_difus_xMu (414) =  ske_difus_xMu (414) + dN6dx * dN2dx + dN6dy * dN2dy + dN6dz * dN2dz
ske_difus_xMu (417) =  ske_difus_xMu (417) + dN6dx * dN3dx + dN6dy * dN3dy + dN6dz * dN3dz
ske_difus_xMu (420) =  ske_difus_xMu (420) + dN6dx * dN4dx + dN6dy * dN4dy + dN6dz * dN4dz
ske_difus_xMu (423) =  ske_difus_xMu (423) + dN6dx * dN5dx + dN6dy * dN5dy + dN6dz * dN5dz
ske_difus_xMu (426) =  ske_difus_xMu (426) + dN6dx * dN6dx + dN6dy * dN6dy + dN6dz * dN6dz
ske_difus_xMu (433) =  ske_difus_xMu (433) + dN7dx * dN1dx + dN7dy * dN1dy + dN7dz * dN1dz
ske_difus_xMu (436) =  ske_difus_xMu (436) + dN7dx * dN2dx + dN7dy * dN2dy + dN7dz * dN2dz
ske_difus_xMu (439) =  ske_difus_xMu (439) + dN7dx * dN3dx + dN7dy * dN3dy + dN7dz * dN3dz
ske_difus_xMu (442) =  ske_difus_xMu (442) + dN7dx * dN4dx + dN7dy * dN4dy + dN7dz * dN4dz
ske_difus_xMu (445) =  ske_difus_xMu (445) + dN7dx * dN5dx + dN7dy * dN5dy + dN7dz * dN5dz
ske_difus_xMu (448) =  ske_difus_xMu (448) + dN7dx * dN6dx + dN7dy * dN6dy + dN7dz * dN6dz
ske_difus_xMu (451) =  ske_difus_xMu (451) + dN7dx * dN7dx + dN7dy * dN7dy + dN7dz * dN7dz
ske_difus_xMu (458) =  ske_difus_xMu (458) + dN7dx * dN1dx + dN7dy * dN1dy + dN7dz * dN1dz
ske_difus_xMu (461) =  ske_difus_xMu (461) + dN7dx * dN2dx + dN7dy * dN2dy + dN7dz * dN2dz
ske_difus_xMu (464) =  ske_difus_xMu (464) + dN7dx * dN3dx + dN7dy * dN3dy + dN7dz * dN3dz
ske_difus_xMu (467) =  ske_difus_xMu (467) + dN7dx * dN4dx + dN7dy * dN4dy + dN7dz * dN4dz
ske_difus_xMu (470) =  ske_difus_xMu (470) + dN7dx * dN5dx + dN7dy * dN5dy + dN7dz * dN5dz
ske_difus_xMu (473) =  ske_difus_xMu (473) + dN7dx * dN6dx + dN7dy * dN6dy + dN7dz * dN6dz
ske_difus_xMu (476) =  ske_difus_xMu (476) + dN7dx * dN7dx + dN7dy * dN7dy + dN7dz * dN7dz
ske_difus_xMu (483) =  ske_difus_xMu (483) + dN7dx * dN1dx + dN7dy * dN1dy + dN7dz * dN1dz
ske_difus_xMu (486) =  ske_difus_xMu (486) + dN7dx * dN2dx + dN7dy * dN2dy + dN7dz * dN2dz
ske_difus_xMu (489) =  ske_difus_xMu (489) + dN7dx * dN3dx + dN7dy * dN3dy + dN7dz * dN3dz
ske_difus_xMu (492) =  ske_difus_xMu (492) + dN7dx * dN4dx + dN7dy * dN4dy + dN7dz * dN4dz
ske_difus_xMu (495) =  ske_difus_xMu (495) + dN7dx * dN5dx + dN7dy * dN5dy + dN7dz * dN5dz
ske_difus_xMu (498) =  ske_difus_xMu (498) + dN7dx * dN6dx + dN7dy * dN6dy + dN7dz * dN6dz
ske_difus_xMu (501) =  ske_difus_xMu (501) + dN7dx * dN7dx + dN7dy * dN7dy + dN7dz * dN7dz
ske_difus_xMu (505) =  ske_difus_xMu (505) + dN8dx * dN1dx + dN8dy * dN1dy + dN8dz * dN1dz
ske_difus_xMu (508) =  ske_difus_xMu (508) + dN8dx * dN2dx + dN8dy * dN2dy + dN8dz * dN2dz
ske_difus_xMu (511) =  ske_difus_xMu (511) + dN8dx * dN3dx + dN8dy * dN3dy + dN8dz * dN3dz
ske_difus_xMu (514) =  ske_difus_xMu (514) + dN8dx * dN4dx + dN8dy * dN4dy + dN8dz * dN4dz
ske_difus_xMu (517) =  ske_difus_xMu (517) + dN8dx * dN5dx + dN8dy * dN5dy + dN8dz * dN5dz
ske_difus_xMu (520) =  ske_difus_xMu (520) + dN8dx * dN6dx + dN8dy * dN6dy + dN8dz * dN6dz
ske_difus_xMu (523) =  ske_difus_xMu (523) + dN8dx * dN7dx + dN8dy * dN7dy + dN8dz * dN7dz
ske_difus_xMu (526) =  ske_difus_xMu (526) + dN8dx * dN8dx + dN8dy * dN8dy + dN8dz * dN8dz
ske_difus_xMu (530) =  ske_difus_xMu (530) + dN8dx * dN1dx + dN8dy * dN1dy + dN8dz * dN1dz
ske_difus_xMu (533) =  ske_difus_xMu (533) + dN8dx * dN2dx + dN8dy * dN2dy + dN8dz * dN2dz
ske_difus_xMu (536) =  ske_difus_xMu (536) + dN8dx * dN3dx + dN8dy * dN3dy + dN8dz * dN3dz
ske_difus_xMu (539) =  ske_difus_xMu (539) + dN8dx * dN4dx + dN8dy * dN4dy + dN8dz * dN4dz
ske_difus_xMu (542) =  ske_difus_xMu (542) + dN8dx * dN5dx + dN8dy * dN5dy + dN8dz * dN5dz
ske_difus_xMu (545) =  ske_difus_xMu (545) + dN8dx * dN6dx + dN8dy * dN6dy + dN8dz * dN6dz
ske_difus_xMu (548) =  ske_difus_xMu (548) + dN8dx * dN7dx + dN8dy * dN7dy + dN8dz * dN7dz
ske_difus_xMu (551) =  ske_difus_xMu (551) + dN8dx * dN8dx + dN8dy * dN8dy + dN8dz * dN8dz
ske_difus_xMu (555) =  ske_difus_xMu (555) + dN8dx * dN1dx + dN8dy * dN1dy + dN8dz * dN1dz
ske_difus_xMu (558) =  ske_difus_xMu (558) + dN8dx * dN2dx + dN8dy * dN2dy + dN8dz * dN2dz
ske_difus_xMu (561) =  ske_difus_xMu (561) + dN8dx * dN3dx + dN8dy * dN3dy + dN8dz * dN3dz
ske_difus_xMu (564) =  ske_difus_xMu (564) + dN8dx * dN4dx + dN8dy * dN4dy + dN8dz * dN4dz
ske_difus_xMu (567) =  ske_difus_xMu (567) + dN8dx * dN5dx + dN8dy * dN5dy + dN8dz * dN5dz
ske_difus_xMu (570) =  ske_difus_xMu (570) + dN8dx * dN6dx + dN8dy * dN6dy + dN8dz * dN6dz
ske_difus_xMu (573) =  ske_difus_xMu (573) + dN8dx * dN7dx + dN8dy * dN7dy + dN8dz * dN7dz
ske_difus_xMu (576) =  ske_difus_xMu (576) + dN8dx * dN8dx + dN8dy * dN8dy + dN8dz * dN8dz
                                                                                                              !
             
      enddo !Fim do loop nos pontos de Gauss FULL INTEGRATION

           ! Begin: Matriz DIFUSIVA (Parcela LAMBDA) integra��o reduzida

             r = 0.5d0
             s = 0.5d0
             t = 0.5d0

            dN1dx = (yz1415*(-1+s+t-s*t)+yz1512*(-1+r+t-t*r)+yz1214*(-1+s+r-r*s))*DetJ_inv   !
            dN1dy = (xz1514*(-1+s+t-s*t)+xz1215*(-1+r+t-t*r)+xz1214*(-1+s+r-r*s))*DetJ_inv  !
            dN1dz = (xy1514*(-1+s+t-s*t)+xy1215*(-1+r+t-t*r)+xy1214*(-1+s+r-r*s))*DetJ_inv  !
                                                                                                                     !
            dN2dx = (yz1415*(1-s-t+s*t)+yz1512*(-r+t*r)+yz1214*(-r+r*s)         )*DetJ_inv  !
            dN2dy = (xz1514*(1-s-t+s*t)+xz1215*(-r+t*r)+xz1214*(-r+r*s)         )*DetJ_inv   !    Primeiras derivadas da funcao de forma
            dN2dz = (xy1514*(1-s-t+s*t)+xy1215*(-r+t*r)+xy1214*(-r+r*s)         )*DetJ_inv  !         avaliadas nos pontos de Gauss
                                                                                                                     !
            dN3dx = (yz1415*(s-s*t)+yz1512*(r-t*r)-yz1214*r*s                   )*DetJ_inv   !
            dN3dy = (xz1514*(s-s*t)+xz1215*(r-t*r)-xz1214*r*s                   )*DetJ_inv   !
            dN3dz = (xy1514*(s-s*t)+xy1215*(r-t*r)-xy1214*r*s                   )*DetJ_inv  !
                                                                                                                    !
            dN4dx = (yz1415*(-s+s*t)+yz1512*(1-r-t+t*r)+yz1214*(-s+r*s)         )*DetJ_inv   !
            dN4dy = (xz1514*(-s+s*t)+xz1215*(1-r-t+t*r)+xz1214*(-s+r*s)         )*DetJ_inv   !
            dN4dz = (xy1514*(-s+s*t)+xy1215*(1-r-t+t*r)+xy1214*(-s+r*s)         )*DetJ_inv  !
                                                                                                                     !
            dN5dx = (yz1415*(-t+s*t)+yz1512*(-t+t*r)+yz1214*(1-s-r+r*s)         )*DetJ_inv   !
            dN5dy = (xz1514*(-t+s*t)+xz1215*(-t+t*r)+xz1214*(1-s-r+r*s)         )*DetJ_inv   !    Primeiras derivadas da funcao de forma
            dN5dz = (xy1514*(-t+s*t)+xy1215*(-t+t*r)+xy1214*(1-s-r+r*s)         )*DetJ_inv  !         avaliadas nos pontos de Gauss
                                                                                                                     !
            dN6dx = (yz1415*(t-s*t)-yz1512*t*r+yz1214*(r-r*s)                   )*DetJ_inv   !
            dN6dy = (xz1514*(t-s*t)-xz1215*t*r+xz1214*(r-r*s)                   )*DetJ_inv   !
            dN6dz = (xy1514*(t-s*t)-xy1215*t*r+xy1214*(r-r*s)                   )*DetJ_inv  !
                                                                                                                     !
            dN7dx = (yz1415*s*t+yz1512*t*r+yz1214*r*s                           )*DetJ_inv   !
            dN7dy = (xz1514*s*t+xz1215*t*r+xz1214*r*s                           )*DetJ_inv   !
            dN7dz = (xy1514*s*t+xy1215*t*r+xy1214*r*s                           )*DetJ_inv  !
                                                                                                                     !
            dN8dx = -(yz1415*s*t+yz1512*(t-t*r)+yz1214*(s-r*s)                  )*DetJ_inv   !
            dN8dy = -(xz1514*s*t+xz1215*(t-t*r)+xz1214*(s-r*s)                  )*DetJ_inv   !    Primeiras derivadas da funcao de forma
            dN8dz = -(xy1514*s*t+xy1215*(t-t*r)+xy1214*(s-r*s)                  )*DetJ_inv  !         avaliadas nos pontos de Gauss
                                                                                                                     !
       
            ske_difus_Lambda (1) =  ske_difus_Lambda (1) + dN1dx * dN1dx
            ske_difus_Lambda (25) =  ske_difus_Lambda (25) + dN1dy * dN1dx           !
            ske_difus_Lambda (26) =  ske_difus_Lambda (26) + dN1dy * dN1dy            !
            ske_difus_Lambda (49) =  ske_difus_Lambda (49) + dN1dz * dN1dx            !
            ske_difus_Lambda (50) =  ske_difus_Lambda (50) + dN1dz * dN1dy            !
            ske_difus_Lambda (51) =  ske_difus_Lambda (51) + dN1dz * dN1dz            ! MATRIZ DIFUSIVA
            ske_difus_Lambda (73) =  ske_difus_Lambda (73) + dN2dx * dN1dx            !
            ske_difus_Lambda (74) =  ske_difus_Lambda (74) + dN2dx * dN1dy            !    PARCELA LAMBDA
            ske_difus_Lambda (75) =  ske_difus_Lambda (75) + dN2dx * dN1dz            !
            ske_difus_Lambda (76) =  ske_difus_Lambda (76) + dN2dx * dN2dx            ! REDUCED INTEGRATION
            ske_difus_Lambda (97) =  ske_difus_Lambda (97) + dN2dy * dN1dx            !
            ske_difus_Lambda (98) =  ske_difus_Lambda (98) + dN2dy * dN1dy            !
            ske_difus_Lambda (99) =  ske_difus_Lambda (99) + dN2dy * dN1dz            !
            ske_difus_Lambda (100) =  ske_difus_Lambda (100) + dN2dy * dN2dx          !
            ske_difus_Lambda (101) =  ske_difus_Lambda (101) + dN2dy * dN2dy          !
            ske_difus_Lambda (121) =  ske_difus_Lambda (121) + dN2dz * dN1dx          !
            ske_difus_Lambda (122) =  ske_difus_Lambda (122) + dN2dz * dN1dy          !
            ske_difus_Lambda (123) =  ske_difus_Lambda (123) + dN2dz * dN1dz          !
            ske_difus_Lambda (124) =  ske_difus_Lambda (124) + dN2dz * dN2dx          !
            ske_difus_Lambda (125) =  ske_difus_Lambda (125) + dN2dz * dN2dy          !
            ske_difus_Lambda (126) =  ske_difus_Lambda (126) + dN2dz * dN2dz          !
            ske_difus_Lambda (145) =  ske_difus_Lambda (145) + dN3dx * dN1dx          !
            ske_difus_Lambda (146) =  ske_difus_Lambda (146) + dN3dx * dN1dy          !
            ske_difus_Lambda (147) =  ske_difus_Lambda (147) + dN3dx * dN1dz          !
            ske_difus_Lambda (148) =  ske_difus_Lambda (148) + dN3dx * dN2dx          !
            ske_difus_Lambda (149) =  ske_difus_Lambda (149) + dN3dx * dN2dy          !
            ske_difus_Lambda (150) =  ske_difus_Lambda (150) + dN3dx * dN2dz          ! MATRIZ DIFUSIVA
            ske_difus_Lambda (151) =  ske_difus_Lambda (151) + dN3dx * dN3dx          !
            ske_difus_Lambda (169) =  ske_difus_Lambda (169) + dN3dy * dN1dx          !    PARCELA LAMBDA
            ske_difus_Lambda (170) =  ske_difus_Lambda (170) + dN3dy * dN1dy          !
            ske_difus_Lambda (171) =  ske_difus_Lambda (171) + dN3dy * dN1dz          ! REDUCED INTEGRATION
            ske_difus_Lambda (172) =  ske_difus_Lambda (172) + dN3dy * dN2dx          !
            ske_difus_Lambda (173) =  ske_difus_Lambda (173) + dN3dy * dN2dy          !
            ske_difus_Lambda (174) =  ske_difus_Lambda (174) + dN3dy * dN2dz          !
            ske_difus_Lambda (175) =  ske_difus_Lambda (175) + dN3dy * dN3dx          !
            ske_difus_Lambda (176) =  ske_difus_Lambda (176) + dN3dy * dN3dy          !
            ske_difus_Lambda (193) =  ske_difus_Lambda (193) + dN3dz * dN1dx          !
            ske_difus_Lambda (194) =  ske_difus_Lambda (194) + dN3dz * dN1dy          !
            ske_difus_Lambda (195) =  ske_difus_Lambda (195) + dN3dz * dN1dz          !
            ske_difus_Lambda (196) =  ske_difus_Lambda (196) + dN3dz * dN2dx          !
            ske_difus_Lambda (197) =  ske_difus_Lambda (197) + dN3dz * dN2dy          !
            ske_difus_Lambda (198) =  ske_difus_Lambda (198) + dN3dz * dN2dz          !
            ske_difus_Lambda (199) =  ske_difus_Lambda (199) + dN3dz * dN3dx          !
            ske_difus_Lambda (200) =  ske_difus_Lambda (200) + dN3dz * dN3dy          !
            ske_difus_Lambda (201) =  ske_difus_Lambda (201) + dN3dz * dN3dz          !
            ske_difus_Lambda (217) =  ske_difus_Lambda (217) + dN4dx * dN1dx          !
            ske_difus_Lambda (218) =  ske_difus_Lambda (218) + dN4dx * dN1dy          !
            ske_difus_Lambda (219) =  ske_difus_Lambda (219) + dN4dx * dN1dz          ! MATRIZ DIFUSIVA
            ske_difus_Lambda (220) =  ske_difus_Lambda (220) + dN4dx * dN2dx          !
            ske_difus_Lambda (221) =  ske_difus_Lambda (221) + dN4dx * dN2dy          !    PARCELA LAMBDA
            ske_difus_Lambda (222) =  ske_difus_Lambda (222) + dN4dx * dN2dz          !
            ske_difus_Lambda (223) =  ske_difus_Lambda (223) + dN4dx * dN3dx          ! REDUCED INTEGRATION
            ske_difus_Lambda (224) =  ske_difus_Lambda (224) + dN4dx * dN3dy          !
            ske_difus_Lambda (225) =  ske_difus_Lambda (225) + dN4dx * dN3dz          !
            ske_difus_Lambda (226) =  ske_difus_Lambda (226) + dN4dx * dN4dx          !
            ske_difus_Lambda (241) =  ske_difus_Lambda (241) + dN4dy * dN1dx          !
            ske_difus_Lambda (242) =  ske_difus_Lambda (242) + dN4dy * dN1dy          !
            ske_difus_Lambda (243) =  ske_difus_Lambda (243) + dN4dy * dN1dz          !
            ske_difus_Lambda (244) =  ske_difus_Lambda (244) + dN4dy * dN2dx          !
            ske_difus_Lambda (245) =  ske_difus_Lambda (245) + dN4dy * dN2dy          !
            ske_difus_Lambda (246) =  ske_difus_Lambda (246) + dN4dy * dN2dz          !
            ske_difus_Lambda (247) =  ske_difus_Lambda (247) + dN4dy * dN3dx          !
            ske_difus_Lambda (248) =  ske_difus_Lambda (248) + dN4dy * dN3dy          !
            ske_difus_Lambda (249) =  ske_difus_Lambda (249) + dN4dy * dN3dz          !
            ske_difus_Lambda (250) =  ske_difus_Lambda (250) + dN4dy * dN4dx          !
            ske_difus_Lambda (251) =  ske_difus_Lambda (251) + dN4dy * dN4dy          !
            ske_difus_Lambda (265) =  ske_difus_Lambda (265) + dN4dz * dN1dx          !
            ske_difus_Lambda (266) =  ske_difus_Lambda (266) + dN4dz * dN1dy          !
            ske_difus_Lambda (267) =  ske_difus_Lambda (267) + dN4dz * dN1dz          ! MATRIZ DIFUSIVA
            ske_difus_Lambda (268) =  ske_difus_Lambda (268) + dN4dz * dN2dx          !
            ske_difus_Lambda (269) =  ske_difus_Lambda (269) + dN4dz * dN2dy          !    PARCELA LAMBDA
            ske_difus_Lambda (270) =  ske_difus_Lambda (270) + dN4dz * dN2dz          !
            ske_difus_Lambda (271) =  ske_difus_Lambda (271) + dN4dz * dN3dx          ! REDUCED INTEGRATION
            ske_difus_Lambda (272) =  ske_difus_Lambda (272) + dN4dz * dN3dy          !
            ske_difus_Lambda (273) =  ske_difus_Lambda (273) + dN4dz * dN3dz          !
            ske_difus_Lambda (274) =  ske_difus_Lambda (274) + dN4dz * dN4dx          !
            ske_difus_Lambda (275) =  ske_difus_Lambda (275) + dN4dz * dN4dy          !
            ske_difus_Lambda (276) =  ske_difus_Lambda (276) + dN4dz * dN4dz          !
            ske_difus_Lambda (289) =  ske_difus_Lambda (289) + dN5dx * dN1dx
            ske_difus_Lambda (290) =  ske_difus_Lambda (290) + dN5dx * dN1dy          !
            ske_difus_Lambda (291) =  ske_difus_Lambda (291) + dN5dx * dN1dz          !
            ske_difus_Lambda (292) =  ske_difus_Lambda (292) + dN5dx * dN2dx          !
            ske_difus_Lambda (293) =  ske_difus_Lambda (293) + dN5dx * dN2dy          !
            ske_difus_Lambda (294) =  ske_difus_Lambda (294) + dN5dx * dN2dz          ! MATRIZ DIFUSIVA
            ske_difus_Lambda (295) =  ske_difus_Lambda (295) + dN5dx * dN3dx          !
            ske_difus_Lambda (296) =  ske_difus_Lambda (296) + dN5dx * dN3dy          !    PARCELA LAMBDA
            ske_difus_Lambda (297) =  ske_difus_Lambda (297) + dN5dx * dN3dz          !
            ske_difus_Lambda (298) =  ske_difus_Lambda (298) + dN5dx * dN4dx          ! REDUCED INTEGRATION
            ske_difus_Lambda (299) =  ske_difus_Lambda (299) + dN5dx * dN4dy          !
            ske_difus_Lambda (300) =  ske_difus_Lambda (300) + dN5dx * dN4dz          !
            ske_difus_Lambda (301) =  ske_difus_Lambda (301) + dN5dx * dN5dx          !
            ske_difus_Lambda (313) =  ske_difus_Lambda (313) + dN5dy * dN1dx          !
            ske_difus_Lambda (314) =  ske_difus_Lambda (314) + dN5dy * dN1dy          !
            ske_difus_Lambda (315) =  ske_difus_Lambda (315) + dN5dy * dN1dz          !
            ske_difus_Lambda (316) =  ske_difus_Lambda (316) + dN5dy * dN2dx          !
            ske_difus_Lambda (317) =  ske_difus_Lambda (317) + dN5dy * dN2dy          !
            ske_difus_Lambda (318) =  ske_difus_Lambda (318) + dN5dy * dN2dz          !
            ske_difus_Lambda (319) =  ske_difus_Lambda (319) + dN5dy * dN3dx          !
            ske_difus_Lambda (320) =  ske_difus_Lambda (320) + dN5dy * dN3dy          !
            ske_difus_Lambda (321) =  ske_difus_Lambda (321) + dN5dy * dN3dz          !
            ske_difus_Lambda (322) =  ske_difus_Lambda (322) + dN5dy * dN4dx          !
            ske_difus_Lambda (323) =  ske_difus_Lambda (323) + dN5dy * dN4dy          !
            ske_difus_Lambda (324) =  ske_difus_Lambda (324) + dN5dy * dN4dz          !
            ske_difus_Lambda (325) =  ske_difus_Lambda (325) + dN5dy * dN5dx          !
            ske_difus_Lambda (326) =  ske_difus_Lambda (326) + dN5dy * dN5dy          ! MATRIZ DIFUSIVA
            ske_difus_Lambda (337) =  ske_difus_Lambda (337) + dN5dz * dN1dx          !
            ske_difus_Lambda (338) =  ske_difus_Lambda (338) + dN5dz * dN1dy          !    PARCELA LAMBDA
            ske_difus_Lambda (339) =  ske_difus_Lambda (339) + dN5dz * dN1dz          !
            ske_difus_Lambda (340) =  ske_difus_Lambda (340) + dN5dz * dN2dx          ! REDUCED INTEGRATION
            ske_difus_Lambda (341) =  ske_difus_Lambda (341) + dN5dz * dN2dy          !
            ske_difus_Lambda (342) =  ske_difus_Lambda (342) + dN5dz * dN2dz          !
            ske_difus_Lambda (343) =  ske_difus_Lambda (343) + dN5dz * dN3dx          !
            ske_difus_Lambda (344) =  ske_difus_Lambda (344) + dN5dz * dN3dy          !
            ske_difus_Lambda (345) =  ske_difus_Lambda (345) + dN5dz * dN3dz          !
            ske_difus_Lambda (346) =  ske_difus_Lambda (346) + dN5dz * dN4dx          !
            ske_difus_Lambda (347) =  ske_difus_Lambda (347) + dN5dz * dN4dy          !
            ske_difus_Lambda (348) =  ske_difus_Lambda (348) + dN5dz * dN4dz          !
            ske_difus_Lambda (349) =  ske_difus_Lambda (349) + dN5dz * dN5dx          !
            ske_difus_Lambda (350) =  ske_difus_Lambda (350) + dN5dz * dN5dy          !
            ske_difus_Lambda (351) =  ske_difus_Lambda (351) + dN5dz * dN5dz          !
            ske_difus_Lambda (361) =  ske_difus_Lambda (361) + dN6dx * dN1dx          !
            ske_difus_Lambda (362) =  ske_difus_Lambda (362) + dN6dx * dN1dy          !
            ske_difus_Lambda (363) =  ske_difus_Lambda (363) + dN6dx * dN1dz          !
            ske_difus_Lambda (364) =  ske_difus_Lambda (364) + dN6dx * dN2dx          !
            ske_difus_Lambda (365) =  ske_difus_Lambda (365) + dN6dx * dN2dy          !
            ske_difus_Lambda (366) =  ske_difus_Lambda (366) + dN6dx * dN2dz          ! MATRIZ DIFUSIVA
            ske_difus_Lambda (367) =  ske_difus_Lambda (367) + dN6dx * dN3dx          !
            ske_difus_Lambda (368) =  ske_difus_Lambda (368) + dN6dx * dN3dy          !    PARCELA LAMBDA
            ske_difus_Lambda (369) =  ske_difus_Lambda (369) + dN6dx * dN3dz          !
            ske_difus_Lambda (370) =  ske_difus_Lambda (370) + dN6dx * dN4dx          ! REDUCED INTEGRATION
            ske_difus_Lambda (371) =  ske_difus_Lambda (371) + dN6dx * dN4dy          !
            ske_difus_Lambda (372) =  ske_difus_Lambda (372) + dN6dx * dN4dz          !
            ske_difus_Lambda (373) =  ske_difus_Lambda (373) + dN6dx * dN5dx          !
            ske_difus_Lambda (374) =  ske_difus_Lambda (374) + dN6dx * dN5dy          !
            ske_difus_Lambda (375) =  ske_difus_Lambda (375) + dN6dx * dN5dz          !
            ske_difus_Lambda (376) =  ske_difus_Lambda (376) + dN6dx * dN6dx          !
            ske_difus_Lambda (385) =  ske_difus_Lambda (385) + dN6dy * dN1dx          !
            ske_difus_Lambda (386) =  ske_difus_Lambda (386) + dN6dy * dN1dy          !
            ske_difus_Lambda (387) =  ske_difus_Lambda (387) + dN6dy * dN1dz          !
            ske_difus_Lambda (388) =  ske_difus_Lambda (388) + dN6dy * dN2dx          !
            ske_difus_Lambda (389) =  ske_difus_Lambda (389) + dN6dy * dN2dy          !
            ske_difus_Lambda (390) =  ske_difus_Lambda (390) + dN6dy * dN2dz          !
            ske_difus_Lambda (391) =  ske_difus_Lambda (391) + dN6dy * dN3dx          !
            ske_difus_Lambda (392) =  ske_difus_Lambda (392) + dN6dy * dN3dy          !
            ske_difus_Lambda (393) =  ske_difus_Lambda (393) + dN6dy * dN3dz          !
            ske_difus_Lambda (394) =  ske_difus_Lambda (394) + dN6dy * dN4dx          !
            ske_difus_Lambda (395) =  ske_difus_Lambda (395) + dN6dy * dN4dy          ! MATRIZ DIFUSIVA
            ske_difus_Lambda (396) =  ske_difus_Lambda (396) + dN6dy * dN4dz          !
            ske_difus_Lambda (397) =  ske_difus_Lambda (397) + dN6dy * dN5dx          !    PARCELA LAMBDA
            ske_difus_Lambda (398) =  ske_difus_Lambda (398) + dN6dy * dN5dy          !
            ske_difus_Lambda (399) =  ske_difus_Lambda (399) + dN6dy * dN5dz          ! REDUCED INTEGRATION
            ske_difus_Lambda (400) =  ske_difus_Lambda (400) + dN6dy * dN6dx          !
            ske_difus_Lambda (401) =  ske_difus_Lambda (401) + dN6dy * dN6dy          !
            ske_difus_Lambda (409) =  ske_difus_Lambda (409) + dN6dz * dN1dx          !
            ske_difus_Lambda (410) =  ske_difus_Lambda (410) + dN6dz * dN1dy          !
            ske_difus_Lambda (411) =  ske_difus_Lambda (411) + dN6dz * dN1dz          !
            ske_difus_Lambda (412) =  ske_difus_Lambda (412) + dN6dz * dN2dx
            ske_difus_Lambda (413) =  ske_difus_Lambda (413) + dN6dz * dN2dy          !
            ske_difus_Lambda (414) =  ske_difus_Lambda (414) + dN6dz * dN2dz          !
            ske_difus_Lambda (415) =  ske_difus_Lambda (415) + dN6dz * dN3dx          !
            ske_difus_Lambda (416) =  ske_difus_Lambda (416) + dN6dz * dN3dy          !
            ske_difus_Lambda (417) =  ske_difus_Lambda (417) + dN6dz * dN3dz          ! MATRIZ DIFUSIVA
            ske_difus_Lambda (418) =  ske_difus_Lambda (418) + dN6dz * dN4dx          !
            ske_difus_Lambda (419) =  ske_difus_Lambda (419) + dN6dz * dN4dy          !    PARCELA LAMBDA
            ske_difus_Lambda (420) =  ske_difus_Lambda (420) + dN6dz * dN4dz          !
            ske_difus_Lambda (421) =  ske_difus_Lambda (421) + dN6dz * dN5dx          ! REDUCED INTEGRATION
            ske_difus_Lambda (422) =  ske_difus_Lambda (422) + dN6dz * dN5dy          !
            ske_difus_Lambda (423) =  ske_difus_Lambda (423) + dN6dz * dN5dz          !
            ske_difus_Lambda (424) =  ske_difus_Lambda (424) + dN6dz * dN6dx          !
            ske_difus_Lambda (425) =  ske_difus_Lambda (425) + dN6dz * dN6dy          !
            ske_difus_Lambda (426) =  ske_difus_Lambda (426) + dN6dz * dN6dz          !
            ske_difus_Lambda (433) =  ske_difus_Lambda (433) + dN7dx * dN1dx          !
            ske_difus_Lambda (434) =  ske_difus_Lambda (434) + dN7dx * dN1dy          !
            ske_difus_Lambda (435) =  ske_difus_Lambda (435) + dN7dx * dN1dz          !
            ske_difus_Lambda (436) =  ske_difus_Lambda (436) + dN7dx * dN2dx          !
            ske_difus_Lambda (437) =  ske_difus_Lambda (437) + dN7dx * dN2dy          !
            ske_difus_Lambda (438) =  ske_difus_Lambda (438) + dN7dx * dN2dz          !
            ske_difus_Lambda (439) =  ske_difus_Lambda (439) + dN7dx * dN3dx          !
            ske_difus_Lambda (440) =  ske_difus_Lambda (440) + dN7dx * dN3dy          !
            ske_difus_Lambda (441) =  ske_difus_Lambda (441) + dN7dx * dN3dz          !
            ske_difus_Lambda (442) =  ske_difus_Lambda (442) + dN7dx * dN4dx          !
            ske_difus_Lambda (443) =  ske_difus_Lambda (443) + dN7dx * dN4dy          !
            ske_difus_Lambda (444) =  ske_difus_Lambda (444) + dN7dx * dN4dz          ! MATRIZ DIFUSIVA
            ske_difus_Lambda (445) =  ske_difus_Lambda (445) + dN7dx * dN5dx          !
            ske_difus_Lambda (446) =  ske_difus_Lambda (446) + dN7dx * dN5dy          !    PARCELA LAMBDA
            ske_difus_Lambda (447) =  ske_difus_Lambda (447) + dN7dx * dN5dz          !
            ske_difus_Lambda (448) =  ske_difus_Lambda (448) + dN7dx * dN6dx          ! REDUCED INTEGRATION
            ske_difus_Lambda (449) =  ske_difus_Lambda (449) + dN7dx * dN6dy          !
            ske_difus_Lambda (450) =  ske_difus_Lambda (450) + dN7dx * dN6dz          !
            ske_difus_Lambda (451) =  ske_difus_Lambda (451) + dN7dx * dN7dx          !
            ske_difus_Lambda (457) =  ske_difus_Lambda (457) + dN7dy * dN1dx          !
            ske_difus_Lambda (458) =  ske_difus_Lambda (458) + dN7dy * dN1dy          !
            ske_difus_Lambda (459) =  ske_difus_Lambda (459) + dN7dy * dN1dz          !
            ske_difus_Lambda (460) =  ske_difus_Lambda (460) + dN7dy * dN2dx          !
            ske_difus_Lambda (461) =  ske_difus_Lambda (461) + dN7dy * dN2dy          !
            ske_difus_Lambda (462) =  ske_difus_Lambda (462) + dN7dy * dN2dz          !
            ske_difus_Lambda (463) =  ske_difus_Lambda (463) + dN7dy * dN3dx          !
            ske_difus_Lambda (464) =  ske_difus_Lambda (464) + dN7dy * dN3dy          !
            ske_difus_Lambda (465) =  ske_difus_Lambda (465) + dN7dy * dN3dz          !
            ske_difus_Lambda (466) =  ske_difus_Lambda (466) + dN7dy * dN4dx          !
            ske_difus_Lambda (467) =  ske_difus_Lambda (467) + dN7dy * dN4dy          !
            ske_difus_Lambda (468) =  ske_difus_Lambda (468) + dN7dy * dN4dz          !
            ske_difus_Lambda (469) =  ske_difus_Lambda (469) + dN7dy * dN5dx          !
            ske_difus_Lambda (470) =  ske_difus_Lambda (470) + dN7dy * dN5dy          ! MATRIZ DIFUSIVA
            ske_difus_Lambda (471) =  ske_difus_Lambda (471) + dN7dy * dN5dz          !
            ske_difus_Lambda (472) =  ske_difus_Lambda (472) + dN7dy * dN6dx          !    PARCELA LAMBDA
            ske_difus_Lambda (473) =  ske_difus_Lambda (473) + dN7dy * dN6dy          !
            ske_difus_Lambda (474) =  ske_difus_Lambda (474) + dN7dy * dN6dz          ! REDUCED INTEGRATION
            ske_difus_Lambda (475) =  ske_difus_Lambda (475) + dN7dy * dN7dx          !
            ske_difus_Lambda (476) =  ske_difus_Lambda (476) + dN7dy * dN7dy          !
            ske_difus_Lambda (481) =  ske_difus_Lambda (481) + dN7dz * dN1dx          !
            ske_difus_Lambda (482) =  ske_difus_Lambda (482) + dN7dz * dN1dy          !
            ske_difus_Lambda (483) =  ske_difus_Lambda (483) + dN7dz * dN1dz          !
            ske_difus_Lambda (484) =  ske_difus_Lambda (484) + dN7dz * dN2dx          !
            ske_difus_Lambda (485) =  ske_difus_Lambda (485) + dN7dz * dN2dy          !
            ske_difus_Lambda (486) =  ske_difus_Lambda (486) + dN7dz * dN2dz          !
            ske_difus_Lambda (487) =  ske_difus_Lambda (487) + dN7dz * dN3dx          !
            ske_difus_Lambda (488) =  ske_difus_Lambda (488) + dN7dz * dN3dy          !
            ske_difus_Lambda (489) =  ske_difus_Lambda (489) + dN7dz * dN3dz          !
            ske_difus_Lambda (490) =  ske_difus_Lambda (490) + dN7dz * dN4dx          !
            ske_difus_Lambda (491) =  ske_difus_Lambda (491) + dN7dz * dN4dy          !
            ske_difus_Lambda (492) =  ske_difus_Lambda (492) + dN7dz * dN4dz          !
            ske_difus_Lambda (493) =  ske_difus_Lambda (493) + dN7dz * dN5dx          !
            ske_difus_Lambda (494) =  ske_difus_Lambda (494) + dN7dz * dN5dy          !
            ske_difus_Lambda (495) =  ske_difus_Lambda (495) + dN7dz * dN5dz          ! MATRIZ DIFUSIVA
            ske_difus_Lambda (496) =  ske_difus_Lambda (496) + dN7dz * dN6dx          !
            ske_difus_Lambda (497) =  ske_difus_Lambda (497) + dN7dz * dN6dy          !    PARCELA LAMBDA
            ske_difus_Lambda (498) =  ske_difus_Lambda (498) + dN7dz * dN6dz          !
            ske_difus_Lambda (499) =  ske_difus_Lambda (499) + dN7dz * dN7dx          ! REDUCED INTEGRATION
            ske_difus_Lambda (500) =  ske_difus_Lambda (500) + dN7dz * dN7dy          !
            ske_difus_Lambda (501) =  ske_difus_Lambda (501) + dN7dz * dN7dz          !
            ske_difus_Lambda (505) =  ske_difus_Lambda (505) + dN8dx * dN1dx          !
            ske_difus_Lambda (506) =  ske_difus_Lambda (506) + dN8dx * dN1dy          !
            ske_difus_Lambda (507) =  ske_difus_Lambda (507) + dN8dx * dN1dz          !
            ske_difus_Lambda (508) =  ske_difus_Lambda (508) + dN8dx * dN2dx
            ske_difus_Lambda (509) =  ske_difus_Lambda (509) + dN8dx * dN2dy          !
            ske_difus_Lambda (510) =  ske_difus_Lambda (510) + dN8dx * dN2dz          !
            ske_difus_Lambda (511) =  ske_difus_Lambda (511) + dN8dx * dN3dx          !
            ske_difus_Lambda (512) =  ske_difus_Lambda (512) + dN8dx * dN3dy          !
            ske_difus_Lambda (513) =  ske_difus_Lambda (513) + dN8dx * dN3dz          ! MATRIZ DIFUSIVA
            ske_difus_Lambda (514) =  ske_difus_Lambda (514) + dN8dx * dN4dx          !
            ske_difus_Lambda (515) =  ske_difus_Lambda (515) + dN8dx * dN4dy          !    PARCELA LAMBDA
            ske_difus_Lambda (516) =  ske_difus_Lambda (516) + dN8dx * dN4dz          !
            ske_difus_Lambda (517) =  ske_difus_Lambda (517) + dN8dx * dN5dx          ! REDUCED INTEGRATION
            ske_difus_Lambda (518) =  ske_difus_Lambda (518) + dN8dx * dN5dy          !
            ske_difus_Lambda (519) =  ske_difus_Lambda (519) + dN8dx * dN5dz          !
            ske_difus_Lambda (520) =  ske_difus_Lambda (520) + dN8dx * dN6dx          !
            ske_difus_Lambda (521) =  ske_difus_Lambda (521) + dN8dx * dN6dy          !
            ske_difus_Lambda (522) =  ske_difus_Lambda (522) + dN8dx * dN6dz          !
            ske_difus_Lambda (523) =  ske_difus_Lambda (523) + dN8dx * dN7dx          !
            ske_difus_Lambda (524) =  ske_difus_Lambda (524) + dN8dx * dN7dy          !
            ske_difus_Lambda (525) =  ske_difus_Lambda (525) + dN8dx * dN7dz          !
            ske_difus_Lambda (526) =  ske_difus_Lambda (526) + dN8dx * dN8dx          !
            ske_difus_Lambda (529) =  ske_difus_Lambda (529) + dN8dy * dN1dx          !
            ske_difus_Lambda (530) =  ske_difus_Lambda (530) + dN8dy * dN1dy          !
            ske_difus_Lambda (531) =  ske_difus_Lambda (531) + dN8dy * dN1dz          !
            ske_difus_Lambda (532) =  ske_difus_Lambda (532) + dN8dy * dN2dx          !
            ske_difus_Lambda (533) =  ske_difus_Lambda (533) + dN8dy * dN2dy          !
            ske_difus_Lambda (534) =  ske_difus_Lambda (534) + dN8dy * dN2dz          !
            ske_difus_Lambda (535) =  ske_difus_Lambda (535) + dN8dy * dN3dx          !
            ske_difus_Lambda (536) =  ske_difus_Lambda (536) + dN8dy * dN3dy          ! MATRIZ DIFUSIVA
            ske_difus_Lambda (537) =  ske_difus_Lambda (537) + dN8dy * dN3dz          !
            ske_difus_Lambda (538) =  ske_difus_Lambda (538) + dN8dy * dN4dx          !    PARCELA LAMBDA
            ske_difus_Lambda (539) =  ske_difus_Lambda (539) + dN8dy * dN4dy          !
            ske_difus_Lambda (540) =  ske_difus_Lambda (540) + dN8dy * dN4dz          ! REDUCED INTEGRATION
            ske_difus_Lambda (541) =  ske_difus_Lambda (541) + dN8dy * dN5dx          !
            ske_difus_Lambda (542) =  ske_difus_Lambda (542) + dN8dy * dN5dy          !
            ske_difus_Lambda (543) =  ske_difus_Lambda (543) + dN8dy * dN5dz          !
            ske_difus_Lambda (544) =  ske_difus_Lambda (544) + dN8dy * dN6dx          !
            ske_difus_Lambda (545) =  ske_difus_Lambda (545) + dN8dy * dN6dy          !
            ske_difus_Lambda (546) =  ske_difus_Lambda (546) + dN8dy * dN6dz          !
            ske_difus_Lambda (547) =  ske_difus_Lambda (547) + dN8dy * dN7dx          !
            ske_difus_Lambda (548) =  ske_difus_Lambda (548) + dN8dy * dN7dy          !
            ske_difus_Lambda (549) =  ske_difus_Lambda (549) + dN8dy * dN7dz          !
            ske_difus_Lambda (550) =  ske_difus_Lambda (550) + dN8dy * dN8dx          !
            ske_difus_Lambda (551) =  ske_difus_Lambda (551) + dN8dy * dN8dy          !
            ske_difus_Lambda (553) =  ske_difus_Lambda (553) + dN8dz * dN1dx          !
            ske_difus_Lambda (554) =  ske_difus_Lambda (554) + dN8dz * dN1dy          !
            ske_difus_Lambda (555) =  ske_difus_Lambda (555) + dN8dz * dN1dz          !
            ske_difus_Lambda (556) =  ske_difus_Lambda (556) + dN8dz * dN2dx          !
            ske_difus_Lambda (557) =  ske_difus_Lambda (557) + dN8dz * dN2dy          !
            ske_difus_Lambda (558) =  ske_difus_Lambda (558) + dN8dz * dN2dz          ! MATRIZ DIFUSIVA
            ske_difus_Lambda (559) =  ske_difus_Lambda (559) + dN8dz * dN3dx          !
            ske_difus_Lambda (560) =  ske_difus_Lambda (560) + dN8dz * dN3dy          !    PARCELA LAMBDA
            ske_difus_Lambda (561) =  ske_difus_Lambda (561) + dN8dz * dN3dz          !
            ske_difus_Lambda (562) =  ske_difus_Lambda (562) + dN8dz * dN4dx          ! REDUCED INTEGRATION
            ske_difus_Lambda (563) =  ske_difus_Lambda (563) + dN8dz * dN4dy          !
            ske_difus_Lambda (564) =  ske_difus_Lambda (564) + dN8dz * dN4dz          !
            ske_difus_Lambda (565) =  ske_difus_Lambda (565) + dN8dz * dN5dx          !
            ske_difus_Lambda (566) =  ske_difus_Lambda (566) + dN8dz * dN5dy          !
            ske_difus_Lambda (567) =  ske_difus_Lambda (567) + dN8dz * dN5dz          !
            ske_difus_Lambda (568) =  ske_difus_Lambda (568) + dN8dz * dN6dx          !
            ske_difus_Lambda (569) =  ske_difus_Lambda (569) + dN8dz * dN6dy          !
            ske_difus_Lambda (570) =  ske_difus_Lambda (570) + dN8dz * dN6dz          !
            ske_difus_Lambda (571) =  ske_difus_Lambda (571) + dN8dz * dN7dx          !
            ske_difus_Lambda (572) =  ske_difus_Lambda (572) + dN8dz * dN7dy          !
            ske_difus_Lambda (573) =  ske_difus_Lambda (573) + dN8dz * dN7dz          !
            ske_difus_Lambda (574) =  ske_difus_Lambda (574) + dN8dz * dN8dx          !
            ske_difus_Lambda (575) =  ske_difus_Lambda (575) + dN8dz * dN8dy          !
            ske_difus_Lambda (576) =  ske_difus_Lambda (576) + dN8dz * dN8dz          !
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
               stiff_difus_Lambda (ima, iel) = ske_difus_Lambda (ima) * xLambda  * vol              !  Multiplica constantes de cada matriz
               stiff_difus_xMu    (ima, iel) = ske_difus_xMu    (ima) * xMu      * vol *  Wpg       !     inclusive a pondera��o de Gauss
           enddo
       
enddo ! end loop nos elementos

write (iout,200) voltot
write (*,200)    voltot

return

100 format ('***(TRIED2D.F90) volume nao positivo p/ o elemento (',i8,')')
200 format (//,' *** Total volume of the Mesh= ', f10.5,' *** '//)

end subroutine
