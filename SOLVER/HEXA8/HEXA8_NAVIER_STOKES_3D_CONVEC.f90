subroutine hexa8_NavierStokes3D_Convec (stiff,fp,up,u,lm,x,y,z,incid,mtype,prop,numnp,nume,&
nummat,nwk,nd,neq1,neqp,nnoel,ncprop,stiff_Convec,stiff_difus_Lambda, stiff_difus_xMu, voltot)

use modtapes , only: iout,ierror
use modmesh  , only: ID
use modvar   , only: ngl, grav
use modloads , only: fbody
    
IMPLICIT NONE
     
!      Global     
INTEGER, INTENT(in)  :: nwk, nume, nd, neq1, numnp, nummat, nnoel, ncprop, neqp
INTEGER, INTENT(in)  :: mtype(nume), lm(nume,nd), incid(nume,nnoel)
REAL*8 , INTENT(in)  :: x(NUMNP), y(NUMNP), z(NUMNP), prop(nummat,ncprop), up(0:neqp), u(0:neq1-1), &
stiff_Difus_Lambda(nd*nd,nume),stiff_Difus_xMu(nd*nd,nume)
REAL*8 , INTENT(out) :: stiff(nd*nd,nume),stiff_Convec(nd*nd,nume), fp(0:neq1-1), voltot

!      Local
INTEGER              :: i, iel, lmaux(nd), ii, jj, kk, ind, jnd, ieq, ima, inoel, ino, iNgl, col
REAL*8               :: x15, y15, z15,x14, y14, z14, x12, y12, z12, DetJ, DetJ_inv, vol
REAL*8               :: rpg(8), spg(8), tpg(8), Wpg, r, s, t, xRho, xMu, xLambda
REAL*8               :: dN1dx, dN2dx, dN3dx, dN4dx, dN5dx, dN6dx, dN7dx, dN8dx, &
                        dN1dy, dN2dy, dN3dy, dN4dy, dN5dy, dN6dy, dN7dy, dN8dy, &
                        dN1dz, dN2dz, dN3dz, dN4dz, dN5dz, dN6dz, dN7dz, dN8dz
REAL*8               :: N1, N2, N3, N4, N5, N6, N7, N8, u1, u2, u3
REAL*8               :: yz1512, xz1514, xy1514, yz1415, xz1215, xy1215, yz1214, xz1214, xy1214

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

fbody(:)   = 0.d0
fp(:)      = 0.d0
stiff(:,:) = 0.d0
voltot    = 0.d0
    
do iel = 1, nume ! loop nos elementos
    
       ue(:)    = 0.d0      
       fe(:)    = 0.d0      
       lmaux(:) = 0      
       
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

       yz1512 = (y15*z12-y12*z15)      ! Caracteristicas
       xz1514 = (x15*z14-x14*z15)      !    Geometricas
       xy1514 = (-x15*y14+x14*y15)     !   Para Jacobiano
       yz1415 = (y14*z15-z14*y15)      !
       xz1215 = (x12*z15-z12*x15)      !
       xy1215 = (-x12*y15+y12*x15)     !
       yz1214 = (y12*z14-y14*z12)      !
       xz1214 = (-x12*z14+z12*x14)     !
       xy1214 = (x12*y14-y12*x14)      !
       
       DetJ = yz1415*x12+xz1514*y12+xy1514*z12     ! Determinante da transformacao em coordenadas naturais
       
       DetJ_inv = 1.d0/DetJ       ! Inverso do determinante da tranformacao
       vol = DetJ                ! Volume do elemento
       voltot = voltot + vol   ! Contador: volume total da malha
           
        if(vol.le.0.d0) then         !
                write(ierror,100)iel  !  Controle, confere que a volume do elemento e positivo
            stop                      !     (nos do elemento ordenados corretamente)
        endif                         !
           
       ske_Convec      (:) = 0.d0 
       vecForcante     (:) = 0.d0
       do i =1, 8  ! loop pontos gauss - FULL INTEGRATION
             r = rpg(i)
             s = spg(i)
             t = tpg(i)
                           
             N1 = 1 - r - s - t + r*s + s*t + t*r - r*s*t   !
             N2 = r - r*s - t*r + r*s*t                     !
             N3 = r*s - r*s*t                               !
             N4 = s - r*s - s*t + r*s*t                     !   Funcoes de forma para o hexa 8
             N5 = t - s*t - t*r + r*s*t                     !
             N6 = t*r - r*s*t                               !
             N7 = r*s*t                                     !
             N8 = s*t - r*s*t                               !

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
       
             u1 = (u_ant(1,1) * N1) + (u_ant(2,1) * N2) + (u_ant(3,1) * N3) + (u_ant(4,1) * N4) + &
             (u_ant(5,1) * N5) +  (u_ant(6,1) *  N6) + (u_ant(7,1) * N7) +  (u_ant(8,1) *  N8)  ! Calcula velocidades
             u2 = (u_ant(1,2) * N1) + (u_ant(2,2) * N2) + (u_ant(3,2) * N3) + (u_ant(4,2) * N4) + &
             (u_ant(5,2) * N5) +  (u_ant(6,2) *  N6) + (u_ant(7,2) * N7) +  (u_ant(8,2) *  N8)  ! Calcula velocidades
             u3 = (u_ant(1,3) * N1) + (u_ant(2,3) * N2) + (u_ant(3,3) * N3) + (u_ant(4,3) * N4) + &
             (u_ant(5,3) * N5) +  (u_ant(6,3) *  N6) + (u_ant(7,3) * N7) +  (u_ant(8,3) *  N8)

        ske_Convec (1) =  ske_Convec (1) + N1 * u1 * dN1dx + N1 * u2 * dN1dy + N1 * u3 * dN1dz
        ske_Convec (4) =  ske_Convec (4) + N2 * u1 * dN1dx + N2 * u2 * dN1dy + N2 * u3 * dN1dz
        ske_Convec (7) =  ske_Convec (7) + N3 * u1 * dN1dx + N3 * u2 * dN1dy + N3 * u3 * dN1dz           !
        ske_Convec (10) =  ske_Convec (10) + N4 * u1 * dN1dx + N4 * u2 * dN1dy + N4 * u3 * dN1dz         !
        ske_Convec (13) =  ske_Convec (13) + N5 * u1 * dN1dx + N5 * u2 * dN1dy + N5 * u3 * dN1dz         !
        ske_Convec (16) =  ske_Convec (16) + N6 * u1 * dN1dx + N6 * u2 * dN1dy + N6 * u3 * dN1dz         ! MATRIZ CONVECTIVA
        ske_Convec (19) =  ske_Convec (19) + N7 * u1 * dN1dx + N7 * u2 * dN1dy + N7 * u3 * dN1dz         !
        ske_Convec (22) =  ske_Convec (22) + N8 * u1 * dN1dx + N8 * u2 * dN1dy + N8 * u3 * dN1dz         !    PARCELA LAMBDA
        ske_Convec (26) =  ske_Convec (26) + N1 * u1 * dN1dx + N1 * u2 * dN1dy + N1 * u3 * dN1dz         !
        ske_Convec (29) =  ske_Convec (29) + N2 * u1 * dN1dx + N2 * u2 * dN1dy + N2 * u3 * dN1dz         ! FULL INTEGRATION
        ske_Convec (32) =  ske_Convec (32) + N3 * u1 * dN1dx + N3 * u2 * dN1dy + N3 * u3 * dN1dz         !
        ske_Convec (35) =  ske_Convec (35) + N4 * u1 * dN1dx + N4 * u2 * dN1dy + N4 * u3 * dN1dz         !
        ske_Convec (38) =  ske_Convec (38) + N5 * u1 * dN1dx + N5 * u2 * dN1dy + N5 * u3 * dN1dz         !
        ske_Convec (41) =  ske_Convec (41) + N6 * u1 * dN1dx + N6 * u2 * dN1dy + N6 * u3 * dN1dz         !
        ske_Convec (44) =  ske_Convec (44) + N7 * u1 * dN1dx + N7 * u2 * dN1dy + N7 * u3 * dN1dz         !
        ske_Convec (47) =  ske_Convec (47) + N8 * u1 * dN1dx + N8 * u2 * dN1dy + N8 * u3 * dN1dz         !
        ske_Convec (51) =  ske_Convec (51) + N1 * u1 * dN1dx + N1 * u2 * dN1dy + N1 * u3 * dN1dz         !
        ske_Convec (54) =  ske_Convec (54) + N2 * u1 * dN1dx + N2 * u2 * dN1dy + N2 * u3 * dN1dz         !
        ske_Convec (57) =  ske_Convec (57) + N3 * u1 * dN1dx + N3 * u2 * dN1dy + N3 * u3 * dN1dz         !
        ske_Convec (60) =  ske_Convec (60) + N4 * u1 * dN1dx + N4 * u2 * dN1dy + N4 * u3 * dN1dz         !
        ske_Convec (63) =  ske_Convec (63) + N5 * u1 * dN1dx + N5 * u2 * dN1dy + N5 * u3 * dN1dz         !
        ske_Convec (66) =  ske_Convec (66) + N6 * u1 * dN1dx + N6 * u2 * dN1dy + N6 * u3 * dN1dz         !
        ske_Convec (69) =  ske_Convec (69) + N7 * u1 * dN1dx + N7 * u2 * dN1dy + N7 * u3 * dN1dz         !
        ske_Convec (72) =  ske_Convec (72) + N8 * u1 * dN1dx + N8 * u2 * dN1dy + N8 * u3 * dN1dz         !
        ske_Convec (73) =  ske_Convec (73) + N1 * u1 * dN2dx + N1 * u2 * dN2dy + N1 * u3 * dN2dz         !
        ske_Convec (76) =  ske_Convec (76) + N2 * u1 * dN2dx + N2 * u2 * dN2dy + N2 * u3 * dN2dz         !
        ske_Convec (79) =  ske_Convec (79) + N3 * u1 * dN2dx + N3 * u2 * dN2dy + N3 * u3 * dN2dz         !      MATRIZ CONVECTIVA
        ske_Convec (82) =  ske_Convec (82) + N4 * u1 * dN2dx + N4 * u2 * dN2dy + N4 * u3 * dN2dz         !
        ske_Convec (85) =  ske_Convec (85) + N5 * u1 * dN2dx + N5 * u2 * dN2dy + N5 * u3 * dN2dz         !    Non Sym
        ske_Convec (88) =  ske_Convec (88) + N6 * u1 * dN2dx + N6 * u2 * dN2dy + N6 * u3 * dN2dz         !
        ske_Convec (91) =  ske_Convec (91) + N7 * u1 * dN2dx + N7 * u2 * dN2dy + N7 * u3 * dN2dz         ! FULL INTEGRATION
        ske_Convec (94) =  ske_Convec (94) + N8 * u1 * dN2dx + N8 * u2 * dN2dy + N8 * u3 * dN2dz         !
        ske_Convec (98) =  ske_Convec (98) + N1 * u1 * dN2dx + N1 * u2 * dN2dy + N1 * u3 * dN2dz         !
        ske_Convec (101) =  ske_Convec (101) + N2 * u1 * dN2dx + N2 * u2 * dN2dy + N2 * u3 * dN2dz       !
        ske_Convec (104) =  ske_Convec (104) + N3 * u1 * dN2dx + N3 * u2 * dN2dy + N3 * u3 * dN2dz       !
        ske_Convec (107) =  ske_Convec (107) + N4 * u1 * dN2dx + N4 * u2 * dN2dy + N4 * u3 * dN2dz       !
        ske_Convec (110) =  ske_Convec (110) + N5 * u1 * dN2dx + N5 * u2 * dN2dy + N5 * u3 * dN2dz       !
        ske_Convec (113) =  ske_Convec (113) + N6 * u1 * dN2dx + N6 * u2 * dN2dy + N6 * u3 * dN2dz       !
        ske_Convec (116) =  ske_Convec (116) + N7 * u1 * dN2dx + N7 * u2 * dN2dy + N7 * u3 * dN2dz       !
        ske_Convec (119) =  ske_Convec (119) + N8 * u1 * dN2dx + N8 * u2 * dN2dy + N8 * u3 * dN2dz       !
        ske_Convec (123) =  ske_Convec (123) + N1 * u1 * dN2dx + N1 * u2 * dN2dy + N1 * u3 * dN2dz       !
        ske_Convec (126) =  ske_Convec (126) + N2 * u1 * dN2dx + N2 * u2 * dN2dy + N2 * u3 * dN2dz       !
        ske_Convec (129) =  ske_Convec (129) + N3 * u1 * dN2dx + N3 * u2 * dN2dy + N3 * u3 * dN2dz       !
        ske_Convec (132) =  ske_Convec (132) + N4 * u1 * dN2dx + N4 * u2 * dN2dy + N4 * u3 * dN2dz       !
        ske_Convec (135) =  ske_Convec (135) + N5 * u1 * dN2dx + N5 * u2 * dN2dy + N5 * u3 * dN2dz       !
        ske_Convec (138) =  ske_Convec (138) + N6 * u1 * dN2dx + N6 * u2 * dN2dy + N6 * u3 * dN2dz       !
        ske_Convec (141) =  ske_Convec (141) + N7 * u1 * dN2dx + N7 * u2 * dN2dy + N7 * u3 * dN2dz       !
        ske_Convec (144) =  ske_Convec (144) + N8 * u1 * dN2dx + N8 * u2 * dN2dy + N8 * u3 * dN2dz       !      MATRIZ CONVECTIVA
        ske_Convec (145) =  ske_Convec (145) + N1 * u1 * dN3dx + N1 * u2 * dN3dy + N1 * u3 * dN3dz       !
        ske_Convec (148) =  ske_Convec (148) + N2 * u1 * dN3dx + N2 * u2 * dN3dy + N2 * u3 * dN3dz       !    Non Sym
        ske_Convec (151) =  ske_Convec (151) + N3 * u1 * dN3dx + N3 * u2 * dN3dy + N3 * u3 * dN3dz       !
        ske_Convec (154) =  ske_Convec (154) + N4 * u1 * dN3dx + N4 * u2 * dN3dy + N4 * u3 * dN3dz       ! FULL INTEGRATION
        ske_Convec (157) =  ske_Convec (157) + N5 * u1 * dN3dx + N5 * u2 * dN3dy + N5 * u3 * dN3dz       !
        ske_Convec (160) =  ske_Convec (160) + N6 * u1 * dN3dx + N6 * u2 * dN3dy + N6 * u3 * dN3dz       !
        ske_Convec (163) =  ske_Convec (163) + N7 * u1 * dN3dx + N7 * u2 * dN3dy + N7 * u3 * dN3dz       !
        ske_Convec (166) =  ske_Convec (166) + N8 * u1 * dN3dx + N8 * u2 * dN3dy + N8 * u3 * dN3dz       !
        ske_Convec (170) =  ske_Convec (170) + N1 * u1 * dN3dx + N1 * u2 * dN3dy + N1 * u3 * dN3dz       !
        ske_Convec (173) =  ske_Convec (173) + N2 * u1 * dN3dx + N2 * u2 * dN3dy + N2 * u3 * dN3dz       !
        ske_Convec (176) =  ske_Convec (176) + N3 * u1 * dN3dx + N3 * u2 * dN3dy + N3 * u3 * dN3dz       !
        ske_Convec (179) =  ske_Convec (179) + N4 * u1 * dN3dx + N4 * u2 * dN3dy + N4 * u3 * dN3dz       !
        ske_Convec (182) =  ske_Convec (182) + N5 * u1 * dN3dx + N5 * u2 * dN3dy + N5 * u3 * dN3dz       !
        ske_Convec (185) =  ske_Convec (185) + N6 * u1 * dN3dx + N6 * u2 * dN3dy + N6 * u3 * dN3dz       !
        ske_Convec (188) =  ske_Convec (188) + N7 * u1 * dN3dx + N7 * u2 * dN3dy + N7 * u3 * dN3dz       !
        ske_Convec (191) =  ske_Convec (191) + N8 * u1 * dN3dx + N8 * u2 * dN3dy + N8 * u3 * dN3dz       !
        ske_Convec (195) =  ske_Convec (195) + N1 * u1 * dN3dx + N1 * u2 * dN3dy + N1 * u3 * dN3dz       !
        ske_Convec (198) =  ske_Convec (198) + N2 * u1 * dN3dx + N2 * u2 * dN3dy + N2 * u3 * dN3dz       !
        ske_Convec (201) =  ske_Convec (201) + N3 * u1 * dN3dx + N3 * u2 * dN3dy + N3 * u3 * dN3dz       !
        ske_Convec (204) =  ske_Convec (204) + N4 * u1 * dN3dx + N4 * u2 * dN3dy + N4 * u3 * dN3dz       !
        ske_Convec (207) =  ske_Convec (207) + N5 * u1 * dN3dx + N5 * u2 * dN3dy + N5 * u3 * dN3dz       !      MATRIZ CONVECTIVA
        ske_Convec (210) =  ske_Convec (210) + N6 * u1 * dN3dx + N6 * u2 * dN3dy + N6 * u3 * dN3dz       !
        ske_Convec (213) =  ske_Convec (213) + N7 * u1 * dN3dx + N7 * u2 * dN3dy + N7 * u3 * dN3dz       !    Non Sym
        ske_Convec (216) =  ske_Convec (216) + N8 * u1 * dN3dx + N8 * u2 * dN3dy + N8 * u3 * dN3dz       !
        ske_Convec (217) =  ske_Convec (217) + N1 * u1 * dN4dx + N1 * u2 * dN4dy + N1 * u3 * dN4dz       ! FULL INTEGRATION
        ske_Convec (220) =  ske_Convec (220) + N2 * u1 * dN4dx + N2 * u2 * dN4dy + N2 * u3 * dN4dz       !
        ske_Convec (223) =  ske_Convec (223) + N3 * u1 * dN4dx + N3 * u2 * dN4dy + N3 * u3 * dN4dz       !
        ske_Convec (226) =  ske_Convec (226) + N4 * u1 * dN4dx + N4 * u2 * dN4dy + N4 * u3 * dN4dz       !
        ske_Convec (229) =  ske_Convec (229) + N5 * u1 * dN4dx + N5 * u2 * dN4dy + N5 * u3 * dN4dz       !
        ske_Convec (232) =  ske_Convec (232) + N6 * u1 * dN4dx + N6 * u2 * dN4dy + N6 * u3 * dN4dz       !
        ske_Convec (235) =  ske_Convec (235) + N7 * u1 * dN4dx + N7 * u2 * dN4dy + N7 * u3 * dN4dz
        ske_Convec (238) =  ske_Convec (238) + N8 * u1 * dN4dx + N8 * u2 * dN4dy + N8 * u3 * dN4dz       !
        ske_Convec (242) =  ske_Convec (242) + N1 * u1 * dN4dx + N1 * u2 * dN4dy + N1 * u3 * dN4dz       !
        ske_Convec (245) =  ske_Convec (245) + N2 * u1 * dN4dx + N2 * u2 * dN4dy + N2 * u3 * dN4dz       !
        ske_Convec (248) =  ske_Convec (248) + N3 * u1 * dN4dx + N3 * u2 * dN4dy + N3 * u3 * dN4dz       !
        ske_Convec (251) =  ske_Convec (251) + N4 * u1 * dN4dx + N4 * u2 * dN4dy + N4 * u3 * dN4dz       !      MATRIZ CONVECTIVA
        ske_Convec (254) =  ske_Convec (254) + N5 * u1 * dN4dx + N5 * u2 * dN4dy + N5 * u3 * dN4dz       !
        ske_Convec (257) =  ske_Convec (257) + N6 * u1 * dN4dx + N6 * u2 * dN4dy + N6 * u3 * dN4dz       !    Non Sym
        ske_Convec (260) =  ske_Convec (260) + N7 * u1 * dN4dx + N7 * u2 * dN4dy + N7 * u3 * dN4dz       !
        ske_Convec (263) =  ske_Convec (263) + N8 * u1 * dN4dx + N8 * u2 * dN4dy + N8 * u3 * dN4dz       ! FULL INTEGRATION
        ske_Convec (267) =  ske_Convec (267) + N1 * u1 * dN4dx + N1 * u2 * dN4dy + N1 * u3 * dN4dz       !
        ske_Convec (270) =  ske_Convec (270) + N2 * u1 * dN4dx + N2 * u2 * dN4dy + N2 * u3 * dN4dz       !
        ske_Convec (273) =  ske_Convec (273) + N3 * u1 * dN4dx + N3 * u2 * dN4dy + N3 * u3 * dN4dz       !
        ske_Convec (276) =  ske_Convec (276) + N4 * u1 * dN4dx + N4 * u2 * dN4dy + N4 * u3 * dN4dz       !
        ske_Convec (279) =  ske_Convec (279) + N5 * u1 * dN4dx + N5 * u2 * dN4dy + N5 * u3 * dN4dz       !
        ske_Convec (282) =  ske_Convec (282) + N6 * u1 * dN4dx + N6 * u2 * dN4dy + N6 * u3 * dN4dz       !
        ske_Convec (285) =  ske_Convec (285) + N7 * u1 * dN4dx + N7 * u2 * dN4dy + N7 * u3 * dN4dz       !
        ske_Convec (288) =  ske_Convec (288) + N8 * u1 * dN4dx + N8 * u2 * dN4dy + N8 * u3 * dN4dz       !
        ske_Convec (289) =  ske_Convec (289) + N1 * u1 * dN5dx + N1 * u2 * dN5dy + N1 * u3 * dN5dz       !
        ske_Convec (292) =  ske_Convec (292) + N2 * u1 * dN5dx + N2 * u2 * dN5dy + N2 * u3 * dN5dz       !
        ske_Convec (295) =  ske_Convec (295) + N3 * u1 * dN5dx + N3 * u2 * dN5dy + N3 * u3 * dN5dz       !
        ske_Convec (298) =  ske_Convec (298) + N4 * u1 * dN5dx + N4 * u2 * dN5dy + N4 * u3 * dN5dz       !
        ske_Convec (301) =  ske_Convec (301) + N5 * u1 * dN5dx + N5 * u2 * dN5dy + N5 * u3 * dN5dz       !
        ske_Convec (304) =  ske_Convec (304) + N6 * u1 * dN5dx + N6 * u2 * dN5dy + N6 * u3 * dN5dz       !
        ske_Convec (307) =  ske_Convec (307) + N7 * u1 * dN5dx + N7 * u2 * dN5dy + N7 * u3 * dN5dz       !
        ske_Convec (310) =  ske_Convec (310) + N8 * u1 * dN5dx + N8 * u2 * dN5dy + N8 * u3 * dN5dz       !
        ske_Convec (314) =  ske_Convec (314) + N1 * u1 * dN5dx + N1 * u2 * dN5dy + N1 * u3 * dN5dz       !      MATRIZ CONVECTIVA
        ske_Convec (317) =  ske_Convec (317) + N2 * u1 * dN5dx + N2 * u2 * dN5dy + N2 * u3 * dN5dz       !
        ske_Convec (320) =  ske_Convec (320) + N3 * u1 * dN5dx + N3 * u2 * dN5dy + N3 * u3 * dN5dz       !    Non Sym
        ske_Convec (323) =  ske_Convec (323) + N4 * u1 * dN5dx + N4 * u2 * dN5dy + N4 * u3 * dN5dz       !
        ske_Convec (326) =  ske_Convec (326) + N5 * u1 * dN5dx + N5 * u2 * dN5dy + N5 * u3 * dN5dz       ! FULL INTEGRATION
        ske_Convec (329) =  ske_Convec (329) + N6 * u1 * dN5dx + N6 * u2 * dN5dy + N6 * u3 * dN5dz       !
        ske_Convec (332) =  ske_Convec (332) + N7 * u1 * dN5dx + N7 * u2 * dN5dy + N7 * u3 * dN5dz       !
        ske_Convec (335) =  ske_Convec (335) + N8 * u1 * dN5dx + N8 * u2 * dN5dy + N8 * u3 * dN5dz       !
        ske_Convec (339) =  ske_Convec (339) + N1 * u1 * dN5dx + N1 * u2 * dN5dy + N1 * u3 * dN5dz       !
        ske_Convec (342) =  ske_Convec (342) + N2 * u1 * dN5dx + N2 * u2 * dN5dy + N2 * u3 * dN5dz       !
        ske_Convec (345) =  ske_Convec (345) + N3 * u1 * dN5dx + N3 * u2 * dN5dy + N3 * u3 * dN5dz       !
        ske_Convec (348) =  ske_Convec (348) + N4 * u1 * dN5dx + N4 * u2 * dN5dy + N4 * u3 * dN5dz       !
        ske_Convec (351) =  ske_Convec (351) + N5 * u1 * dN5dx + N5 * u2 * dN5dy + N5 * u3 * dN5dz       !
        ske_Convec (354) =  ske_Convec (354) + N6 * u1 * dN5dx + N6 * u2 * dN5dy + N6 * u3 * dN5dz       !
        ske_Convec (357) =  ske_Convec (357) + N7 * u1 * dN5dx + N7 * u2 * dN5dy + N7 * u3 * dN5dz       !
        ske_Convec (360) =  ske_Convec (360) + N8 * u1 * dN5dx + N8 * u2 * dN5dy + N8 * u3 * dN5dz       !
        ske_Convec (361) =  ske_Convec (361) + N1 * u1 * dN6dx + N1 * u2 * dN6dy + N1 * u3 * dN6dz       !
        ske_Convec (364) =  ske_Convec (364) + N2 * u1 * dN6dx + N2 * u2 * dN6dy + N2 * u3 * dN6dz       !
        ske_Convec (367) =  ske_Convec (367) + N3 * u1 * dN6dx + N3 * u2 * dN6dy + N3 * u3 * dN6dz       !
        ske_Convec (370) =  ske_Convec (370) + N4 * u1 * dN6dx + N4 * u2 * dN6dy + N4 * u3 * dN6dz       !
        ske_Convec (373) =  ske_Convec (373) + N5 * u1 * dN6dx + N5 * u2 * dN6dy + N5 * u3 * dN6dz       !
        ske_Convec (376) =  ske_Convec (376) + N6 * u1 * dN6dx + N6 * u2 * dN6dy + N6 * u3 * dN6dz       !      MATRIZ CONVECTIVA
        ske_Convec (379) =  ske_Convec (379) + N7 * u1 * dN6dx + N7 * u2 * dN6dy + N7 * u3 * dN6dz       !
        ske_Convec (382) =  ske_Convec (382) + N8 * u1 * dN6dx + N8 * u2 * dN6dy + N8 * u3 * dN6dz       !    Non Sym
        ske_Convec (386) =  ske_Convec (386) + N1 * u1 * dN6dx + N1 * u2 * dN6dy + N1 * u3 * dN6dz       !
        ske_Convec (389) =  ske_Convec (389) + N2 * u1 * dN6dx + N2 * u2 * dN6dy + N2 * u3 * dN6dz       ! FULL INTEGRATION
        ske_Convec (392) =  ske_Convec (392) + N3 * u1 * dN6dx + N3 * u2 * dN6dy + N3 * u3 * dN6dz       !
        ske_Convec (395) =  ske_Convec (395) + N4 * u1 * dN6dx + N4 * u2 * dN6dy + N4 * u3 * dN6dz       !
        ske_Convec (398) =  ske_Convec (398) + N5 * u1 * dN6dx + N5 * u2 * dN6dy + N5 * u3 * dN6dz       !
        ske_Convec (401) =  ske_Convec (401) + N6 * u1 * dN6dx + N6 * u2 * dN6dy + N6 * u3 * dN6dz       !
        ske_Convec (404) =  ske_Convec (404) + N7 * u1 * dN6dx + N7 * u2 * dN6dy + N7 * u3 * dN6dz       !
        ske_Convec (407) =  ske_Convec (407) + N8 * u1 * dN6dx + N8 * u2 * dN6dy + N8 * u3 * dN6dz       !
        ske_Convec (411) =  ske_Convec (411) + N1 * u1 * dN6dx + N1 * u2 * dN6dy + N1 * u3 * dN6dz       !
        ske_Convec (414) =  ske_Convec (414) + N2 * u1 * dN6dx + N2 * u2 * dN6dy + N2 * u3 * dN6dz       !
        ske_Convec (417) =  ske_Convec (417) + N3 * u1 * dN6dx + N3 * u2 * dN6dy + N3 * u3 * dN6dz       !
        ske_Convec (420) =  ske_Convec (420) + N4 * u1 * dN6dx + N4 * u2 * dN6dy + N4 * u3 * dN6dz       !
        ske_Convec (423) =  ske_Convec (423) + N5 * u1 * dN6dx + N5 * u2 * dN6dy + N5 * u3 * dN6dz       !
        ske_Convec (426) =  ske_Convec (426) + N6 * u1 * dN6dx + N6 * u2 * dN6dy + N6 * u3 * dN6dz       !
        ske_Convec (429) =  ske_Convec (429) + N7 * u1 * dN6dx + N7 * u2 * dN6dy + N7 * u3 * dN6dz       !
        ske_Convec (432) =  ske_Convec (432) + N8 * u1 * dN6dx + N8 * u2 * dN6dy + N8 * u3 * dN6dz       !
        ske_Convec (433) =  ske_Convec (433) + N1 * u1 * dN7dx + N1 * u2 * dN7dy + N1 * u3 * dN7dz       !
        ske_Convec (436) =  ske_Convec (436) + N2 * u1 * dN7dx + N2 * u2 * dN7dy + N2 * u3 * dN7dz       !
        ske_Convec (439) =  ske_Convec (439) + N3 * u1 * dN7dx + N3 * u2 * dN7dy + N3 * u3 * dN7dz       !      MATRIZ CONVECTIVA
        ske_Convec (442) =  ske_Convec (442) + N4 * u1 * dN7dx + N4 * u2 * dN7dy + N4 * u3 * dN7dz       !
        ske_Convec (445) =  ske_Convec (445) + N5 * u1 * dN7dx + N5 * u2 * dN7dy + N5 * u3 * dN7dz       !    Non Sym
        ske_Convec (448) =  ske_Convec (448) + N6 * u1 * dN7dx + N6 * u2 * dN7dy + N6 * u3 * dN7dz       !
        ske_Convec (451) =  ske_Convec (451) + N7 * u1 * dN7dx + N7 * u2 * dN7dy + N7 * u3 * dN7dz       ! FULL INTEGRATION
        ske_Convec (454) =  ske_Convec (454) + N8 * u1 * dN7dx + N8 * u2 * dN7dy + N8 * u3 * dN7dz       !
        ske_Convec (458) =  ske_Convec (458) + N1 * u1 * dN7dx + N1 * u2 * dN7dy + N1 * u3 * dN7dz       !
        ske_Convec (461) =  ske_Convec (461) + N2 * u1 * dN7dx + N2 * u2 * dN7dy + N2 * u3 * dN7dz       !
        ske_Convec (464) =  ske_Convec (464) + N3 * u1 * dN7dx + N3 * u2 * dN7dy + N3 * u3 * dN7dz       !
        ske_Convec (467) =  ske_Convec (467) + N4 * u1 * dN7dx + N4 * u2 * dN7dy + N4 * u3 * dN7dz       !
        ske_Convec (470) =  ske_Convec (470) + N5 * u1 * dN7dx + N5 * u2 * dN7dy + N5 * u3 * dN7dz
        ske_Convec (473) =  ske_Convec (473) + N6 * u1 * dN7dx + N6 * u2 * dN7dy + N6 * u3 * dN7dz       !
        ske_Convec (476) =  ske_Convec (476) + N7 * u1 * dN7dx + N7 * u2 * dN7dy + N7 * u3 * dN7dz       !
        ske_Convec (479) =  ske_Convec (479) + N8 * u1 * dN7dx + N8 * u2 * dN7dy + N8 * u3 * dN7dz       !
        ske_Convec (483) =  ske_Convec (483) + N1 * u1 * dN7dx + N1 * u2 * dN7dy + N1 * u3 * dN7dz       !
        ske_Convec (486) =  ske_Convec (486) + N2 * u1 * dN7dx + N2 * u2 * dN7dy + N2 * u3 * dN7dz       !      MATRIZ CONVECTIVA
        ske_Convec (489) =  ske_Convec (489) + N3 * u1 * dN7dx + N3 * u2 * dN7dy + N3 * u3 * dN7dz       !
        ske_Convec (492) =  ske_Convec (492) + N4 * u1 * dN7dx + N4 * u2 * dN7dy + N4 * u3 * dN7dz       !    Non Sym
        ske_Convec (495) =  ske_Convec (495) + N5 * u1 * dN7dx + N5 * u2 * dN7dy + N5 * u3 * dN7dz       !
        ske_Convec (498) =  ske_Convec (498) + N6 * u1 * dN7dx + N6 * u2 * dN7dy + N6 * u3 * dN7dz       ! FULL INTEGRATION
        ske_Convec (501) =  ske_Convec (501) + N7 * u1 * dN7dx + N7 * u2 * dN7dy + N7 * u3 * dN7dz       !
        ske_Convec (504) =  ske_Convec (504) + N8 * u1 * dN7dx + N8 * u2 * dN7dy + N8 * u3 * dN7dz       !
        ske_Convec (505) =  ske_Convec (505) + N1 * u1 * dN8dx + N1 * u2 * dN8dy + N1 * u3 * dN8dz       !
        ske_Convec (508) =  ske_Convec (508) + N2 * u1 * dN8dx + N2 * u2 * dN8dy + N2 * u3 * dN8dz       !
        ske_Convec (511) =  ske_Convec (511) + N3 * u1 * dN8dx + N3 * u2 * dN8dy + N3 * u3 * dN8dz       !
        ske_Convec (514) =  ske_Convec (514) + N4 * u1 * dN8dx + N4 * u2 * dN8dy + N4 * u3 * dN8dz       !
        ske_Convec (517) =  ske_Convec (517) + N5 * u1 * dN8dx + N5 * u2 * dN8dy + N5 * u3 * dN8dz       !
        ske_Convec (520) =  ske_Convec (520) + N6 * u1 * dN8dx + N6 * u2 * dN8dy + N6 * u3 * dN8dz       !
        ske_Convec (523) =  ske_Convec (523) + N7 * u1 * dN8dx + N7 * u2 * dN8dy + N7 * u3 * dN8dz       !
        ske_Convec (526) =  ske_Convec (526) + N8 * u1 * dN8dx + N8 * u2 * dN8dy + N8 * u3 * dN8dz       !
        ske_Convec (530) =  ske_Convec (530) + N1 * u1 * dN8dx + N1 * u2 * dN8dy + N1 * u3 * dN8dz       !
        ske_Convec (533) =  ske_Convec (533) + N2 * u1 * dN8dx + N2 * u2 * dN8dy + N2 * u3 * dN8dz       !
        ske_Convec (536) =  ske_Convec (536) + N3 * u1 * dN8dx + N3 * u2 * dN8dy + N3 * u3 * dN8dz       !
        ske_Convec (539) =  ske_Convec (539) + N4 * u1 * dN8dx + N4 * u2 * dN8dy + N4 * u3 * dN8dz       !
        ske_Convec (542) =  ske_Convec (542) + N5 * u1 * dN8dx + N5 * u2 * dN8dy + N5 * u3 * dN8dz       !
        ske_Convec (545) =  ske_Convec (545) + N6 * u1 * dN8dx + N6 * u2 * dN8dy + N6 * u3 * dN8dz       !
        ske_Convec (548) =  ske_Convec (548) + N7 * u1 * dN8dx + N7 * u2 * dN8dy + N7 * u3 * dN8dz       !      MATRIZ CONVECTIVA
        ske_Convec (551) =  ske_Convec (551) + N8 * u1 * dN8dx + N8 * u2 * dN8dy + N8 * u3 * dN8dz       !
        ske_Convec (555) =  ske_Convec (555) + N1 * u1 * dN8dx + N1 * u2 * dN8dy + N1 * u3 * dN8dz       !    Non Sym
        ske_Convec (558) =  ske_Convec (558) + N2 * u1 * dN8dx + N2 * u2 * dN8dy + N2 * u3 * dN8dz       !
        ske_Convec (561) =  ske_Convec (561) + N3 * u1 * dN8dx + N3 * u2 * dN8dy + N3 * u3 * dN8dz       ! FULL INTEGRATION
        ske_Convec (564) =  ske_Convec (564) + N4 * u1 * dN8dx + N4 * u2 * dN8dy + N4 * u3 * dN8dz       !
        ske_Convec (567) =  ske_Convec (567) + N5 * u1 * dN8dx + N5 * u2 * dN8dy + N5 * u3 * dN8dz       !
        ske_Convec (570) =  ske_Convec (570) + N6 * u1 * dN8dx + N6 * u2 * dN8dy + N6 * u3 * dN8dz       !
        ske_Convec (573) =  ske_Convec (573) + N7 * u1 * dN8dx + N7 * u2 * dN8dy + N7 * u3 * dN8dz       !
        ske_Convec (576) =  ske_Convec (576) + N8 * u1 * dN8dx + N8 * u2 * dN8dy + N8 * u3 * dN8dz       !
       
             !!!!!!!!!!!!!!!! Vetor forcante de corpo (GRAVIDADE)
             
             !vecForcante(1)  = vecForcante(1)   + N1 * grav !(1)
             !vecForcante(2)  = vecForcante(2)   + N1 * grav !(2)
              vecForcante(3)  = vecForcante(3)   + N1 * grav !(3)
             !vecForcante(4)  = vecForcante(4)   + N2 * grav !(1)
             !vecForcante(5)  = vecForcante(5)   + N2 * grav !(2)
              vecForcante(6)  = vecForcante(6)   + N2 * grav !(3)
             !vecForcante(7)  = vecForcante(7)   + N3 * grav !(1)
             !vecForcante(8)  = vecForcante(8)   + N3 * grav !(2)
              vecForcante(9)  = vecForcante(9)   + N3 * grav !(3)
             !vecForcante(10) = vecForcante(10)  + N4 * grav !(1)
             !vecForcante(11) = vecForcante(11)  + N4 * grav !(2)
              vecForcante(12) = vecForcante(12)  + N4 * grav !(3)
             !vecForcante(13) = vecForcante(13)  + N5 * grav !(1)
             !vecForcante(14) = vecForcante(14)  + N5 * grav !(2)
              vecForcante(15) = vecForcante(15)  + N5 * grav !(3)
             !vecForcante(16) = vecForcante(16)  + N6 * grav !(1)
             !vecForcante(17) = vecForcante(17)  + N6 * grav !(2)
              vecForcante(18) = vecForcante(18)  + N6 * grav !(3)
             !vecForcante(19) = vecForcante(16)  + N7 * grav !(1)
             !vecForcante(20) = vecForcante(17)  + N7 * grav !(2)
              vecForcante(21) = vecForcante(18)  + N7 * grav !(3)
             !vecForcante(22) = vecForcante(16)  + N8 * grav !(1)
             !vecForcante(23) = vecForcante(17)  + N8 * grav !(2)
              vecForcante(24) = vecForcante(18)  + N8 * grav !(3)
             
       enddo !Fim do loop nos pontos de Gauss FULL INTEGRATION
           
       !!!!! Multiplica constantes
       do ima = 1, nd * nd
           stiff_convec       (ima, iel) = ske_convec       (ima) * xRho     * vol *  Wpg       !       no caso de Full Integration
       enddo
       do ind = 1, nd
           vecForcante (ind) = vecForcante(ind) * xRho * vol * Wpg
       enddo
       
       !!!!! Monta Matriz difusiva para multiplicar o termo for�ante
       do ima = 1, nd*nd                                                     
           stiff       (ima, iel) =  stiff_difus_Lambda (ima, iel) + stiff_difus_xMu (ima, iel) + stiff_convec (ima, iel)   !!!CONVECCAO NAO ESTA SOMANDO NO TERMO FOR�ANTE
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
           stiff (ima, iel) =  stiff (ima, iel)  !+ stiff_convec (ima, iel)      !   soma das parcelas convectiva e difusivas  !!! CONVECCAO NAO ESTA SOMANDO NO TERMO FORCANTE
       enddo    
       
enddo ! end loop nos elementos

!write (iout,200) voltot
!write (*,200)    voltot

deallocate (sk )
deallocate (fe)
deallocate (ue)
deallocate (ske_Convec)



deallocate (vecforcante)
deallocate (u_ant)

100 format ('***(TRIED2D.F90) volume nao positivo p/ o elemento (',i8,')')
!200 format (//,' *** Total Volume of the Mesh= ', f10.5,' *** '//)

return
end subroutine
