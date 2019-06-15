subroutine Calc_Pressure_Hexa8(Pressure,TFUNC,N,Time)

use    modsolver, only: u,Forca, Fint, sig11, sig12, sig13, sig22, sig23, sig33
use    modloads , only: up
use    modvar   , only: nume, nnoel, numnp, nummat,nd,ncprop,neq,neqp,ngl,RefPressure,elrefp
use    modmesh  , only: incid, id, x, y, z, prop, mtype, lm

      implicit none
    
!   Global
REAL*8               :: Pressure(nume)
      
!   Local
INTEGER              :: i, iel, inoel, ino, ieq(ngl), jeq,N, ieldivUmax, ieldivUmin, ind
REAL*8               :: x15, y15, z15,x14, y14, z14, x12, y12, z12, DetJ, DetJ_inv, xLambda, xMu, TFUNC, vol,voltot, Time
REAL*8               :: r, s, t, dNx(nnoel), dNy(nnoel), Vel_inoel(ngl), Wpg
REAL*8               :: dUxdx, dUxdy, dUxdz, dUydx, dUydy, dUydz, dUzdx, dUzdy, dUzdz
REAL*8               :: dNdx(nnoel),dNdy(nnoel),dNdz(nnoel)
REAL*8               :: dUx, dUy, divUmax, divUmin
REAL*8               :: yz1512, xz1514, xy1514, yz1415, xz1215, xy1215, yz1214, xz1214, xy1214
character*70 :: Filename1

FILENAME1 = './OUT1/div.OUT'
if (N .eq. 1) then
    OPEN  (UNIT=22,FILE=FILENAME1, form = 'formatted')  
    write (22,'(a, 1e12.4E3)') 'Tempo, Pressao minima, Elemento de div(u) maximo, &
    Pressao maxima, Elemento de div(u) minimo, Lambda = ', prop(mtype(1),6)
else
    OPEN (UNIT=22,FILE=FILENAME1, position='append', Status= 'unknown') 
endif
    
!   Pontos de Gauss:

          !Wpg = 0.12500d0

          !rpg(1) =   0.2113248654
          !rpg(2) =   0.7886751346
          !rpg(3) =   0.7886751346
          !rpg(4) =   0.2113248654
          !rpg(5) =   0.2113248654
          !rpg(6) =   0.7886751346
          !rpg(7) =   0.7886751346
          !rpg(8) =   0.2113248654

          !spg(1) =   0.2113248654
          !spg(2) =   0.2113248654
          !spg(3) =   0.7886751346
          !spg(4) =   0.7886751346
          !spg(5) =   0.2113248654
          !spg(6) =   0.2113248654
          !spg(7) =   0.7886751346
          !spg(8) =   0.7886751346

          !tpg(1) =   0.2113248654
          !tpg(2) =   0.2113248654
          !tpg(3) =   0.2113248654
          !tpg(4) =   0.2113248654
          !tpg(5) =   0.7886751346
          !tpg(6) =   0.7886751346
          !tpg(7) =   0.7886751346
          !tpg(8) =   0.7886751346

divUmax = 0.d0
divUmin = 0.d0
forca(:,:) = 0.d0
     
do iel = 1, nume ! loop nos elementos
         
     xLambda = prop(mtype(iel),6)       ! Coeficiente de incompressibilidade
     xMu     = prop(mtype(iel),5)       ! Viscosidade din�mica
     
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
       
     dUxdx = 0.d0
     dUxdy = 0.d0
     dUxdz = 0.d0
     dUydx = 0.d0
     dUydy = 0.d0
     dUydz = 0.d0
     dUzdx = 0.d0
     dUzdy = 0.d0
     dUzdz = 0.d0
     
     ! do i =1, 3  ! loop pontos gauss

             r = 0.5d0
             s = 0.5d0
             t = 0.5d0

            dNdx(1) = (yz1415*(-1+s+t-s*t)+yz1512*(-1+r+t-t*r)+yz1214*(-1+s+r-r*s))*DetJ_inv  !
            dNdy(1) = (xz1514*(-1+s+t-s*t)+xz1215*(-1+r+t-t*r)+xz1214*(-1+s+r-r*s))*DetJ_inv  !
            dNdz(1) = (xy1514*(-1+s+t-s*t)+xy1215*(-1+r+t-t*r)+xy1214*(-1+s+r-r*s))*DetJ_inv  !
                                                                                                                       !
            dNdx(2) = (yz1415*(1-s-t+s*t)+yz1512*(-r+t*r)+yz1214*(-r+r*s)         )*DetJ_inv  !
            dNdy(2) = (xz1514*(1-s-t+s*t)+xz1215*(-r+t*r)+xz1214*(-r+r*s)         )*DetJ_inv  !    Primeiras derivadas da funcao de forma
            dNdz(2) = (xy1514*(1-s-t+s*t)+xy1215*(-r+t*r)+xy1214*(-r+r*s)         )*DetJ_inv  !         avaliadas nos pontos de Gauss
                                                                                                                       !
            dNdx(3) = (yz1415*(s-s*t)+yz1512*(r-t*r)-yz1214*r*s                   )*DetJ_inv  !
            dNdy(3) = (xz1514*(s-s*t)+xz1215*(r-t*r)-xz1214*r*s                   )*DetJ_inv  !
            dNdz(3) = (xy1514*(s-s*t)+xy1215*(r-t*r)-xy1214*r*s                   )*DetJ_inv  !
                                                                                                                      !
            dNdx(4) = (yz1415*(-s+s*t)+yz1512*(1-r-t+t*r)+yz1214*(-s+r*s)         )*DetJ_inv  !
            dNdy(4) = (xz1514*(-s+s*t)+xz1215*(1-r-t+t*r)+xz1214*(-s+r*s)         )*DetJ_inv  !
            dNdz(4) = (xy1514*(-s+s*t)+xy1215*(1-r-t+t*r)+xy1214*(-s+r*s)         )*DetJ_inv  !
                                                                                                                       !
            dNdx(5) = (yz1415*(-t+s*t)+yz1512*(-t+t*r)+yz1214*(1-s-r+r*s)         )*DetJ_inv  !
            dNdy(5) = (xz1514*(-t+s*t)+xz1215*(-t+t*r)+xz1214*(1-s-r+r*s)         )*DetJ_inv  !    Primeiras derivadas da funcao de forma
            dNdz(5) = (xy1514*(-t+s*t)+xy1215*(-t+t*r)+xy1214*(1-s-r+r*s)         )*DetJ_inv  !         avaliadas nos pontos de Gauss
                                                                                                                       !
            dNdx(6) = (yz1415*(t-s*t)-yz1512*t*r+yz1214*(r-r*s)                   )*DetJ_inv  !
            dNdy(6) = (xz1514*(t-s*t)-xz1215*t*r+xz1214*(r-r*s)                   )*DetJ_inv  !
            dNdz(6) = (xy1514*(t-s*t)-xy1215*t*r+xy1214*(r-r*s)                   )*DetJ_inv  !
                                                                                                                       !
            dNdx(7) = (yz1415*s*t+yz1512*t*r+yz1214*r*s                           )*DetJ_inv  !
            dNdy(7) = (xz1514*s*t+xz1215*t*r+xz1214*r*s                           )*DetJ_inv  !
            dNdz(7) = (xy1514*s*t+xy1215*t*r+xy1214*r*s                           )*DetJ_inv  !
                                                                                                                       !
            dNdx(8) = -(yz1415*s*t+yz1512*(t-t*r)+yz1214*(s-r*s)                  )*DetJ_inv  !
            dNdy(8) = -(xz1514*s*t+xz1215*(t-t*r)+xz1214*(s-r*s)                  )*DetJ_inv  !    Primeiras derivadas da funcao de forma
            dNdz(8) = -(xy1514*s*t+xy1215*(t-t*r)+xy1214*(s-r*s)                  )*DetJ_inv  !         avaliadas nos pontos de Gauss
                                                                                                                     !                                                                                                     !
             do inoel = 1, nnoel  ! Loop nos n�s do elemento
                  ino = incid (iel,inoel)   ! Indice do n� a qual se refere este Loop
                  ieq(1) = id(1,ino) ! Indice da equa��o referente ao primeiro grau de liberdade deste n�
                  ieq(2) = id(2,ino) ! Indice da equa��o referente ao segundo grau de liberdade deste n�
                  ieq(3) = id(3,ino) ! Indice da equa��o referente ao terceiro grau de liberdade deste n�
                  if (ieq(1).ge.0) then     ! Grau de liberdade n�o prescrito
	                  vel_inoel(1) = u(ieq(1))
                  else                      ! Grau de liberdade prescrito
                      vel_inoel(1) = up(-ieq(1)) * TFUNC
                  endif
                  if (ieq(2).ge.0) then
	                  vel_inoel(2) = u(ieq(2))
                  else
                      vel_inoel(2) = up(-ieq(2)) * TFUNC
                  endif
                  if (ieq(3).ge.0) then
                      vel_inoel(3) = u(ieq(3))
                  else
                      vel_inoel(3) = up(-ieq(3)) * TFUNC
                  endif
                 dUxdx =  dUxdx + dNdx(inoel) * vel_inoel(1)
                 dUxdy =  dUxdy + dNdy(inoel) * vel_inoel(1)
                 dUxdz =  dUxdz + dNdz(inoel) * vel_inoel(1)
                 dUydx =  dUydx + dNdx(inoel) * vel_inoel(2)
                 dUydy =  dUydy + dNdy(inoel) * vel_inoel(2)
                 dUydz =  dUydz + dNdz(inoel) * vel_inoel(2)
                 dUzdx =  dUzdx + dNdx(inoel) * vel_inoel(3)
                 dUzdy =  dUzdy + dNdy(inoel) * vel_inoel(3)
                 dUzdz =  dUzdz + dNdz(inoel) * vel_inoel(3)
             enddo
     ! enddo
     
             Pressure(iel) = - xLambda * (dUxdx + dUydy + dUzdz) ! * Wpg
             sig11(iel) =  2 * xMu * (dUxdx)
             sig12(iel) =      xMu * (dUxdy + dUydx)
             sig13(iel) =      xMu * (dUxdz + dUzdx)
             sig22(iel) =  2 * xMu * (dUydy)
             sig23(iel) =      xMu * (dUydz + dUzdy)
             sig33(iel) =  2 * xMu * (dUzdz)


            Fint(1 ,iel) = ( dNdx(1) * sig11(iel) + dNdy(1) * sig12(iel) + dNdz(1) * sig13(iel) ) * vol
            Fint(2 ,iel) = ( dNdx(1) * sig12(iel) + dNdy(1) * sig22(iel) + dNdz(1) * sig23(iel) ) * vol
            Fint(3 ,iel) = ( dNdx(1) * sig13(iel) + dNdy(1) * sig23(iel) + dNdz(1) * sig33(iel) ) * vol

            Fint(4 ,iel) = ( dNdx(2) * sig11(iel) + dNdy(2) * sig12(iel) + dNdz(2) * sig13(iel) ) * vol
            Fint(5 ,iel) = ( dNdx(2) * sig12(iel) + dNdy(2) * sig22(iel) + dNdz(2) * sig23(iel) ) * vol
            Fint(6 ,iel) = ( dNdx(2) * sig13(iel) + dNdy(2) * sig23(iel) + dNdz(2) * sig33(iel) ) * vol

            Fint(7 ,iel) = ( dNdx(3) * sig11(iel) + dNdy(3) * sig12(iel) + dNdz(3) * sig13(iel) ) * vol
            Fint(8 ,iel) = ( dNdx(3) * sig12(iel) + dNdy(3) * sig22(iel) + dNdz(3) * sig23(iel) ) * vol
            Fint(9 ,iel) = ( dNdx(3) * sig13(iel) + dNdy(3) * sig23(iel) + dNdz(3) * sig33(iel) ) * vol

            Fint(10,iel) = ( dNdx(4) * sig11(iel) + dNdy(4) * sig12(iel) + dNdz(4) * sig13(iel) ) * vol
            Fint(11,iel) = ( dNdx(4) * sig12(iel) + dNdy(4) * sig22(iel) + dNdz(4) * sig23(iel) ) * vol
            Fint(12,iel) = ( dNdx(4) * sig13(iel) + dNdy(4) * sig23(iel) + dNdz(4) * sig33(iel) ) * vol

            Fint(13,iel) = ( dNdx(5) * sig11(iel) + dNdy(5) * sig12(iel) + dNdz(5) * sig13(iel) ) * vol
            Fint(14,iel) = ( dNdx(5) * sig12(iel) + dNdy(5) * sig22(iel) + dNdz(5) * sig23(iel) ) * vol
            Fint(15,iel) = ( dNdx(5) * sig13(iel) + dNdy(5) * sig23(iel) + dNdz(5) * sig33(iel) ) * vol

            Fint(16,iel) = ( dNdx(6) * sig11(iel) + dNdy(6) * sig12(iel) + dNdz(6) * sig13(iel) ) * vol
            Fint(17,iel) = ( dNdx(6) * sig12(iel) + dNdy(6) * sig22(iel) + dNdz(6) * sig23(iel) ) * vol
            Fint(18,iel) = ( dNdx(6) * sig13(iel) + dNdy(6) * sig23(iel) + dNdz(6) * sig33(iel) ) * vol

            Fint(19,iel) = ( dNdx(7) * sig11(iel) + dNdy(7) * sig12(iel) + dNdz(7) * sig13(iel) ) * vol
            Fint(20,iel) = ( dNdx(7) * sig12(iel) + dNdy(7) * sig22(iel) + dNdz(7) * sig23(iel) ) * vol
            Fint(21,iel) = ( dNdx(7) * sig13(iel) + dNdy(7) * sig23(iel) + dNdz(7) * sig33(iel) ) * vol

            Fint(22,iel) = ( dNdx(8) * sig11(iel) + dNdy(8) * sig12(iel) + dNdz(8) * sig13(iel) ) * vol
            Fint(23,iel) = ( dNdx(8) * sig12(iel) + dNdy(8) * sig22(iel) + dNdz(8) * sig23(iel) ) * vol
            Fint(24,iel) = ( dNdx(8) * sig13(iel) + dNdy(8) * sig23(iel) + dNdz(8) * sig33(iel) ) * vol

           do inoel =1, nnoel
               ino = incid(iel,inoel)
               ind = ((inoel-1) * ngl + 1)
               Forca(ino,1) = Forca(ino,1) + Fint(ind,iel)
               ind = ((inoel-1) * ngl + 2)
               Forca(ino,2) = Forca(ino,2) + Fint(ind,iel)
               ind = ((inoel-1) * ngl + 3)
               Forca(ino,3) = Forca(ino,3) + Fint(ind,iel)
           enddo

             if ((dUxdx + dUydy + dUzdz) .gt. divUmax) then
                 divUmax = (dUxdx + dUydy + dUzdz)
                 ielDivUMax = iel
             endif
             if ((dUxdx + dUydy + dUzdz) .lt. divUmin) then
                 divUmin = (dUxdx + dUydy + dUzdz)
                 ielDivUMin = iel
             endif
     
enddo

write(22,'(1e12.4E3,a,1e12.4E3,1i6,a,1e12.4E3,1i6)') Time, ' ',-divUmax*xlambda, ielDivUMax, &
' ',-divUmin*xlambda, ielDivUMin

if (N .eq. 1) then
    RefPressure = Pressure(elrefp)
endif

Do iel=1, nume
    Pressure(iel) = Pressure(iel) - RefPressure
    sig11(iel) = sig11(iel) - Pressure(iel)
    sig22(iel) = sig22(iel) - Pressure(iel)
    sig33(iel) = sig33(iel) - Pressure(iel)
Enddo

return
end
