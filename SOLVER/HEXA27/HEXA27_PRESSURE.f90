subroutine Calc_Pressure_Hexa27(Pressure,TFUNC,N,Time)

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
REAL*8               :: Vel_inoel(ngl), r, s, t
REAL*8               :: dUxdx, dUxdy, dUxdz, dUydx, dUydy, dUydz, dUzdx, dUzdy, dUzdz
REAL*8               :: dNdX(nnoel,ngl),dNdR(nnoel,ngl)
REAL*8               :: dUx, dUy, divUmax, divUmin
REAL*8               :: Jinv(3,3)
character*70 :: Filename1

FILENAME1 = './OUT1/div.OUT'
if (N .eq. 1) then
    OPEN  (UNIT=22,FILE=FILENAME1, form = 'formatted')  
    write (22,'(a, 1e12.4E3)') 'Tempo, Pressao minima, Elemento de div(u) maximo, &
    Pressao maxima, Elemento de div(u) minimo, Lambda = ', prop(mtype(1),6)
else
    OPEN (UNIT=22,FILE=FILENAME1, position='append', Status= 'unknown') 
endif

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

       DetJ = (y14*z15-z14*y15)*x12+(x15*z14-x14*z15)*y12+(-x15*y14+x14*y15)*z12     ! Determinante da transformacao em coordenadas naturais

       DetJ_inv = 1.d0/DetJ       ! Inverso do determinante da tranformacao
       vol = DetJ                ! Volume do elemento
       voltot = voltot + vol   ! Contador: volume total da malha

      Jinv(1,1) = 0.1D1 * detJ_inv * (y14 * z15 - z14 * y15)
      Jinv(1,2) = -0.1D1 * detJ_inv * (-y15 * z12 + y12 * z15)
      Jinv(1,3) = 0.1D1 * detJ_inv * (y12 * z14 - y14 * z12)
      Jinv(2,1) = -0.1D1 * detJ_inv * (-x15 * z14 + x14 * z15)
      Jinv(2,2) = 0.1D1 * detJ_inv * (x12 * z15 - x15 * z12)
      Jinv(2,3) = -0.1D1 * detJ_inv * (x12 * z14 - z12 * x14)
      Jinv(3,1) = 0.1D1 * detJ_inv * (-x15 * y14 + x14 * y15)
      Jinv(3,2) = -0.1D1 * detJ_inv * (x12 * y15 - x15 * y12)
      Jinv(3,3) = 0.1D1 * detJ_inv * (x12 * y14 - x14 * y12)
       
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

            dNdR(1,1)=-0.3D1-0.27D2*s*t+0.18D2*t**2*s+0.16D2*s**2*t**2*r-0.24D2*t**2*r*s+0.18D2*s**2*t-0.24D2*s**2*r*t+0.36D2*r*s*t-0.12D2*s**2*t**2-0.12D2*r*s-0.6D1*s**2+0.8D1*t**2*r-0.6D1*t**2+0.9D1*t-0.12D2*t*r+0.9D1*s+0.8D1*s**2*r+0.4D1*r

            dNdR(1,2)=-0.3D1-0.27D2*t*r+0.18D2*t**2*r+0.16D2*r**2*t**2*s-0.12D2*r**2*t**2+0.36D2*r*s*t-0.24D2*r**2*s*t+0.18D2*r**2*t-0.24D2*t**2*r*s-0.6D1*r**2-0.12D2*r*s+0.8D1*t**2*s-0.6D1*t**2+0.9D1*r+0.9D1*t+0.8D1*r**2*s-0.12D2*s*t+0.4D1*s

            dNdR(1,3)=-0.3D1-0.27D2*r*s+0.36D2*r*s*t+0.16D2*r**2*s**2*t-0.24D2*r**2*s*t+0.18D2*s**2*r-0.12D2*r**2*s**2+0.18D2*r**2*s-0.24D2*s**2*r*t+0.8D1*r**2*t-0.12D2*t*r+0.9D1*r+0.8D1*s**2*t-0.6D1*r**2-0.12D2*s*t+0.9D1*s-0.6D1*s**2+0.4D1*t

            dNdR(2,1)=-0.1D1-0.9D1*s*t+0.6D1*t**2*s-0.24D2*t**2*r*s+0.6D1*s**2*t-0.24D2*s**2*r*t+0.36D2*r*s*t-0.4D1*s**2*t**2+0.16D2*s**2*t**2*r-0.12D2*r*s-0.2D1*s**2+0.8D1*t**2*r-0.2D1*t**2+0.3D1*t-0.12D2*t*r+0.3D1*s+0.8D1*s**2*r+0.4D1*r

            dNdR(2,2)=-0.9D1*t*r+0.6D1*t**2*r-0.12D2*r**2*t**2+0.12D2*r*s*t-0.24D2*r**2*s*t+0.18D2*r**2*t-0.8D1*t**2*r*s+0.16D2*r**2*t**2*s-0.6D1*r**2-0.4D1*r*s+0.3D1*r+0.8D1*r**2*s

            dNdR(2,3)=-0.9D1*r*s+0.12D2*r*s*t-0.24D2*r**2*s*t+0.6D1*s**2*r-0.12D2*r**2*s**2+0.18D2*r**2*s-0.8D1*s**2*r*t+0.16D2*r**2*s**2*t+0.8D1*r**2*t-0.4D1*t*r+0.3D1*r-0.6D1*r**2

            dNdR(3,1)=-0.3D1*s*t+0.2D1*t**2*s-0.8D1*t**2*r*s+0.6D1*s**2*t-0.24D2*s**2*r*t+0.12D2*r*s*t-0.4D1*s**2*t**2+0.16D2*s**2*t**2*r-0.4D1*r*s-0.2D1*s**2+0.1D1*s+0.8D1*s**2*r

            dNdR(3,2)=-0.3D1*t*r+0.2D1*t**2*r-0.4D1*r**2*t**2+0.12D2*r*s*t-0.24D2*r**2*s*t+0.6D1*r**2*t-0.8D1*t**2*r*s+0.16D2*r**2*t**2*s-0.2D1*r**2-0.4D1*r*s+0.1D1*r+0.8D1*r**2*s

            dNdR(3,3)=-0.3D1*r*s+0.4D1*r*s*t-0.8D1*r**2*s*t+0.6D1*s**2*r-0.12D2*r**2*s**2+0.6D1*r**2*s-0.8D1*s**2*r*t+0.16D2*r**2*s**2*t

            dNdR(4,1)=-0.9D1*s*t+0.6D1*t**2*s-0.8D1*t**2*r*s+0.18D2*s**2*t-0.24D2*s**2*r*t+0.12D2*r*s*t-0.12D2*s**2*t**2+0.16D2*s**2*t**2*r-0.4D1*r*s-0.6D1*s**2+0.3D1*s+0.8D1*s**2*r

            dNdR(4,2)=-0.1D1-0.9D1*t*r+0.6D1*t**2*r-0.4D1*r**2*t**2+0.36D2*r*s*t-0.24D2*r**2*s*t+0.6D1*r**2*t-0.24D2*t**2*r*s+0.16D2*r**2*t**2*s-0.2D1*r**2-0.12D2*r*s+0.8D1*t**2*s-0.2D1*t**2+0.3D1*r+0.3D1*t+0.8D1*r**2*s-0.12D2*s*t+0.4D1*s

            dNdR(4,3)=-0.9D1*r*s+0.12D2*r*s*t-0.8D1*r**2*s*t+0.18D2*s**2*r-0.12D2*r**2*s**2+0.6D1*r**2*s-0.24D2*s**2*r*t+0.16D2*r**2*s**2*t+0.8D1*s**2*t-0.4D1*s*t+0.3D1*s-0.6D1*s**2

            dNdR(5,1)=-0.9D1*s*t+0.18D2*t**2*s-0.24D2*t**2*r*s+0.6D1*s**2*t-0.8D1*s**2*r*t+0.12D2*r*s*t-0.12D2*s**2*t**2+0.16D2*s**2*t**2*r+0.8D1*t**2*r-0.6D1*t**2+0.3D1*t-0.4D1*t*r

            dNdR(5,2)=-0.9D1*t*r+0.18D2*t**2*r-0.12D2*r**2*t**2+0.12D2*r*s*t-0.8D1*r**2*s*t+0.6D1*r**2*t-0.24D2*t**2*r*s+0.16D2*r**2*t**2*s+0.8D1*t**2*s-0.6D1*t**2+0.3D1*t-0.4D1*s*t

            dNdR(5,3)=-0.1D1-0.9D1*r*s+0.36D2*r*s*t-0.24D2*r**2*s*t+0.6D1*s**2*r-0.4D1*r**2*s**2+0.6D1*r**2*s-0.24D2*s**2*r*t+0.16D2*r**2*s**2*t+0.8D1*r**2*t-0.12D2*t*r+0.3D1*r+0.8D1*s**2*t-0.2D1*r**2-0.12D2*s*t+0.3D1*s-0.2D1*s**2+0.4D1*t

            dNdR(6,1)=-0.3D1*s*t+0.6D1*t**2*s-0.24D2*t**2*r*s+0.2D1*s**2*t-0.8D1*s**2*r*t+0.12D2*r*s*t-0.4D1*s**2*t**2+0.16D2*s**2*t**2*r+0.8D1*t**2*r-0.2D1*t**2+0.1D1*t-0.4D1*t*r

            dNdR(6,2)=-0.3D1*t*r+0.6D1*t**2*r-0.12D2*r**2*t**2+0.4D1*r*s*t-0.8D1*r**2*s*t+0.6D1*r**2*t-0.8D1*t**2*r*s+0.16D2*r**2*t**2*s

            dNdR(6,3)=-0.3D1*r*s+0.12D2*r*s*t-0.24D2*r**2*s*t+0.2D1*s**2*r-0.4D1*r**2*s**2+0.6D1*r**2*s-0.8D1*s**2*r*t+0.16D2*r**2*s**2*t+0.8D1*r**2*t-0.4D1*t*r+0.1D1*r-0.2D1*r**2

            dNdR(7,1)=-0.1D1*s*t+0.2D1*t**2*s-0.8D1*t**2*r*s+0.2D1*s**2*t-0.8D1*s**2*r*t+0.4D1*r*s*t-0.4D1*s**2*t**2+0.16D2*s**2*t**2*r

            dNdR(7,2)=-0.1D1*t*r+0.2D1*t**2*r-0.4D1*r**2*t**2+0.4D1*r*s*t-0.8D1*r**2*s*t+0.2D1*r**2*t-0.8D1*t**2*r*s+0.16D2*r**2*t**2*s

            dNdR(7,3)=-0.1D1*r*s+0.4D1*r*s*t-0.8D1*r**2*s*t+0.2D1*s**2*r-0.4D1*r**2*s**2+0.2D1*r**2*s-0.8D1*s**2*r*t+0.16D2*r**2*s**2*t

            dNdR(8,1)=-0.3D1*s*t+0.6D1*t**2*s-0.8D1*t**2*r*s+0.6D1*s**2*t-0.8D1*s**2*r*t+0.4D1*r*s*t-0.12D2*s**2*t**2+0.16D2*s**2*t**2*r

            dNdR(8,2)=-0.3D1*t*r+0.6D1*t**2*r-0.4D1*r**2*t**2+0.12D2*r*s*t-0.8D1*r**2*s*t+0.2D1*r**2*t-0.24D2*t**2*r*s+0.16D2*r**2*t**2*s+0.8D1*t**2*s-0.2D1*t**2+0.1D1*t-0.4D1*s*t

            dNdR(8,3)=-0.3D1*r*s+0.12D2*r*s*t-0.8D1*r**2*s*t+0.6D1*s**2*r-0.4D1*r**2*s**2+0.2D1*r**2*s-0.24D2*s**2*r*t+0.16D2*r**2*s**2*t+0.8D1*s**2*t-0.4D1*s*t+0.1D1*s-0.2D1*s**2

            dNdR(9,1)=0.36D2*s*t-0.36D2*t**2*s+0.48D2*t**2*r*s-0.24D2*s**2*t+0.32D2*s**2*r*t-0.48D2*r*s*t+0.24D2*s**2*t**2-0.32D2*s**2*t**2*r-0.16D2*t**2*r+0.12D2*t**2-0.12D2*t+0.16D2*t*r

            dNdR(9,2)=0.36D2*t*r-0.36D2*t**2*r+0.24D2*r**2*t**2-0.48D2*r*s*t+0.32D2*r**2*s*t-0.24D2*r**2*t+0.48D2*t**2*r*s-0.32D2*r**2*t**2*s-0.16D2*t**2*s+0.12D2*t**2-0.12D2*t+0.16D2*s*t

            dNdR(9,3)=0.4D1+0.36D2*r*s-0.72D2*r*s*t+0.48D2*r**2*s*t-0.24D2*s**2*r+0.16D2*r**2*s**2-0.24D2*r**2*s+0.48D2*s**2*r*t-0.32D2*r**2*s**2*t-0.16D2*r**2*t+0.24D2*t*r-0.12D2*r-0.16D2*s**2*t+0.8D1*r**2+0.24D2*s*t-0.12D2*s+0.8D1*s**2-0.8D1*t

            dNdR(10,1)=0.12D2*s*t-0.12D2*t**2*s+0.48D2*t**2*r*s-0.8D1*s**2*t+0.32D2*s**2*r*t-0.48D2*r*s*t+0.8D1*s**2*t**2-0.32D2*s**2*t**2*r-0.16D2*t**2*r+0.4D1*t**2-0.4D1*t+0.16D2*t*r

            dNdR(10,2)=0.12D2*t*r-0.12D2*t**2*r+0.24D2*r**2*t**2-0.16D2*r*s*t+0.32D2*r**2*s*t-0.24D2*r**2*t+0.16D2*t**2*r*s-0.32D2*r**2*t**2*s

            dNdR(10,3)=0.12D2*r*s-0.24D2*r*s*t+0.48D2*r**2*s*t-0.8D1*s**2*r+0.16D2*r**2*s**2-0.24D2*r**2*s+0.16D2*s**2*r*t-0.32D2*r**2*s**2*t-0.16D2*r**2*t+0.8D1*t*r-0.4D1*r+0.8D1*r**2

            dNdR(11,1)=0.4D1*s*t-0.4D1*t**2*s+0.16D2*t**2*r*s-0.8D1*s**2*t+0.32D2*s**2*r*t-0.16D2*r*s*t+0.8D1*s**2*t**2-0.32D2*s**2*t**2*r

            dNdR(11,2)=0.4D1*t*r-0.4D1*t**2*r+0.8D1*r**2*t**2-0.16D2*r*s*t+0.32D2*r**2*s*t-0.8D1*r**2*t+0.16D2*t**2*r*s-0.32D2*r**2*t**2*s

            dNdR(11,3)=0.4D1*r*s-0.8D1*r*s*t+0.16D2*r**2*s*t-0.8D1*s**2*r+0.16D2*r**2*s**2-0.8D1*r**2*s+0.16D2*s**2*r*t-0.32D2*r**2*s**2*t

            dNdR(12,1)=0.12D2*s*t-0.12D2*t**2*s+0.16D2*t**2*r*s-0.24D2*s**2*t+0.32D2*s**2*r*t-0.16D2*r*s*t+0.24D2*s**2*t**2-0.32D2*s**2*t**2*r

            dNdR(12,2)=0.12D2*t*r-0.12D2*t**2*r+0.8D1*r**2*t**2-0.48D2*r*s*t+0.32D2*r**2*s*t-0.8D1*r**2*t+0.48D2*t**2*r*s-0.32D2*r**2*t**2*s-0.16D2*t**2*s+0.4D1*t**2-0.4D1*t+0.16D2*s*t

            dNdR(12,3)=0.12D2*r*s-0.24D2*r*s*t+0.16D2*r**2*s*t-0.24D2*s**2*r+0.16D2*r**2*s**2-0.8D1*r**2*s+0.48D2*s**2*r*t-0.32D2*r**2*s**2*t-0.16D2*s**2*t+0.8D1*s*t-0.4D1*s+0.8D1*s**2

            dNdR(13,1)=-0.48D2*s*t+0.32D2*t**2*s-0.64D2*t**2*r*s+0.48D2*s**2*t-0.96D2*s**2*r*t+0.96D2*r*s*t-0.32D2*s**2*t**2+0.64D2*s**2*t**2*r-0.32D2*r*s-0.16D2*s**2+0.16D2*s+0.32D2*s**2*r

            dNdR(13,2)=-0.48D2*t*r+0.32D2*t**2*r-0.32D2*r**2*t**2+0.96D2*r*s*t-0.96D2*r**2*s*t+0.48D2*r**2*t-0.64D2*t**2*r*s+0.64D2*r**2*t**2*s-0.16D2*r**2-0.32D2*r*s+0.16D2*r+0.32D2*r**2*s

            dNdR(13,3)=-0.48D2*r*s+0.64D2*r*s*t-0.64D2*r**2*s*t+0.48D2*s**2*r-0.48D2*r**2*s**2+0.48D2*r**2*s-0.64D2*s**2*r*t+0.64D2*r**2*s**2*t

            dNdR(14,1)=0.64D2*s*t-0.64D2*t**2*s+0.128D3*t**2*r*s-0.64D2*s**2*t+0.128D3*s**2*r*t-0.128D3*r*s*t+0.64D2*s**2*t**2-0.128D3*s**2*t**2*r

            dNdR(14,2)=0.64D2*t*r-0.64D2*t**2*r+0.64D2*r**2*t**2-0.128D3*r*s*t+0.128D3*r**2*s*t-0.64D2*r**2*t+0.128D3*t**2*r*s-0.128D3*r**2*t**2*s

            dNdR(14,3)=0.64D2*r*s-0.128D3*r*s*t+0.128D3*r**2*s*t-0.64D2*s**2*r+0.64D2*r**2*s**2-0.64D2*r**2*s+0.128D3*s**2*r*t-0.128D3*r**2*s**2*t

            dNdR(15,1)=-0.16D2*s*t+0.32D2*t**2*s-0.64D2*t**2*r*s+0.16D2*s**2*t-0.32D2*s**2*r*t+0.32D2*r*s*t-0.32D2*s**2*t**2+0.64D2*s**2*t**2*r

            dNdR(15,2)=-0.16D2*t*r+0.32D2*t**2*r-0.32D2*r**2*t**2+0.32D2*r*s*t-0.32D2*r**2*s*t+0.16D2*r**2*t-0.64D2*t**2*r*s+0.64D2*r**2*t**2*s

            dNdR(15,3)=-0.16D2*r*s+0.64D2*r*s*t-0.64D2*r**2*s*t+0.16D2*s**2*r-0.16D2*r**2*s**2+0.16D2*r**2*s-0.64D2*s**2*r*t+0.64D2*r**2*s**2*t

            dNdR(16,1)=0.4D1+0.36D2*s*t-0.24D2*t**2*s+0.48D2*t**2*r*s-0.24D2*s**2*t+0.48D2*s**2*r*t-0.72D2*r*s*t+0.16D2*s**2*t**2-0.32D2*s**2*t**2*r+0.24D2*r*s+0.8D1*s**2-0.16D2*t**2*r+0.8D1*t**2-0.12D2*t+0.24D2*t*r-0.12D2*s-0.16D2*s**2*r-0.8D1*r

            dNdR(16,2)=0.36D2*t*r-0.24D2*t**2*r+0.24D2*r**2*t**2-0.48D2*r*s*t+0.48D2*r**2*s*t-0.36D2*r**2*t+0.32D2*t**2*r*s-0.32D2*r**2*t**2*s+0.12D2*r**2+0.16D2*r*s-0.12D2*r-0.16D2*r**2*s

            dNdR(16,3)=0.36D2*r*s-0.48D2*r*s*t+0.48D2*r**2*s*t-0.24D2*s**2*r+0.24D2*r**2*s**2-0.36D2*r**2*s+0.32D2*s**2*r*t-0.32D2*r**2*s**2*t-0.16D2*r**2*t+0.16D2*t*r-0.12D2*r+0.12D2*r**2

            dNdR(17,1)=0.12D2*s*t-0.8D1*t**2*s+0.32D2*t**2*r*s-0.12D2*s**2*t+0.48D2*s**2*r*t-0.48D2*r*s*t+0.8D1*s**2*t**2-0.32D2*s**2*t**2*r+0.16D2*r*s+0.4D1*s**2-0.4D1*s-0.16D2*s**2*r

            dNdR(17,2)=0.12D2*t*r-0.8D1*t**2*r+0.16D2*r**2*t**2-0.24D2*r*s*t+0.48D2*r**2*s*t-0.24D2*r**2*t+0.16D2*t**2*r*s-0.32D2*r**2*t**2*s+0.8D1*r**2+0.8D1*r*s-0.4D1*r-0.16D2*r**2*s

            dNdR(17,3)=0.12D2*r*s-0.16D2*r*s*t+0.32D2*r**2*s*t-0.12D2*s**2*r+0.24D2*r**2*s**2-0.24D2*r**2*s+0.16D2*s**2*r*t-0.32D2*r**2*s**2*t

            dNdR(18,1)=0.12D2*s*t-0.8D1*t**2*s+0.16D2*t**2*r*s-0.24D2*s**2*t+0.48D2*s**2*r*t-0.24D2*r*s*t+0.16D2*s**2*t**2-0.32D2*s**2*t**2*r+0.8D1*r*s+0.8D1*s**2-0.4D1*s-0.16D2*s**2*r

            dNdR(18,2)=0.12D2*t*r-0.8D1*t**2*r+0.8D1*r**2*t**2-0.48D2*r*s*t+0.48D2*r**2*s*t-0.12D2*r**2*t+0.32D2*t**2*r*s-0.32D2*r**2*t**2*s+0.4D1*r**2+0.16D2*r*s-0.4D1*r-0.16D2*r**2*s

            dNdR(18,3)=0.12D2*r*s-0.16D2*r*s*t+0.16D2*r**2*s*t-0.24D2*s**2*r+0.24D2*r**2*s**2-0.12D2*r**2*s+0.32D2*s**2*r*t-0.32D2*r**2*s**2*t

            dNdR(19,1)=0.36D2*s*t-0.24D2*t**2*s+0.32D2*t**2*r*s-0.36D2*s**2*t+0.48D2*s**2*r*t-0.48D2*r*s*t+0.24D2*s**2*t**2-0.32D2*s**2*t**2*r+0.16D2*r*s+0.12D2*s**2-0.12D2*s-0.16D2*s**2*r

            dNdR(19,2)=0.4D1+0.36D2*t*r-0.24D2*t**2*r+0.16D2*r**2*t**2-0.72D2*r*s*t+0.48D2*r**2*s*t-0.24D2*r**2*t+0.48D2*t**2*r*s-0.32D2*r**2*t**2*s+0.8D1*r**2+0.24D2*r*s-0.16D2*t**2*s+0.8D1*t**2-0.12D2*r-0.12D2*t-0.16D2*r**2*s+0.24D2*s*t-0.8D1*s

            dNdR(19,3)=0.36D2*r*s-0.48D2*r*s*t+0.32D2*r**2*s*t-0.36D2*s**2*r+0.24D2*r**2*s**2-0.24D2*r**2*s+0.48D2*s**2*r*t-0.32D2*r**2*s**2*t-0.16D2*s**2*t+0.16D2*s*t-0.12D2*s+0.12D2*s**2

            dNdR(20,1)=-0.48D2*s*t+0.48D2*t**2*s-0.96D2*t**2*r*s+0.32D2*s**2*t-0.64D2*s**2*r*t+0.96D2*r*s*t-0.32D2*s**2*t**2+0.64D2*s**2*t**2*r+0.32D2*t**2*r-0.16D2*t**2+0.16D2*t-0.32D2*t*r

            dNdR(20,2)=-0.48D2*t*r+0.48D2*t**2*r-0.48D2*r**2*t**2+0.64D2*r*s*t-0.64D2*r**2*s*t+0.48D2*r**2*t-0.64D2*t**2*r*s+0.64D2*r**2*t**2*s

            dNdR(20,3)=-0.48D2*r*s+0.96D2*r*s*t-0.96D2*r**2*s*t+0.32D2*s**2*r-0.32D2*r**2*s**2+0.48D2*r**2*s-0.64D2*s**2*r*t+0.64D2*r**2*s**2*t+0.32D2*r**2*t-0.32D2*t*r+0.16D2*r-0.16D2*r**2

            dNdR(21,1)=-0.16D2*s*t+0.16D2*t**2*s-0.64D2*t**2*r*s+0.16D2*s**2*t-0.64D2*s**2*r*t+0.64D2*r*s*t-0.16D2*s**2*t**2+0.64D2*s**2*t**2*r

            dNdR(21,2)=-0.16D2*t*r+0.16D2*t**2*r-0.32D2*r**2*t**2+0.32D2*r*s*t-0.64D2*r**2*s*t+0.32D2*r**2*t-0.32D2*t**2*r*s+0.64D2*r**2*t**2*s

            dNdR(21,3)=-0.16D2*r*s+0.32D2*r*s*t-0.64D2*r**2*s*t+0.16D2*s**2*r-0.32D2*r**2*s**2+0.32D2*r**2*s-0.32D2*s**2*r*t+0.64D2*r**2*s**2*t

            dNdR(22,1)=-0.16D2*s*t+0.16D2*t**2*s-0.32D2*t**2*r*s+0.32D2*s**2*t-0.64D2*s**2*r*t+0.32D2*r*s*t-0.32D2*s**2*t**2+0.64D2*s**2*t**2*r

            dNdR(22,2)=-0.16D2*t*r+0.16D2*t**2*r-0.16D2*r**2*t**2+0.64D2*r*s*t-0.64D2*r**2*s*t+0.16D2*r**2*t-0.64D2*t**2*r*s+0.64D2*r**2*t**2*s

            dNdR(22,3)=-0.16D2*r*s+0.32D2*r*s*t-0.32D2*r**2*s*t+0.32D2*s**2*r-0.32D2*r**2*s**2+0.16D2*r**2*s-0.64D2*s**2*r*t+0.64D2*r**2*s**2*t

            dNdR(23,1)=-0.48D2*s*t+0.48D2*t**2*s-0.64D2*t**2*r*s+0.48D2*s**2*t-0.64D2*s**2*r*t+0.64D2*r*s*t-0.48D2*s**2*t**2+0.64D2*s**2*t**2*r

            dNdR(23,2)=-0.48D2*t*r+0.48D2*t**2*r-0.32D2*r**2*t**2+0.96D2*r*s*t-0.64D2*r**2*s*t+0.32D2*r**2*t-0.96D2*t**2*r*s+0.64D2*r**2*t**2*s+0.32D2*t**2*s-0.16D2*t**2+0.16D2*t-0.32D2*s*t

            dNdR(23,3)=-0.48D2*r*s+0.96D2*r*s*t-0.64D2*r**2*s*t+0.48D2*s**2*r-0.32D2*r**2*s**2+0.32D2*r**2*s-0.96D2*s**2*r*t+0.64D2*r**2*s**2*t+0.32D2*s**2*t-0.32D2*s*t+0.16D2*s-0.16D2*s**2

            dNdR(24,1)=0.12D2*s*t-0.24D2*t**2*s+0.48D2*t**2*r*s-0.8D1*s**2*t+0.16D2*s**2*r*t-0.24D2*r*s*t+0.16D2*s**2*t**2-0.32D2*s**2*t**2*r-0.16D2*t**2*r+0.8D1*t**2-0.4D1*t+0.8D1*t*r

            dNdR(24,2)=0.12D2*t*r-0.24D2*t**2*r+0.24D2*r**2*t**2-0.16D2*r*s*t+0.16D2*r**2*s*t-0.12D2*r**2*t+0.32D2*t**2*r*s-0.32D2*r**2*t**2*s

            dNdR(24,3)=0.12D2*r*s-0.48D2*r*s*t+0.48D2*r**2*s*t-0.8D1*s**2*r+0.8D1*r**2*s**2-0.12D2*r**2*s+0.32D2*s**2*r*t-0.32D2*r**2*s**2*t-0.16D2*r**2*t+0.16D2*t*r-0.4D1*r+0.4D1*r**2

            dNdR(25,1)=0.4D1*s*t-0.8D1*t**2*s+0.32D2*t**2*r*s-0.4D1*s**2*t+0.16D2*s**2*r*t-0.16D2*r*s*t+0.8D1*s**2*t**2-0.32D2*s**2*t**2*r

            dNdR(25,2)=0.4D1*t*r-0.8D1*t**2*r+0.16D2*r**2*t**2-0.8D1*r*s*t+0.16D2*r**2*s*t-0.8D1*r**2*t+0.16D2*t**2*r*s-0.32D2*r**2*t**2*s

            dNdR(25,3)=0.4D1*r*s-0.16D2*r*s*t+0.32D2*r**2*s*t-0.4D1*s**2*r+0.8D1*r**2*s**2-0.8D1*r**2*s+0.16D2*s**2*r*t-0.32D2*r**2*s**2*t

            dNdR(26,1)=0.4D1*s*t-0.8D1*t**2*s+0.16D2*t**2*r*s-0.8D1*s**2*t+0.16D2*s**2*r*t-0.8D1*r*s*t+0.16D2*s**2*t**2-0.32D2*s**2*t**2*r

            dNdR(26,2)=0.4D1*t*r-0.8D1*t**2*r+0.8D1*r**2*t**2-0.16D2*r*s*t+0.16D2*r**2*s*t-0.4D1*r**2*t+0.32D2*t**2*r*s-0.32D2*r**2*t**2*s

            dNdR(26,3)=0.4D1*r*s-0.16D2*r*s*t+0.16D2*r**2*s*t-0.8D1*s**2*r+0.8D1*r**2*s**2-0.4D1*r**2*s+0.32D2*s**2*r*t-0.32D2*r**2*s**2*t

            dNdR(27,1)=0.12D2*s*t-0.24D2*t**2*s+0.32D2*t**2*r*s-0.12D2*s**2*t+0.16D2*s**2*r*t-0.16D2*r*s*t+0.24D2*s**2*t**2-0.32D2*s**2*t**2*r

            dNdR(27,2)=0.12D2*t*r-0.24D2*t**2*r+0.16D2*r**2*t**2-0.24D2*r*s*t+0.16D2*r**2*s*t-0.8D1*r**2*t+0.48D2*t**2*r*s-0.32D2*r**2*t**2*s-0.16D2*t**2*s+0.8D1*t**2-0.4D1*t+0.8D1*s*t

            dNdR(27,3)=0.12D2*r*s-0.48D2*r*s*t+0.32D2*r**2*s*t-0.12D2*s**2*r+0.8D1*r**2*s**2-0.8D1*r**2*s+0.48D2*s**2*r*t-0.32D2*r**2*s**2*t-0.16D2*s**2*t+0.16D2*s*t-0.4D1*s+0.4D1*s**2

              dNdX(1,1) = Jinv(1,1) * dNdR(1,1) + Jinv(1,2) * dNdR(1,2) + Jinv(1,3) * dNdR(1,3)
              dNdX(1,2) = Jinv(2,1) * dNdR(1,1) + Jinv(2,2) * dNdR(1,2) + Jinv(2,3) * dNdR(1,3)
              dNdX(1,3) = Jinv(3,1) * dNdR(1,1) + Jinv(3,2) * dNdR(1,2) + Jinv(3,3) * dNdR(1,3)

              dNdX(2,1) = Jinv(1,1) * dNdR(2,1) + Jinv(1,2) * dNdR(2,2) + Jinv(1,3) * dNdR(2,3)
              dNdX(2,2) = Jinv(2,1) * dNdR(2,1) + Jinv(2,2) * dNdR(2,2) + Jinv(2,3) * dNdR(2,3)
              dNdX(2,3) = Jinv(3,1) * dNdR(2,1) + Jinv(3,2) * dNdR(2,2) + Jinv(3,3) * dNdR(2,3)

              dNdX(3,1) = Jinv(1,1) * dNdR(3,1) + Jinv(1,2) * dNdR(3,2) + Jinv(1,3) * dNdR(3,3)
              dNdX(3,2) = Jinv(2,1) * dNdR(3,1) + Jinv(2,2) * dNdR(3,2) + Jinv(2,3) * dNdR(3,3)
              dNdX(3,3) = Jinv(3,1) * dNdR(3,1) + Jinv(3,2) * dNdR(3,2) + Jinv(3,3) * dNdR(3,3)

              dNdX(4,1) = Jinv(1,1) * dNdR(4,1) + Jinv(1,2) * dNdR(4,2) + Jinv(1,3) * dNdR(4,3)
              dNdX(4,2) = Jinv(2,1) * dNdR(4,1) + Jinv(2,2) * dNdR(4,2) + Jinv(2,3) * dNdR(4,3)
              dNdX(4,3) = Jinv(3,1) * dNdR(4,1) + Jinv(3,2) * dNdR(4,2) + Jinv(3,3) * dNdR(4,3)

              dNdX(5,1) = Jinv(1,1) * dNdR(5,1) + Jinv(1,2) * dNdR(5,2) + Jinv(1,3) * dNdR(5,3)
              dNdX(5,2) = Jinv(2,1) * dNdR(5,1) + Jinv(2,2) * dNdR(5,2) + Jinv(2,3) * dNdR(5,3)
              dNdX(5,3) = Jinv(3,1) * dNdR(5,1) + Jinv(3,2) * dNdR(5,2) + Jinv(3,3) * dNdR(5,3)

              dNdX(6,1) = Jinv(1,1) * dNdR(6,1) + Jinv(1,2) * dNdR(6,2) + Jinv(1,3) * dNdR(6,3)
              dNdX(6,2) = Jinv(2,1) * dNdR(6,1) + Jinv(2,2) * dNdR(6,2) + Jinv(2,3) * dNdR(6,3)
              dNdX(6,3) = Jinv(3,1) * dNdR(6,1) + Jinv(3,2) * dNdR(6,2) + Jinv(3,3) * dNdR(6,3)

              dNdX(7,1) = Jinv(1,1) * dNdR(7,1) + Jinv(1,2) * dNdR(7,2) + Jinv(1,3) * dNdR(7,3)
              dNdX(7,2) = Jinv(2,1) * dNdR(7,1) + Jinv(2,2) * dNdR(7,2) + Jinv(2,3) * dNdR(7,3)
              dNdX(7,3) = Jinv(3,1) * dNdR(7,1) + Jinv(3,2) * dNdR(7,2) + Jinv(3,3) * dNdR(7,3)

              dNdX(8,1) = Jinv(1,1) * dNdR(8,1) + Jinv(1,2) * dNdR(8,2) + Jinv(1,3) * dNdR(8,3)
              dNdX(8,2) = Jinv(2,1) * dNdR(8,1) + Jinv(2,2) * dNdR(8,2) + Jinv(2,3) * dNdR(8,3)
              dNdX(8,3) = Jinv(3,1) * dNdR(8,1) + Jinv(3,2) * dNdR(8,2) + Jinv(3,3) * dNdR(8,3)

              dNdX(9 ,1) = Jinv(1,1) * dNdR(9 ,1) + Jinv(1,2) * dNdR(9 ,2) + Jinv(1,3) * dNdR(9 ,3)
              dNdX(9 ,2) = Jinv(2,1) * dNdR(9 ,1) + Jinv(2,2) * dNdR(9 ,2) + Jinv(2,3) * dNdR(9 ,3)
              dNdX(9 ,3) = Jinv(3,1) * dNdR(9 ,1) + Jinv(3,2) * dNdR(9 ,2) + Jinv(3,3) * dNdR(9 ,3)

              dNdX(10,1) = Jinv(1,1) * dNdR(10,1) + Jinv(1,2) * dNdR(10,2) + Jinv(1,3) * dNdR(10,3)
              dNdX(10,2) = Jinv(2,1) * dNdR(10,1) + Jinv(2,2) * dNdR(10,2) + Jinv(2,3) * dNdR(10,3)
              dNdX(10,3) = Jinv(3,1) * dNdR(10,1) + Jinv(3,2) * dNdR(10,2) + Jinv(3,3) * dNdR(10,3)

              dNdX(11,1) = Jinv(1,1) * dNdR(11,1) + Jinv(1,2) * dNdR(11,2) + Jinv(1,3) * dNdR(11,3)
              dNdX(11,2) = Jinv(2,1) * dNdR(11,1) + Jinv(2,2) * dNdR(11,2) + Jinv(2,3) * dNdR(11,3)
              dNdX(11,3) = Jinv(3,1) * dNdR(11,1) + Jinv(3,2) * dNdR(11,2) + Jinv(3,3) * dNdR(11,3)

              dNdX(12,1) = Jinv(1,1) * dNdR(12,1) + Jinv(1,2) * dNdR(12,2) + Jinv(1,3) * dNdR(12,3)
              dNdX(12,2) = Jinv(2,1) * dNdR(12,1) + Jinv(2,2) * dNdR(12,2) + Jinv(2,3) * dNdR(12,3)
              dNdX(12,3) = Jinv(3,1) * dNdR(12,1) + Jinv(3,2) * dNdR(12,2) + Jinv(3,3) * dNdR(12,3)

              dNdX(13,1) = Jinv(1,1) * dNdR(13,1) + Jinv(1,2) * dNdR(13,2) + Jinv(1,3) * dNdR(13,3)
              dNdX(13,2) = Jinv(2,1) * dNdR(13,1) + Jinv(2,2) * dNdR(13,2) + Jinv(2,3) * dNdR(13,3)
              dNdX(13,3) = Jinv(3,1) * dNdR(13,1) + Jinv(3,2) * dNdR(13,2) + Jinv(3,3) * dNdR(13,3)

              dNdX(14,1) = Jinv(1,1) * dNdR(14,1) + Jinv(1,2) * dNdR(14,2) + Jinv(1,3) * dNdR(14,3)
              dNdX(14,2) = Jinv(2,1) * dNdR(14,1) + Jinv(2,2) * dNdR(14,2) + Jinv(2,3) * dNdR(14,3)
              dNdX(14,3) = Jinv(3,1) * dNdR(14,1) + Jinv(3,2) * dNdR(14,2) + Jinv(3,3) * dNdR(14,3)

              dNdX(15,1) = Jinv(1,1) * dNdR(15,1) + Jinv(1,2) * dNdR(15,2) + Jinv(1,3) * dNdR(15,3)
              dNdX(15,2) = Jinv(2,1) * dNdR(15,1) + Jinv(2,2) * dNdR(15,2) + Jinv(2,3) * dNdR(15,3)
              dNdX(15,3) = Jinv(3,1) * dNdR(15,1) + Jinv(3,2) * dNdR(15,2) + Jinv(3,3) * dNdR(15,3)

              dNdX(16,1) = Jinv(1,1) * dNdR(16,1) + Jinv(1,2) * dNdR(16,2) + Jinv(1,3) * dNdR(16,3)
              dNdX(16,2) = Jinv(2,1) * dNdR(16,1) + Jinv(2,2) * dNdR(16,2) + Jinv(2,3) * dNdR(16,3)
              dNdX(16,3) = Jinv(3,1) * dNdR(16,1) + Jinv(3,2) * dNdR(16,2) + Jinv(3,3) * dNdR(16,3)

              dNdX(17,1) = Jinv(1,1) * dNdR(17,1) + Jinv(1,2) * dNdR(17,2) + Jinv(1,3) * dNdR(17,3)
              dNdX(17,2) = Jinv(2,1) * dNdR(17,1) + Jinv(2,2) * dNdR(17,2) + Jinv(2,3) * dNdR(17,3)
              dNdX(17,3) = Jinv(3,1) * dNdR(17,1) + Jinv(3,2) * dNdR(17,2) + Jinv(3,3) * dNdR(17,3)

              dNdX(18,1) = Jinv(1,1) * dNdR(18,1) + Jinv(1,2) * dNdR(18,2) + Jinv(1,3) * dNdR(18,3)
              dNdX(18,2) = Jinv(2,1) * dNdR(18,1) + Jinv(2,2) * dNdR(18,2) + Jinv(2,3) * dNdR(18,3)
              dNdX(18,3) = Jinv(3,1) * dNdR(18,1) + Jinv(3,2) * dNdR(18,2) + Jinv(3,3) * dNdR(18,3)

              dNdX(19,1) = Jinv(1,1) * dNdR(19,1) + Jinv(1,2) * dNdR(19,2) + Jinv(1,3) * dNdR(19,3)
              dNdX(19,2) = Jinv(2,1) * dNdR(19,1) + Jinv(2,2) * dNdR(19,2) + Jinv(2,3) * dNdR(19,3)
              dNdX(19,3) = Jinv(3,1) * dNdR(19,1) + Jinv(3,2) * dNdR(19,2) + Jinv(3,3) * dNdR(19,3)

              dNdX(20,1) = Jinv(1,1) * dNdR(20,1) + Jinv(1,2) * dNdR(20,2) + Jinv(1,3) * dNdR(20,3)
              dNdX(20,2) = Jinv(2,1) * dNdR(20,1) + Jinv(2,2) * dNdR(20,2) + Jinv(2,3) * dNdR(20,3)
              dNdX(20,3) = Jinv(3,1) * dNdR(20,1) + Jinv(3,2) * dNdR(20,2) + Jinv(3,3) * dNdR(20,3)

              dNdX(21,1) = Jinv(1,1) * dNdR(21,1) + Jinv(1,2) * dNdR(21,2) + Jinv(1,3) * dNdR(21,3)
              dNdX(21,2) = Jinv(2,1) * dNdR(21,1) + Jinv(2,2) * dNdR(21,2) + Jinv(2,3) * dNdR(21,3)
              dNdX(21,3) = Jinv(3,1) * dNdR(21,1) + Jinv(3,2) * dNdR(21,2) + Jinv(3,3) * dNdR(21,3)

              dNdX(22,1) = Jinv(1,1) * dNdR(22,1) + Jinv(1,2) * dNdR(22,2) + Jinv(1,3) * dNdR(22,3)
              dNdX(22,2) = Jinv(2,1) * dNdR(22,1) + Jinv(2,2) * dNdR(22,2) + Jinv(2,3) * dNdR(22,3)
              dNdX(22,3) = Jinv(3,1) * dNdR(22,1) + Jinv(3,2) * dNdR(22,2) + Jinv(3,3) * dNdR(22,3)

              dNdX(23,1) = Jinv(1,1) * dNdR(23,1) + Jinv(1,2) * dNdR(23,2) + Jinv(1,3) * dNdR(23,3)
              dNdX(23,2) = Jinv(2,1) * dNdR(23,1) + Jinv(2,2) * dNdR(23,2) + Jinv(2,3) * dNdR(23,3)
              dNdX(23,3) = Jinv(3,1) * dNdR(23,1) + Jinv(3,2) * dNdR(23,2) + Jinv(3,3) * dNdR(23,3)

              dNdX(24,1) = Jinv(1,1) * dNdR(24,1) + Jinv(1,2) * dNdR(24,2) + Jinv(1,3) * dNdR(24,3)
              dNdX(24,2) = Jinv(2,1) * dNdR(24,1) + Jinv(2,2) * dNdR(24,2) + Jinv(2,3) * dNdR(24,3)
              dNdX(24,3) = Jinv(3,1) * dNdR(24,1) + Jinv(3,2) * dNdR(24,2) + Jinv(3,3) * dNdR(24,3)

              dNdX(25,1) = Jinv(1,1) * dNdR(25,1) + Jinv(1,2) * dNdR(25,2) + Jinv(1,3) * dNdR(25,3)
              dNdX(25,2) = Jinv(2,1) * dNdR(25,1) + Jinv(2,2) * dNdR(25,2) + Jinv(2,3) * dNdR(25,3)
              dNdX(25,3) = Jinv(3,1) * dNdR(25,1) + Jinv(3,2) * dNdR(25,2) + Jinv(3,3) * dNdR(25,3)

              dNdX(26,1) = Jinv(1,1) * dNdR(26,1) + Jinv(1,2) * dNdR(26,2) + Jinv(1,3) * dNdR(26,3)
              dNdX(26,2) = Jinv(2,1) * dNdR(26,1) + Jinv(2,2) * dNdR(26,2) + Jinv(2,3) * dNdR(26,3)
              dNdX(26,3) = Jinv(3,1) * dNdR(26,1) + Jinv(3,2) * dNdR(26,2) + Jinv(3,3) * dNdR(26,3)

              dNdX(27,1) = Jinv(1,1) * dNdR(27,1) + Jinv(1,2) * dNdR(27,2) + Jinv(1,3) * dNdR(27,3)
              dNdX(27,2) = Jinv(2,1) * dNdR(27,1) + Jinv(2,2) * dNdR(27,2) + Jinv(2,3) * dNdR(27,3)
              dNdX(27,3) = Jinv(3,1) * dNdR(27,1) + Jinv(3,2) * dNdR(27,2) + Jinv(3,3) * dNdR(27,3)
                                                                                                          !
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
                 dUxdx =  dUxdx + dNdX(inoel,1) * vel_inoel(1)
                 dUxdy =  dUxdy + dNdX(inoel,2) * vel_inoel(1)
                 dUxdz =  dUxdz + dNdX(inoel,3) * vel_inoel(1)
                 dUydx =  dUydx + dNdX(inoel,1) * vel_inoel(2)
                 dUydy =  dUydy + dNdX(inoel,2) * vel_inoel(2)
                 dUydz =  dUydz + dNdX(inoel,3) * vel_inoel(2)
                 dUzdx =  dUzdx + dNdX(inoel,1) * vel_inoel(3)
                 dUzdy =  dUzdy + dNdX(inoel,2) * vel_inoel(3)
                 dUzdz =  dUzdz + dNdX(inoel,3) * vel_inoel(3)
             enddo
     ! enddo
     
             Pressure(iel) = - xLambda * (dUxdx + dUydy + dUzdz) ! * Wpg
             sig11(iel) =  2 * xMu * (dUxdx)
             sig12(iel) =      xMu * (dUxdy + dUydx)
             sig13(iel) =      xMu * (dUxdz + dUzdx)
             sig22(iel) =  2 * xMu * (dUydy)
             sig23(iel) =      xMu * (dUydz + dUzdy)
             sig33(iel) =  2 * xMu * (dUzdz)


            do inoel = 1, nnoel
                ind = ((inoel-1) * ngl + 1)
                Fint(ind ,iel)   = ( dNdX(inoel,1) * sig11(iel) + dNdX(inoel,2) * sig12(iel) + dNdX(inoel,3) * sig13(iel) ) * vol
                Fint(ind+1 ,iel) = ( dNdX(inoel,1) * sig12(iel) + dNdX(inoel,2) * sig22(iel) + dNdX(inoel,3) * sig23(iel) ) * vol
                Fint(ind+2 ,iel) = ( dNdX(inoel,1) * sig13(iel) + dNdX(inoel,2) * sig23(iel) + dNdX(inoel,3) * sig33(iel) ) * vol
            enddo

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
