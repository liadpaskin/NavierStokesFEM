subroutine Calc_Pressure_quad4(Pressure,TFUNC,N,Time)

use    modsolver, only: u,Forca, Fint, sig11, sig12, sig13, sig22, sig23, sig33
use    modloads , only: up
use    modvar   , only: nume, nnoel, numnp, nummat,nd,ncprop,neq,neqp,ngl,RefPressure,elrefp
use    modmesh  , only: incid, id, x, y, prop, mtype, lm

      implicit none
    
!   Global
REAL*8               :: Pressure(nume)
      
!   Local
INTEGER              :: i, iel, inoel, ino, ieq(ngl), jeq,N, ieldivUmax, ieldivUmin, ind


REAL*8               :: x12, x14, y12, y14, DetJ, DetJ_inv, area, xLambda, xMu, TFUNC,areatot, Time
REAL*8               :: r, s, dNx(nnoel), dNy(nnoel), Vel_inoel(ngl), Wpg
REAL*8               :: dUxdx, dUxdy, dUxdz, dUydx, dUydy, dUydz, dUzdx, dUzdy, dUzdz
REAL*8               :: dNdx(nnoel),dNdy(nnoel),dNdz(nnoel)
REAL*8               :: dUx, dUy, divUmax, divUmin
character*70 :: Filename1

FILENAME1 = './OUT1/div.OUT'
if (N .eq. 1) then
    OPEN  (UNIT=22,FILE=FILENAME1, form = 'formatted')  
    write (22,'(a, 1e12.4E3)') 'Tempo, Pressao minima, Elemento de div(u) maximo, &
    Pressao maxima, Elemento de div(u) minimo, Lambda = ', prop(mtype(1),6)
else
    OPEN (UNIT=22,FILE=FILENAME1, position='append', Status= 'unknown') 
endif
    
!rpg(1) = 0.666666666666667d0 ! Coordenadas ponto 1 (2/3; 1/6)
!spg(1) = 0.166666666666667d0 !
!rpg(2) = 0.166666666666667d0 ! Coordenadas ponto 2 (1/6; 2/3)
!spg(2) = 0.666666666666667d0
!rpg(3) = 0.166666666666667d0 ! Coordenadas ponto 3 (1/6; 1/6)
!spg(3) = 0.166666666666667d0
!Wpg    = 0.333333333333333d0 ! Fator de ponderação dos pontos de Gauss (1/3)

divUmax = 0.d0
divUmin = 0.d0
forca(:,:) = 0.d0
     
do iel = 1, nume ! loop nos elementos
         
     xLambda = prop(mtype(iel),6)       ! Coeficiente de incompressibilidade
     xMu     = prop(mtype(iel),5)       ! Viscosidade din�mica
     
       x12 = x(incid(iel,2)) - x(incid(iel,1))  !
       x14 = x(incid(iel,4)) - x(incid(iel,1))  ! Vetores que definem o quadrilatero
       y12 = y(incid(iel,2)) - y(incid(iel,1))  !
       y14 = y(incid(iel,4)) - y(incid(iel,1))  ! Vetores que definem o quadrilatero

       DetJ = x12*y14-x14*y12

       DetJ_inv = 1.d0/DetJ       ! Inverso do determinante da tranformação
       area = DetJ        ! Area do elemento
       areatot = areatot + area   ! Contador: Area total da malha

       
         dUxdx = 0.d0
         dUxdy = 0.d0
         dUydx = 0.d0
         dUydy = 0.d0
     
     ! do i =1, 3  ! loop pontos gauss

             r = 0.5d0 ! rpg(i)
             s = 0.5d0 ! spg(i)

            dNdx(1) = (y14  *(-1+s)-y12 *(-1+r) )*detJ_Inv
            dNdy(1) = (-x14 *(-1+s)+x12 *(-1+r) )*detJ_Inv
            dNdx(2) = (y14  *(1-s) +y12 *r      )*detJ_Inv
            dNdy(2) = (-x14 *(1-s) -x12 *r      )*detJ_Inv
            dNdx(3) = (y14  *s     -y12 *r      )*detJ_Inv
            dNdy(3) = (-x14 *s     +x12 *r      )*detJ_Inv
            dNdx(4) = (-y14 *s     -y12 *(1-r)  )*detJ_Inv
            dNdy(4) = (x14  *s     +x12 *(1-r)  )*detJ_Inv
                                                                                                        !
             do inoel = 1, nnoel  ! Loop nos n�s do elemento
                  ino = incid (iel,inoel)   ! Indice do n� a qual se refere este Loop
                  ieq(1) = id(1,ino) ! Indice da equa��o referente ao primeiro grau de liberdade deste n�
                  ieq(2) = id(2,ino) ! Indice da equa��o referente ao segundo grau de liberdade deste n�
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
                 dUxdx =  dUxdx + dNdx(inoel) * vel_inoel(1)
                 dUxdy =  dUxdy + dNdy(inoel) * vel_inoel(1)
                 dUydx =  dUydx + dNdx(inoel) * vel_inoel(2)
                 dUydy =  dUydy + dNdy(inoel) * vel_inoel(2)
             enddo
     ! enddo
     
             Pressure(iel) = - xLambda * (dUxdx + dUydy) ! * Wpg
             sig11(iel) =  2 * xMu * (dUxdx)
             sig12(iel) =      xMu * (dUxdy + dUydx)
             sig22(iel) =  2 * xMu * (dUydy)


           Fint( 1,iel) = ( dNdx(1) * sig11(iel) + dNdy(1) * sig12(iel) ) * area
           Fint( 2,iel) = ( dNdy(1) * sig22(iel) + dNdx(1) * sig12(iel) ) * area
           Fint( 3,iel) = ( dNdx(2) * sig11(iel) + dNdy(2) * sig12(iel) ) * area
           Fint( 4,iel) = ( dNdy(2) * sig22(iel) + dNdx(2) * sig12(iel) ) * area
           Fint( 5,iel) = ( dNdx(3) * sig11(iel) + dNdy(3) * sig12(iel) ) * area
           Fint( 6,iel) = ( dNdy(3) * sig22(iel) + dNdx(3) * sig12(iel) ) * area
           Fint( 7,iel) = ( dNdx(4) * sig11(iel) + dNdy(4) * sig12(iel) ) * area
           Fint( 8,iel) = ( dNdy(4) * sig22(iel) + dNdx(4) * sig12(iel) ) * area

           do inoel =1, nnoel
               ino = incid(iel,inoel)
               ind = ((inoel-1) * ngl + 1)
               Forca(ino,1) = Forca(ino,1) + Fint(ind,iel)
               ind = ((inoel-1) * ngl + 2)
               Forca(ino,2) = Forca(ino,2) + Fint(ind,iel)
           enddo

             if ((dUxdx + dUydy) .gt. divUmax) then
                 divUmax = (dUxdx + dUydy)
                 ielDivUMax = iel
             endif
             if ((dUxdx + dUydy) .lt. divUmin) then
                 divUmin = (dUxdx + dUydy)
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
Enddo

return
end
