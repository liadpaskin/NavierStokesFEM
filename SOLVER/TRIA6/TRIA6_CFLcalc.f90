subroutine tria6_CFLcalc (N, T, passoutput, dT)

 use    modvar,    only : nume, nnoel, ngl
 use    modsolver, only : u, elPe
 use    modmesh,   only : x, y, incid, id, prop, mtype
 use    modloads,  only : up   
 use    modtim,    only : CFLmax, nflag
    
IMPLICIT NONE

real*8       :: CFLmin, Wpg, x23, x13, y23, y13, x12, y12, uCFL(nnoel,ngl),&
 velCFL, largCFL, u1, u2, CFL, T, xMu, xRho, Pe, PeMax, dT

integer      :: iel, inoel, ino, ingl, ieq, N, ielCFLmax, ielCFLmin, ielPeMax, passoutput
character*70 :: Filename1

FILENAME1 = './OUT1/CFL.OUT'
if (N .eq. 1) then
    OPEN  (UNIT=22,FILE=FILENAME1, form = 'formatted')  
    write (22,'(a)') 'Tempo, CFL maximo, Elemento de CFL maximo, Pecle maximo, Elemento de Pecle maximo,'
else
    OPEN (UNIT=22,FILE=FILENAME1, position='append', Status= 'unknown') 
endif

CFLmax = 0.d0
PeMax = 0.d0

do iel = 1, nume ! loop nos elementos       

       x13 = x(incid(iel,1)) - x(incid(iel,3))  !
       x23 = x(incid(iel,2)) - x(incid(iel,3))  ! Vetores que definem o elemento triangular
       x12 = x(incid(iel,1)) - x(incid(iel,2))  !
       y13 = y(incid(iel,1)) - y(incid(iel,3))  !
       y23 = y(incid(iel,2)) - y(incid(iel,3))  !
       y12 = y(incid(iel,1)) - y(incid(iel,2))  !

        do inoel = 1, nnoel                              !
            ino  = INCID(iel,inoel)                  !
            do ingl = 1, ngl                           !
                ieq = ID(ingl,ino)                   !
       	        if (ieq.gt.0) then                   !  Armazena velocidades do passo anterior
       	           uCFL(inoel,ingl) = u(ieq)         !
                else                                 !
                   uCFL(inoel,ingl) =  up(-ieq)      !
                endif                                !
            enddo                                    !
        enddo                                        !

        u1 = 0d0
        u2 = 0d0
        Wpg    = 0.333333333333333d0
        do inoel = 1, nnoel
                u1 = u1 + uCFL(inoel,1)  ! Calcula velocidades
                u2 = u2 + uCFL(inoel,2) !  m�dias
        enddo
        u1 = u1 * Wpg ! Calcula velocidades
        u2 = u2 * Wpg !  m�dias
        
        velCFL = ((u1*u1+u2*u2)**0.5d0)
        largCFL = (((x13*x13+y13*y13)**0.5d0) + ((x23*x23+y23*y23)**0.5d0) + ((x12*x12+y12*y12)**0.5d0)) * Wpg

        CFL = velCFL * Dt / largCFL
        
        if (CFL .gt. CFLmax) then
            CFLmax = CFL
            ielCFLmax = iel
        endif
        
        xRho    = prop(mtype(iel),4)       ! Densidade
        xMu     = prop(mtype(iel),5)       ! Viscosidade din�mica
        Pe  = velCFL * largCFL / (2*xMu/xRho)
        
        elPe(iel) = Pe
        
        if (Pe .gt. Pemax) then
            Pemax = Pe
            ielPEmax = iel
        endif
        
enddo ! end loop nos elementos

IF ((MOD(N-1,nflag) .EQ. 0) .and. (passoutput .ne. 1))  THEN    
    write(22,'(2e12.4E3,1i6,1e12.4E3,1i6)') T, CFLmax, ielCFLmax, Pemax, ielPemax
    write(*   ,'(a,1e12.4E3)') ' CFL   maximo = ' , CFLmax
    write(*   ,'(a,1e12.4E3)') ' Pecle maximo = ' , Pemax
ENDIF

CLOSE (UNIT=22) 

return

end subroutine
