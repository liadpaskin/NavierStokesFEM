     SUBROUTINE SOLVER
	
!     LEITURA DOS DADOS DA MALHA E DEFINICAO DE VARIOS ARRANJOS

      use    modvar
      use    modsolver
      use    modmesh
      use    modloads
      use    modtapes
      use    modtim 


      IMPLICIT NONE

      integer n,kensight,inel,ind
      real*8 tempo_cpu, xtempo, T,FUNC,tempo_cpu1,tempo_cpu2,tempo_cpu3 !, fp1(neq) , fp2(neq)

      tempo_cpu = 0.d0
      kensight = 0
      
      if ((etype .eq. 'NavierStokes3D') .or. (etype .eq. 'NavierStokes2D')) then
           call NavierStokes3D_Solver
      else
          write (*,*) 'Este solver e especifico para NavierStokes3D.&
           Favor checar o etype.'
      endif

      RETURN
    END          
