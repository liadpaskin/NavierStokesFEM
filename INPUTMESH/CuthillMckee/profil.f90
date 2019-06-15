      subroutine profil (jdiag, ix, nnomax, nnos, nelm, nad ,nume1)
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
! .                                                                 .
! .                                                    20/05/88     .
! .   - descricao :                                                 .
! .                                                                 .
! .     determina o perfil da matriz de rigidez.                    .
! .                                                                 .
! .                                                                 .
! .   - parametros de entrada :                                     .
! .                                                                 .
! .     nnos = numero de nos da malha                               .
! .     nelm = numero de elementos                                  .
! .     ix   = incidencia dos elementos                             .
! .     nnomax = numero max. de nos por elemento                    .
! .                                                                 .
! .   - parametros de saida :                                       .
! .                                                                 .
! .     jdiag = vetor apontador do perfil da matriz                 .
! .     nad  = numero de coeficientes do perfil da matriz           .
! .                                                                 .
! .   - chamada por:                                                .
! .                                                                 .
! .                                                                 .
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
    
      dimension jdiag(nnos), ix(nnomax,nelm)


!.... compute column heights
      do 300 n = nume1, nelm
         do 200 i = 1, nnomax
           ii = ix(i,n)
           if(ii.eq.0) go to 200
           do 100 j = i, nnomax
              jj = ix(j,n)
              if(jj.eq.0) go to 100
              m = max0 (ii,jj)
              jdiag(m) = max0 (jdiag(m),iabs(ii-jj))
  100      continue
  200    continue
  300 continue
!.... compute diagonal pointers for profile
      nad = 1
      jdiag(1) = 1
      if(nnos.eq.1) return
      do 400 n = 2, nnos
         jdiag(n) = jdiag(n) + jdiag(n-1) + 1
  400 continue
      nad = jdiag(nnos)
      return
      end