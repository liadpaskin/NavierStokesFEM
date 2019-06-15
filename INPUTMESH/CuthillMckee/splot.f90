      subroutine splot (ine , iperm , invp ,nnomax, nnos, nelm , i6,nume1,nnoelm)

!c . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
!c .                                                                 .
!c .                                                    20/04/88     .
!c .   - descricao :                                                 .
!c .                                                                 .
!c .      grava o arquivo para plotter.                              .
!c .                                                                 .
!c .                                                                 .
!c .   - chamada por : contr                                         .
!c .                                                                 .
!c .                                                                 .
!c .   - parametros de entrada :                                     .
!c .                                                                 .
!c .      ine    =  incidencias dos elementos                        .
!c .      x      =  coordenadas dos nos                              .
!c .      ncor   =  dimensao                                         .
!c .      nnoelm =  numero de nos por elemento                       .
!c .      nnos   =  numero de nos                                    .
!c .      nelm   =  numero de elementos                              .
!c .      noini  =  numero do no inicial                             .
!c .      nelini =  numero do elemento inicial                       .
!c .                                                                 .
!c .                                                                 .
!c . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
       character*80 a, b, tit
       integer ine(nnomax,nelm),iperm(i6), invp(nnos)


      do 100 i = 1, nnos
       j = iperm(i)
       invp(j) = i
  100 continue

      do 300 i = nume1, nelm
         ! itype = 4
         ! if (itype.eq.1)                 nnoelm = 2
         ! if (itype.eq.2 .or. itype.eq.3) nnoelm = 3
         ! if (itype.eq.4)                 nnoelm = 4
         ! if (itype.eq.5)                 nnoelm = 6
         ! if (itype.eq.6 .or. itype.eq.7) nnoelm = 8
         ! if (itype.eq.8)                 nnoelm = 20
         do 200 j = 1, nnoelm
            k = ine(j,i)
            ine(j,i) = invp(k)
  200    continue
  300 continue


       !do 400 i = nume1, nelm
       !   !itype = 4
       !   !b = '(10i8)'
       !   !if (itype .eq. 1) then
       !   !   nnoelm = 2
       !   !elseif (itype .eq. 2) then
       !   !   nnoelm = 3
       !   !elseif (itype .eq. 3) then
       !   !   nnoelm = 3
       !   !elseif (itype .eq. 4) then
       !   !   nnoelm = 4
       !   !elseif (itype .eq. 5) then
       !   !   nnoelm = 6
       !   !elseif (itype .eq. 6) then
       !   !   nnoelm = 8
       !   !elseif (itype .eq. 7) then
       !   !   nnoelm = 8
       !   !elseif (itype .eq. 8) then
       !   !   nnoelm = 20
       !   !   b = '(2i5,8i5,/,5x,12i5)'
       !   !endif
  !400 !continue

      return
      end