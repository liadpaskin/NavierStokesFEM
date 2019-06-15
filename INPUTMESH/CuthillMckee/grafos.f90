      subroutine grafos (ndn,ip,ips,noviz,nelm,nnomax,nnos,i3,i6,nume1,nnoelm)
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
! .                                                                 .
! .                                                                 .
! .                                                                 .
! .                                                                 .
! .                                                                 .
! .                                                                 .
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
      integer     ::   nnos, ndn(nnomax,nelm), ip(nnos+1), ips(i6), noviz(i6)


      do 400 i = nume1, nelm
         !itype = 4
         !if (itype.eq.1)                 nnoelm = 2
         !if (itype.eq.2 .or. itype.eq.3) nnoelm = 3
         !if (itype.eq.4)                 nnoelm = 4
         !if (itype.eq.5)                 nnoelm = 6
         !if (itype.eq.6 .or. itype.eq.7) nnoelm = 8
         !if (itype.eq.8)                 nnoelm = 20
         do 300 j = 1, nnoelm
            noj = ndn(j,i)
            do 200 k = 1, nnoelm
               nok = ndn(k,i)
               if (nok .eq. noj) go to 100
               ipos = ip(noj) + ips(noj)
               ips(noj) = ips(noj) + 1
               noviz(ipos) = nok
  100          continue
  200       continue
  300    continue
  400 continue
      return
      end