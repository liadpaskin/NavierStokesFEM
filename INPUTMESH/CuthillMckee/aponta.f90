      subroutine aponta (ndn,ip,nnos,nelm,nnomax,i5,nnoelm)
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
! .                                                                 .
! .                                                                 .
! .                                                                 .
! .                                                                 .
! .                                                                 .
! .                                                                 .
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
      dimension ndn(nnomax,nelm), ip(nnos+1)

      ip(1) = 1
      do 200 i = 1, nelm
         !itype = 4
         !if (itype.eq.1)                 nnoelm = 2
         !if (itype.eq.2 .or. itype.eq.3) nnoelm = 3
         !if (itype.eq.4)                 nnoelm = 4
         !if (itype.eq.5)                 nnoelm = 6
         !if (itype.eq.6 .or. itype.eq.7) nnoelm = 8
         !if (itype.eq.8)                 nnoelm = 20
         do 100 j = 1, nnoelm
            no = ndn(j,i)
            ip(no+1) = ip(no+1) + nnoelm - 1
  100    continue
  200 continue
      do 300 i = 2, nnos+1
         ip(i) = ip(i) + ip(i-1)
  300 continue
      i5 = ip(nnos+1) - 1
      return
      end