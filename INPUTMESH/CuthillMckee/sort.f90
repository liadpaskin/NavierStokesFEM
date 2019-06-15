      subroutine sort (ip,ips,noviz,nnos,i3,i4,i5,i6)
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
! .                                                                 .
! .                                                                 .
! .                                                                 .
! .                                                                 .
! .                                                                 .
! .                                                                 .
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
      dimension ip(nnos+1), ips(i6), noviz(i6)


      do 200 i = 1, nnos
         l = ips(i)
         ipasso = 1
         itroca = 1
   50    continue
         if (ipasso .lt. l .and. itroca .eq. 1) then
            itroca = 0
            ind = ip(i)
            do 100 j = 1, l-ipasso
               if (noviz(ind+1) .lt. noviz(ind)) then
                  iaux = noviz(ind)
                  noviz(ind) = noviz(ind+1)
                  noviz(ind+1) = iaux
                  itroca = 1
               endif
               ind = ind + 1
  100       continue
            ipasso = ipasso + 1
            go to 50
         endif
  200 continue
      do 400 i = 1, nnos
         l = ip(i)
         k = ip(i+1)-2
         do 300 j = l, k
            if (noviz(j) .ne. 0)    no = noviz(j)
            if (noviz(j+1) .eq. no) then
               noviz(j+1) = 0
               ips(i) = ips(i) - 1
            endif
  300    continue
  400 continue
      do 500 i = 1, nnos
         ip(i+1) = ip(i) + ips(i)
  500 continue
      j = 1
      do 600 i = 1, i5-i4
         if (noviz(i) .ne. 0) then
            ips(j) = noviz(i)
            j = j+1
         endif
  600 continue
      i4 = i3 + ip(nnos+1) -1
      return
      end