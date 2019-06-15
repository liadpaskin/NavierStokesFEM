      subroutine rcm ( root, xadj, adjncy, mask, perm, ccsize, deg, neq , i6 )
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
! .                                                                 .
! .   reverse cuthill-mckee ordering                                .
! .                                                                 .
! .   purpose - rcm numbers a connected component specified         .
! .      by mask and root, using the rcm algorithm.                 .
! .      the numbering is to be started at the node root.           .
! .                                                                 .
! .   input parameters -                                            .
! .      (xadj, adjncy) - adjacency structure pair for the graph.   .
! .      root - is the node that defines the connected component    .
! .             and it is used as the starting node for the rcm     .
! .             ordering.                                           .
! .                                                                 .
! .   updated parameters -                                          .
! .      mask - only those nodes with nonzero input mask values     .
! .             are considered by the routine. the nodes numbered   .
! .             by rcm will have their mask values set to zero.     .
! .                                                                 .
! .   output parameters -                                           .
! .      perm - will contain the rcm ordering. level structure      .
! .      ccsize - is the size of the connected component that       .
! .               has been numbered by rcm.                         .
! .                                                                 .
! .   working parameters -                                          .
! .      deg - is a temporary vector used to hold the degree of     .
! .            the nodes in the section graph specified by mask     .
! .            and root.                                            .
! .                                                                 .
! .   program subroutines - degree.                                 .
! .                                                                 .
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
      integer adjncy(i6), deg(neq), mask(0:neq), perm(neq), xadj(neq+1)
      integer ccsize, fnbr, i, j, jstop, jstrt, k, l, lbegin,lnbr, lperm, lvlend, nbr, node,  root

!     ----------------------------------------------
!     find the degrees of the nodes in the component
!     specified by mask and root.
!     ----------------------------------------------
      call degree ( root, xadj, adjncy, mask, deg, ccsize, perm , neq, i6 )
      mask(root) = 0
      if ( ccsize .le. 1 ) return
      lvlend = 0
      lnbr = 1
!     ----------------------------------------------
!     lbegin and lvlend point to the beginning and
!     the end of the current level respectively.
!     ----------------------------------------------
  100 lbegin = lvlend + 1
      lvlend = lnbr
      do 600 i = lbegin, lvlend
!        ------------------------------
!        for each node in current level
!        ------------------------------
         node = perm(i)
         jstrt = xadj(node)
         jstop = xadj(node+1) - 1
!        -----------------------------------------
!        find the unnumbered neighbors of node.
!        fnbr and lnbr point to the first and last
!        unnumbered neighbors respectively of the
!        current node in perm.
!        -----------------------------------------
         fnbr = lnbr + 1
         do 200 j = jstrt, jstop
            nbr = adjncy(j)
            if ( mask(nbr) .eq. 0 ) go to 150
               lnbr = lnbr + 1
               mask(nbr) = 0
               perm(lnbr) = nbr
  150       continue
  200    continue
         if ( fnbr .ge. lnbr ) go to 550
!           -------------------------------------------
!           sort the neighbors of node in increasing
!           order by degree. linear insertion is used.
!           -------------------------------------------
            k = fnbr
  300       l = k
               k = k + 1
               nbr = perm(k)
  400          if ( l .lt. fnbr ) go to 500
                  lperm = perm(l)
                  if ( deg(lperm) .le. deg(nbr) ) go to 500
                     perm(l+1) = lperm
                     l = l - 1
                     go to 400
  500          perm(l+1) = nbr
               if ( k .lt. lnbr ) go to 300
  550    continue
  600 continue
      if ( lnbr .gt. lvlend ) go to 100
!     ---------------------------------------
!     we now have the cuthill mckee ordering.
!     reverse it below ...
!     ---------------------------------------
      k = ccsize/2
      l = ccsize
      do 700 i = 1, k
         lperm = perm(l)
         perm(l) = perm(i)
         perm(i) = lperm
         l = l - 1
  700 continue
      return
      end
