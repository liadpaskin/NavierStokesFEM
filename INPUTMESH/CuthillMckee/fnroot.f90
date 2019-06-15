      subroutine fnroot ( root, xadj, adjncy, mask, nlvl, xls, ls  , neq , i6)
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
! .                                                                 .
! .   find pseudo-peripheral node                                   .
! .                                                                 .
! .   purpose - fnroot implements a modified version of the         .
! .      scheme by gibbs, pole, and stockmeyer to find pseudo-      .
! .      peripheral nodes. it determines such a node for the        .
! .      section subgraph specified by mask and root.               .
! .                                                                 .
! .   input parameters -                                            .
! .      (xadj, adjncy) - adjacency structure pair for the graph.   .
! .      mask - specifies a section subgraph. nodes for which       .
! .             mask is zero are ignored by fnroot.                 .
! .                                                                 .
! .   updated parameters -                                          .
! .      root - on input, it (along with mask) defines the          .
! .             component for which a pseudo-peripheral node is     .
! .             to be found. on output, it is the node obtained.    .
! .                                                                 .
! .   output parameters -                                           .
! .      nlvl - is the number of levels in the level structure      .
! .             rooted at the node root.                            .
! .      (xls,ls) - the level structure array pair containing       .
! .                 the level structure found.                      .
! .                                                                 .
! .   program subroutines - rootls.                                 .
! .                                                                 .
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
      integer adjncy(i6), ls(neq), mask(0:neq), xls(neq), xadj(neq+1)
      integer ccsize, j, jstrt, k, kstop, kstrt, mindeg, nabor, ndeg, nlvl, node, nunlvl, root
!
!     ---------------------------------------------
!     determine the level structure rooted at root.
!     ---------------------------------------------
     call rootls ( root, xadj, adjncy, mask, nlvl, xls, ls,neq ,i6)
      ccsize = xls(nlvl+1) - 1
      if (nlvl .eq. 1 .or. nlvl .eq. ccsize) return
!     ----------------------------------------------------
!     pick a node with minimum degree from the last level.
!     ----------------------------------------------------
  100 jstrt = xls(nlvl)
      mindeg = ccsize
      root = ls(jstrt)
      if (ccsize .eq. jstrt) go to 400
         do 300 j = jstrt, ccsize
            node = ls(j)
            ndeg = 0
            kstrt = xadj(node)
            kstop = xadj(node+1) - 1
            do 200 k = kstrt, kstop
               nabor = adjncy(k)
               if (mask(nabor) .gt. 0) ndeg = ndeg + 1
  200       continue
            if ( ndeg .ge. mindeg ) go to 250
               root = node
               mindeg = ndeg
  250       continue
  300    continue
!     ---------------------------------------
!     and generate its rooted level structure.
!     ---------------------------------------
  400 call rootls ( root, xadj, adjncy, mask, nunlvl, xls, ls,neq, i6 )
      if ( nunlvl .le. nlvl ) return
         nlvl = nunlvl
         if ( nlvl .lt. ccsize ) go to 100
      return
      end