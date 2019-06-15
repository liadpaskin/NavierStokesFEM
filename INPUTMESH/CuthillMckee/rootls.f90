     subroutine rootls ( root, xadj, adjncy, mask, nlvl, xls, ls , neq,i6)
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
! .                                                                 .
! .   rooted level structure                                        .
! .                                                                 .
! .   purpose - rootls generates the level structure rooted         .
! .      at the input node called root. only those nodes for        .
! .      which mask is nonzero wil be considered.                   .
! .                                                                 .
! .   input parameters -                                            .
! .      root - the node at which the level structure is to         .
! .             be rooted.                                          .
! .      (xadj, adjncy) - adjacency structure pair for the          .
! .             given graph.                                        .
! .      mask - specifies a section subgraph. nodes for which       .
! .             mask is zero are ignored.                           .
! .                                                                 .
! .   output parameters -                                           .
! .      nlvl - is the number of levels in the level structure.     .
! .      (xls,ls) - array pair for the rooted level structure.      .
! .                                                                 .
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
      integer adjncy(i6), ls(neq), mask(0:neq), xls(neq), xadj(neq+1)
      integer i, j, jstop, jstrt, lbegin, ccsize, lvlend, lvsize, nbr, nlvl, node, root


!     -------------------
!     initialization ...
!     -------------------
      mask(root) = 0
      ls(1) = root
      nlvl = 0
      lvlend = 0
      ccsize = 1
!     ------------------------------------------------------
!     lbegin is the pointer to the beginning of the current
!     level, and lvlend points to the end of this level.
!     ------------------------------------------------------
  200 lbegin = lvlend + 1
      lvlend = ccsize
      nlvl = nlvl + 1
      xls(nlvl) = lbegin
!     ------------------------------------------------------
!     generate the next level by finding all the masked
!     neighbors of nodes in the current level.
!     ------------------------------------------------------
      do 400 i = lbegin, lvlend
         node = ls(i)
         jstrt = xadj(node)
         jstop = xadj(node+1) - 1
         if ( jstop .lt. jstrt ) go to 350
            do 300 j = jstrt, jstop
               nbr = adjncy(j)
               if ( mask(nbr) .eq. 0 ) go to 250
                  ccsize = ccsize + 1
                  ls(ccsize) = nbr
                  mask(nbr) = 0
  250          continue
  300       continue
  350    continue
  400 continue
!     ------------------------------------------------------
!     compute the current level width.
!     if it is nonzero, generate the next level.
!     ------------------------------------------------------
      lvsize = ccsize - lvlend
      if ( lvsize .gt. 0 ) go to 200
!     --------------------------------------------------------
!     reset mask to one for the nodes in the level structure.
!     --------------------------------------------------------
      xls(nlvl+1) = lvlend + 1
      do 500 i = 1, ccsize
         node = ls(i)
         mask(node) = 1
  500 continue
      return
      end