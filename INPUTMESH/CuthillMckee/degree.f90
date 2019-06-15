     subroutine degree ( root, xadj, adjncy, mask, deg, ccsize, ls ,neq,i6)
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
! .                                                                 .
! .   degree in masked component                                    .
! .                                                                 .
! .   purpose - degree computes the degrees of the nodes in         .
! .      the connected component specified by mask and root.        .
! .      nodes for which mask is zero are ignored.                  .
! .                                                                 .
! .   input parameters -                                            .
! .      root - is the input node that defines the component.       .
! .      (xadj, adjncy) - adjacency structure pair.                 .
! .      mask - specifies a section subgraph.                       .
! .                                                                 .
! .   output parameters -                                           .
! .      deg - array containing the degrees of the nodes in         .
! .            the component.                                       .
! .      ccsize - size of the component specified by mask and root. .
! .                                                                 .
! .   working parameters -                                          .
! .      ls - a temporary vector used to store the nodes of the     .
! .           component level by level.                             .
! .                                                                 .
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
   
      integer adjncy(i6), deg(neq), ls(neq), mask(0:neq), xadj(neq+1)
      integer i, j, jstop, jstrt, lbegin, ccsize, lvlend, lvsize, nbr, ideg, node, root

!     --------------------------------------------------
!     initialization ...
!     the array xadj is used as a temporary marker to
!     indicate which nodes have been considered so far.
!     --------------------------------------------------
      ls(1) = root
      xadj(root) = -xadj(root)
      lvlend = 0
      ccsize = 1
!     ------------------------------------------------------
!     lbegin is the pointer to the beginning of the current
!     level, and lvlend points to the end of this level.
!     ------------------------------------------------------
  100 lbegin = lvlend + 1
      lvlend = ccsize
!     ------------------------------------------------------
!     find the degrees of nodes in the current level,
!     and at the same time, generate the next level.
!     ------------------------------------------------------
      do 400 i = lbegin, lvlend
         node = ls(i)
         jstrt = -xadj(node)
         jstop =  iabs(xadj(node+1)) - 1
         ideg = 0
         if ( jstop .lt. jstrt ) go to 300
            do 200 j = jstrt, jstop
               nbr = adjncy(j)
               if ( mask(nbr) .eq. 0 ) go to 150
                  ideg = ideg + 1
                  if ( xadj(nbr) .lt. 0 ) go to 150
                     xadj(nbr) = -xadj(nbr)
                     ccsize = ccsize + 1
                     ls(ccsize) = nbr
  150          continue
  200       continue
  300    deg(node) = ideg
  400 continue
!     ------------------------------------------------------
!     compute the current level width.
!     if it is nonzero, generate the next level.
!     ------------------------------------------------------
      lvsize = ccsize - lvlend
      if ( lvsize .gt. 0 ) go to 100
!     --------------------------------------------------------
!     reset xadj to its correct sign and return.
!     --------------------------------------------------------
      do 500 i = 1, ccsize
         node = ls(i)
         xadj(node) = -xadj(node)
  500 continue
      return
      end