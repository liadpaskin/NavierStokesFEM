      subroutine genrcm ( neq, xadj, adjncy, perm, mask, xls , i6, nnos)
!. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
! .                                                                 .
! .   general reverse cuthill mckee                                 .
! .                                                                 .
! .   purpose - genrcm finds the reverse cuthill-mckee ordering     .
! .      for a general graph. for each connected component in       .
! .      the graph, genrcm obtains the ordering by calling the      .
! .      subroutine rcm.                                            .
! .                                                                 .
! .   input parameters -                                            .
! .      neq - number of nodes                                      .
! .      (xadj, adjncy) - array pair containing the adjacency       .
! .             structure of the graph of the matrix.               .
! .                                                                 .
! .   output parameters -                                           .
! .      perm - vector that contains the rcm ordering.              .
! .                                                                 .
! .   working parameters -                                          .
! .      mask - is used to mark variables that have been            .
! .             numbered during the ordering process. it            .
! .             is initialized to 1, and set to zero as each        .
! .             node is numbered.                                   .
! .      xls - the index vector for a level structure. the level    .
! .            structure is stored in the currently unused spaces   .
! .            in the permutation vector perm.                      .
! .                                                                 .
! .   program subroutines - fnroot, rcm.                            .
! .                                                                 .
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
      integer adjncy(i6), mask(0:neq), perm(i6), xls(neq)
      integer xadj(neq+1), ccsize, i, neq, nlvl, num, root

      perm(:) = 0

      do  i = 1, neq
         mask(i) = 1
  100 continue
      enddo
      num = 1
      do i = 1, neq
!        -----------------------------------
!        for each masked connected component
!        -----------------------------------
         if (mask(i) .eq. 0) go to 200
            root = i
!           -------------------------------------------------------
!           first find a pseudo peripheral node root.
!           note that the level structure found by fnroot is
!           stored starting at perm(num). then rcm is called
!           to order the component using root as the starting node.
!         -------------------------------------------------------
      call fnroot ( root, xadj, adjncy, mask, nlvl, xls, perm(num), neq , i6 )
      call rcm (root, xadj, adjncy, mask, perm(num), ccsize, xls , neq , i6 )
      num = num + ccsize
      if (num .gt. neq) return
  200 continue
      enddo

      return
      end
