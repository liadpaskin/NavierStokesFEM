      subroutine matvec ( p,ap,ve1,ve2,neq,ielblk,nsize,nd,nblock,numel,s,lm )
!c
!c.... routine to compute EBE matrix-vector multiplication
!c
      real*8    ap  ( 0 : neq ) ,p   ( 0 : neq ),ve1 ( nsize , nd  ),ve2 ( nsize , nd  ),s   (  nd  , nd , numel )

      integer             lm     ( numel , nd   ),ielblk ( 1 )

!c
!c.... loop over element blocks
!c
      ielm  = 1

      do iblock = 1 , nblock

        nvec = ielblk ( iblock )
        iel0 = ielm - 1
        p(0) = 0.d0
!c
!c.... gather global components of 'p' to 've1'
!c
        call local ( p,ve1,lm,numel,nd,nvec,neq,iel0,nsize,'localize' )
!c
!c.... EBE matrix-vector multiplication for all elements in this block
!c
        ve2(:,:) = 0.d0

        do j = 1 , nd
          do i = 1 , nd
            do ie = 1 , nvec
              ve2(ie,i) = ve2(ie,i) + s(i,j,iel0+ie) * ve1(ie,j)
            enddo
          enddo
        enddo

!c
!c.... scatter local components of 've2' to 'ap'
!c
        ap (:) = 0.d0
        call local (ap,ve2,lm,numel,nd,nvec,neq,iel0,nsize,'add glob' )
        ap (0) = 0.d0

!c
!c.... increment element block counter
!c
        ielm = ielm + nvec
!c
!c.... end loop over the blocks
!c
      enddo                    ! iblock
!c
!c.... return
!c
      return
      end
