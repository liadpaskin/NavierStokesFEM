      subroutine local ( global,rlocal,lm,numel,nd,nvec,neq,iel0,nsize,task)
!c
!c.... routine to perform local/global computations
!c
      character*8        task

      real*8   rlocal ( nsize , nd ),global ( 0 : neq    )

      integer            lm     ( numel , nd   )
!c
!c.... 'LOCALIZE'
!c
      if ( task .eq. 'localize' )       then


!c     inclusao de if
          do kd = 1 , nd
            do ie = 1 , nvec
              if (lm (iel0+ie,kd  ).ge.0) then
                rlocal ( ie,kd   ) = global ( lm (iel0+ie,kd  ) )
              endif
		  enddo       !  kd
          enddo         !  ie

          return

      endif
!c
!c.... 'GLOBALIZ'
!c
      if ( task .eq. 'globaliz' )       then

          do kd = 1 , nd
!c:dir ignore recrdeps
!            
!c     inclusao de if
		  
		  do ie = 1 , nvec
              if ( lm (iel0+ie,kd  ).ge.0) then
                global ( lm (iel0+ie,kd  ) ) = rlocal ( ie,kd   )
              endif
		  enddo       !  kd
          enddo       !  ie

          return

      endif
!c
!c.... 'ADD GLOB'
!c
      if ( task .eq. 'add glob' )       then

          do kd = 1 , nd
!c:dir ignore recrdeps
!c     inclusao de if            
		  do ie = 1 , nvec
	        if ( lm(iel0+ie,kd  ).ge.0) then
              global ( lm(iel0+ie,kd  )) = global ( lm(iel0+ie,kd  ) ) + rlocal ( ie , kd   )
	        endif
            enddo       !  kd
          enddo         !  ie

          return

      endif
!c
!c.... end
!c
      end
