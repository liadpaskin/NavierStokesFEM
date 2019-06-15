 subroutine getNALHS(ngl,nnoel,nume,neq,lm,nALHS)
        
    implicit none
    
    integer, intent(in) :: ngl,nnoel,nume,neq, lm(nume,ngl*nnoel)
    integer*8, intent(out):: nALHS
    
    integer*8 :: iel, jnd, innoel, ingl
    integer*8, allocatable :: lm_parab(:,:,:),lm_aux(:,:), idiag(:)
    
    allocate (idiag(neq), lm_aux(ngl*nnoel,nume), lm_parab(ngl,nnoel,nume))
    Lm_aux  (:,  :) = 0
    Lm_Parab(:,:,:) = 0
    idiag(:)        = 0
    
    do iel=1,nume                         !
       do jnd=1, ngl*nnoel                      !
           lm_aux(jnd,iel) = lm(iel,jnd)  ! Transforma lm(nume,nd) -> lm(nd,nume)
       enddo                              !                         = lm (ned,nen,numel)
    enddo 
    do iel=1,nume                                                           !
       do innoel=1, nnoel                                                   !
           do ingl = 1, ngl                                                 !
                lm_parab(ingl,innoel,iel) = lm_aux((innoel-1)*ngl+ingl,iel) ! Transforma lm(nd,nume) -> lm (ned,nen,numel)
           enddo                                                            !       
       enddo                                                                !
    enddo                                                                   !
    call colht(idiag,lm_Parab,ngl,nnoel,nume,neq) ! Forma o idiag
    call diag(idiag,neq,nALHS)  
    
return
end
