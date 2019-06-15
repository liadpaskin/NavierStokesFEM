subroutine matrixTrans2 (matrix1,i1,j1,matrix2)
    
    ! Retorna atraves de matrix2 a transposta da matriz contida no parametro matrix1.
    ! i1 e j1 sao, respectivamente, n. de linhas e de colunas da matriz original.
    
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: i1, j1  
    INTEGER, INTENT(IN)  :: matrix1(i1,j1)
    INTEGER, SAVE       :: i,j
    INTEGER              :: matrix2(j1,i1)
    
    matrix2 = 0.d0
    
    do i = 1,i1
        do j = 1,j1
            matrix2(j,i) = matrix1(i,j)
        enddo
    enddo

return
endsubroutine