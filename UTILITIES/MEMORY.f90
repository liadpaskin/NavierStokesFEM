     SUBROUTINE MEMORY (N,IFLAG,XVALUE)
     ! AVALIA EM MEGABYTES O TAMANHO DE UM ARRANJO DE N POSICOES
     ! iflag = 1 (inteiro)
     ! iflag = 1 (real*8)
     IMPLICIT NONE
     
     REAL*8 XVALUE,xmega
     INTEGER,INTENT(IN) :: N, IFLAG
     
     xMEGA = 1.0e6
     
     IF (IFLAG.EQ.1) THEN
      
     XVALUE = N / xMEGA
     
     ELSEIF (IFLAG.EQ.2) THEN
     
     XVALUE = N * 8.d0 / xMEGA
     
     ENDIF
     
     RETURN
     END 