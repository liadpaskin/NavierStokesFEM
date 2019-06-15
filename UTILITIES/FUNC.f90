!    ---------------
      FUNCTION FUNC  ( T )
!    ---------------

      use    modtim ,    only: ifunc,nptf,TIMEF,a_FUNC,b_FUNC,omega_FUNC,phase_FUNC, TFTIM , TFVAL, KPFOUN

      IMPLICIT REAL*8 (A-H,O-Z)

!-----------------------------------------------------------------------
!
!     FUNCAO DE CARGA NO TEMPO
!
!            T     --> INSTANTE DE TEMPO EM QUESTAO
!            IFUNC --> TIPO DA FUNCAO
!
!                      = 1     HEAVYSIDE ( 1.0D0 , TEND )
!                      = 2     SENOIDAL :  A + B SIN ( WT + PHI )
!                      = 3     TIME-HISTORY
!
!            TEND  --> TEMPO FINAL DA ATUACAO DO CARREGAMENTO
!
!-----------------------------------------------------------------------

     ! COMMON / TIME /    ifunc,nptf,tend,a,b,omega,phase
     ! COMMON / TIMEV /   TFTIM(200) , TFVAL(200) , KPFOUN


      FUNC = 0.D0
      
      KPFOUN = 1

      IF ( T .LE. TIMEF ) THEN
         IF ( IFUNC .EQ. 1 ) THEN
              FUNC = 1.0D0
         
         ELSEIF ( IFUNC .EQ. 2 ) THEN
              FUNC = A_FUNC + B_FUNC * SIN ( OMEGA_FUNC*T + PHASE_FUNC )
         
         ELSEIF ( IFUNC .EQ. 3 ) THEN
         
              DO 100 I = KPFOUN , NPTF-1
                 IF ( TFTIM(I).LE.T .AND. T.LE.TFTIM(I+1) ) THEN
                      FUNC = TFVAL(I) + ( TFVAL(I+1) - TFVAL(I) ) / ( TFTIM(I+1) - TFTIM(I) ) * ( T - TFTIM(I))
                      KPFOUN = I
                      RETURN
                 ENDIF
  100         CONTINUE
         ENDIF
      ENDIF


      RETURN
      END
