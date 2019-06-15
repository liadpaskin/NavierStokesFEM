!    -------------------------------------------------------------------
      SUBROUTINE CONTRL
!    -------------------------------------------------------------------
!
!     CHAMADA POR : PROGRAMA PRINCIPAL
!     RETORNA     : INFORMACOES DE CONTROLE E DADOS GERAIS
!
      use    modtapes, only : iin1 , iin2, iout
      use    modvar,   only : NUMNP,NUME,NUMMAT,NNOEL,NGL, ND, HED, NOPLOT, NGLPLOT, NDAUX ,&
       isolver_type , ireordering , kmax , lmax , etol , grav , niter, kmax1 , nelblk, nsize, etype, nbody,elrefp
      use    modtim,   only : TIMEF , DT,  NIMP,NFLAG, DAMP1, DAMP2 , ALPHA , BETA , GAMMA ,&
       IFUNC , NPTF,  A_FUNC,  B_FUNC,  OMEGA_FUNC ,  PHASE_FUNC, TFTIM , TFVAL,NSTEPS , it1, clockcount1

      IMPLICIT REAL*8 (A-H,O-Z)

      CHARACTER*70 CDATE ,CTIME

      CALL Date_and_Time(CDATE ,CTIME)


      READ (IIN1,1000)
      READ (IIN1,1000) HED
      READ (IIN1,1000)

      READ (IIN1,1000)
      READ (IIN1,*) NUMNP,NUME,NUMMAT,NNOEL,NGL,NOPLOT,NGLPLOT, iSOLVER_TYPE , iREORDERING ,&
       kmax , lmax , etol , niter, nbody, elrefP
      READ (IIN1,1000)


      ND = NNOEL * NGL
      NDAUX = (nd*nd + nd)/2
  
      READ (IIN1,1000)
      READ (IIN1,*) TIMEF , DT,   NIMP, DAMP1, DAMP2 , ALPHA , BETA , GAMMA , IFUNC , NPTF, &
       A_FUNC,  B_FUNC,  OMEGA_FUNC ,  PHASE_FUNC , grav
      READ (IIN1,1000)
      
      write(*,*) huge(nsteps)
      NSTEPS = int(TIMEF/DT)
      NFLAG  = int(NSTEPS/NIMP)
      if (nflag.lt.1) nflag = 1

      
      IF ((IFUNC.EQ.3).AND.(NPTF.GT.0)) THEN
      WRITE (IOUT,*) 
      WRITE (IOUT,*) ' PONTOS DA FUNCAO F(t)'
      READ (IIN1,1000)
  
         DO I = 1, NPTF
         
         READ(IIN1,*) TFTIM(I) , TFVAL(I)
         WRITE(IOUT,*) TFTIM(I) , TFVAL(I)
         ENDDO 
         
         READ (IIN1,1000)

      ELSE
      WRITE (IOUT,*) ' NPTF = 0, ADOTANDO IFUNC = 1 '
      WRITE (IOUT,*)
      IFUNC = 1
         
      ENDIF
      
      CALL SYSTEM_CLOCK(it1, count_rate = clockcount1) 
      
      WRITE (IOUT,2100) HED,CDATE,CTIME
      WRITE (IOUT,2200) NUMNP,NUME,NUMMAT,NNOEL,NGL

      !WRITE (IOUT,2300) TIMEF , DT,  NIMP, DAMP1, DAMP2 , ALPHA , BETA , GAMMA , IFUNC , NPTF
      
      !if (isolver_type .eq. 3) then
          kmax1   = kmax + 1
          nsize   = nume
          nelblk  = 1
      !endif    
      
      READ (IIN1,1000) 
      READ (IIN1,1000) etype
      READ (IIN1,1000) 
      
      RETURN
    
 1000 FORMAT (A)
 
 2100 FORMAT (' TITULO  : ',A,// &
             ,' DATA    : ',A,/  &
             ,' HORA    : ',A,//)

 2200 FORMAT ( &
      ' INFORMACOES DE CONTROLE  -  DADOS GERAIS '            ,///,  &
      ' NUMERO DE PONTOS NODAIS . . . . . . . . . (NUMNP)  = ',I8//, &
      ' NUMERO DE ELEMENTOS . . . . . . . . . . . (NUME)   = ',I8//, &
      ' NUMERO DE TIPOS DE MATERIAIS  . . . . . . (NUMMAT) = ',I8//, &
      ' NUMERO DE NOS POR ELEMENTO  . . . . . . . (NNOEL)  = ',I8//, &      
      ' NUMERO DE GRAUS DE LIBERDADE POR NO . . . (NGL)    = ',I8/)  

 2300 FORMAT ( &
      ' INFORMACOES DE CONTROLE  -  DADOS AN�LISE '            ,///,  &
      ' TEMPO DO TRANSIENTE     . . . . . . . . . (TIMEF)  = ',E12.5//, &
      ' PASSO DE TEMPO      . . . . . . . . . . . (DT)     = ',E12.5//, &
      ' NUMERO DE IMPRESSOES  DE TEMPO  . . . . . (NIMP)   = ',E12.5//, &
      ' ALPHA DO AMORTECIMENTO        . . . . . . (DAMP1)  = ',E12.5//, &
      ' BETA DO AMORTECIMENTO   . . . . . . . . . (DAMP2)  = ',E12.5//, &      
      ' ALPHA DO ALGORITMO 1a ORDEM   . . . . . . (ALPHA)  = ',E12.5//, &
      ' BETA  DO ALGORITMO 2a ORDEM   . . . . . . (BETA)   = ',E12.5//, &
      ' GAMMA DO ALGORITMO 2a ORDEM   . . . . . . (GAMMA)  = ',E12.5//, &
      ' TIPO DA F(t)    . . . . . . . . . . . . . (IFUNC)  = ',I8//, &
      ' NUMERO DE PONTOS DA F(t)    . . . . . . . (NPTF)   = ',I8/)  
      END


