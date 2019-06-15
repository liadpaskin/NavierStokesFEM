	  SUBROUTINE MATSET (PROP,NUMMAT,NCPROP,IIN1,IOUT,IERROR,etype,isolver_type)
      
      use modvar, only: grav
      use modwave, only: CompOnda, Amp, Depht, numOnda, angFreq

      IMPLICIT REAL*8 (A-H,O-Z)

      dimension    PROP(NUMMAT,NCPROP)
      character*14 etype

      WRITE(IOUT,2010) 
	
!     LEITURA DE PROPRIEDADES DO MATERIAL
      
      ! READ (IIN1,1000) 
      ! READ (IIN1,1000) etype
      ! READ (IIN1,1000) 

   WRITE(IOUT,2020) etype
   WRITE(*,2020) etype
   
   if ((etype .eq. 'NavierStokes2D') .or. (etype .eq. 'NavierStokes3D'))then
       READ (IIN1,1000)
       DO I = 1, NUMMAT
           READ(IIN1 ,*)  N, VISCOS, POISS, RHO, (PROP(n,ii),ii=5,10)
           IF (N .GT. NUMMAT .OR. N .LE. 0) THEN
              WRITE (*,300) N
              WRITE (IERROR,300) N
             read (*,*)
             STOP
           ENDIF
           IF (VISCOS .eq. 0.d0) THEN
              WRITE (*,400) N,VISCOS,poiss,RHO
              WRITE (IERROR,400) N,VISCOS,poiss,RHO
              read (*,*)
              !STOP
           ENDIF
           xmi = Viscos
           !if (xmi .lt. 1E-06) xmi = 1.d0
           EYOUNG  = (2.d0*(1.d0+POISS)) * xmi
           xlambda = (EYOUNG*POISS)/((1.d0+POISS)*(1.d0-2.d0*POISS))
           !     ARMAZENA PROPRIEDADES DO MATERIAL
           PROP(N,1) = EYOUNG
           PROP(N,2) = POISS
           PROP(N,3) = 0.d0
           PROP(N,4) = RHO
           PROP(N,5) = xmi
           PROP(N,6) = xlambda 
           prop(N,11) = xlambda + 2.d0*xmi
           prop(N,12) = xlambda
           prop(N,13) = xlambda + 2.d0*xmi
           prop(N,14) = xmi
           WRITE (IOUT,2030) I,EYOUNG,POISS,1.d0,RHO
           !if (xmi .ne. Viscos) then 
           !    xmi = Viscos
           !    PROP(N,5) = xmi
           !endif
       ENDDO
   else
       write (*,*) 'Este solver e especifico para NavierStokes3D. &
       ou NavierStokes2D; Favor checar o etype.'
   endif

  
 1000 FORMAT (A)
 1010 FORMAT (I8,7F10.0)
 2010 FORMAT (/,' DEFINICAO DOS MATERIAIS ' ,/) 
 
 2020 FORMAT (/,' TIPO DE PROBLEMA ' ,(a),/) 
                  
 2030 FORMAT (' MATERIAL NO. ',I5,//,                                    &
              '   MOD. ELASTICIDADE  . . . ( EYOUNG) = ',1PE15.8,/,      &  
              '   POISSON  . . . . . . . . ( POISS ) = ',1PE15.8,/,      &
              '   ESPESSURA    . . . . . . ( BULK  ) = ',1PE15.8,/,      &
              '   MASSA ESPECIFICA . . . . ( RHO   ) = ',1PE15.8,/)

 2040 FORMAT (' MATERIAL NO. ',I5,//,                                    &
              '   MOD. ELASTICIDADE  . . . ( EYOUNG) = ',1PE15.8,/,      &  
              '   POISSON  . . . . . . . . ( POISS ) = ',1PE15.8,/,      &
              '   MASSA ESPECIFICA . . . . ( RHO   ) = ',1PE15.8,/)
  
 2050 FORMAT (' MATERIAL NO. ',I5,//,                                    &
              '   CONTUTIVIDADE TERMICA. . ( KX    ) = ',1PE15.8,/,      &  
              '   CONDUTIVIDADE TERMICA. . ( KY    ) = ',1PE15.8,/,      &
              '   CALOR ESPEC�FICO     . . ( CE    ) = ',1PE15.8,/,      &
              '   MASSA ESPECIFICA . . . . ( RHO   ) = ',1PE15.8,/)

 2060 FORMAT (' MATERIAL NO. ',I5,//,                                    &
              '   CONTUTIVIDADE TERMICA. . ( KX    ) = ',1PE15.8,/,      &  
              '   MASSA ESPECIFICA . . . . ( RHO   ) = ',1PE15.8,/)
  
                
  300 FORMAT (' *** (MATSET.F90) ERRO DE DADOS: GRUPO DE MATERIAL ',         &
      		'INV�LIDO (',I5,')')

  400 FORMAT (' *** (MATSET.F90) ERRO : Prop nula (Eyoung, poiss, thic)',i5,3e12.5 )
      
  500 FORMAT (' *** (MATSET.F90) ERRO : Prop nula (Eyoung, poiss)',i5,2e12.5 )

  600 FORMAT (' *** (MATSET.F90) ERRO : Prop nula (Kxx)',i5,1e12.5 )

RETURN


END 
