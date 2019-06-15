 subroutine  ensg_tensor2_elem (tipo,tens11,tens12,tens22,nume,etype,istep)
	 use modtapes,     only: filer, iblank, iplt

     implicit none
     
     character*14 	etype
     
     real*8 ::   tens11(nume), tens22(nume), tens12(nume)
     
     character*70 filename
     character*1  tipo
     integer :: iel, istep, nume

     
     write (filename(1:4),'(i4.4)') istep

     write(filename,'(a)') './OUT2/'//filer(1:iblank-1)//'.'//tipo//filename(1:4)
 
     open (unit=iplt,file=filename,form='formatted')
     
     write(*,'(a)') ' Gravando '//filename(1:50)
     
     write (iplt,'(a)') 'Per_elem symmetric tensor values'
     
     write (iplt,'(a/a)') 'part 1','tria6' 
     
     do iel = 1, nume
        write (iplt,'(1p,6e12.5)') tens11(iel), tens22(iel), 0.d0, tens12(iel), 0.d0, 0.d0!(esca(1,i),i=1,nume)
     enddo
     
     close (unit=iplt)
     
     return
end
    
    
    
    
    
    
    
    
    
    
    subroutine  ensg_scalar_elem (tipo,esca,nume,etype,istep)
	 use modtapes,     only: filer, iblank, iplt


          
     implicit real*8 (a-h,o-z)
     
     character*14 	etype

     
     dimension   esca     (1,nume) 
     
     character*70 filename
     character*3  tipo

     
     write (filename(1:4),'(i4.4)') istep

     write(filename,'(a)') './OUT2/'//filer(1:iblank-1)//'.'//tipo//filename(1:4)
 
     open (unit=iplt,file=filename,form='formatted')
     
     write(*,'(a)') ' Gravando '//filename(1:50)
     


     write (iplt,'(a,i5,1x,a)') 'Ensight Escalar passo ',istep

      if     ( (etype.eq.'ept  ') .or. (etype.eq.'epd ').or. (etype.eq.'epdmohrc '))    then     ! *************** ELASTICIDADE 2D ***************

        write (iplt,'(a/a)') 'part 1','tria3'

      elseif (  etype.eq.'calor2d       ' .or. etype.eq.'convec2D       ')                                      then     ! ***************  CALOR 2D  ******************** 
    
        write (iplt,'(a/a)') 'part 1','tria3' 
 
      elseif (  etype.eq.'elasticidade3d')                                      then     ! *************** ELASTICIDADE 3D *************** 
     
        write (iplt,'(a/a)') 'part 1','tetra4' 
        
      elseif (  etype.eq.'truss3d')                                      then     ! *************** TRELI�A 3D *************** 
     
        write (iplt,'(a/a)') 'part 1','bar2'  
        
      elseif (  (etype.eq.'NavierStokes2D') .or. (etype.eq.'GravNavStok2D')  .or. (etype .eq. 'NSWaveMake2D')&
       .or. (etype .eq. 'Stokes2D'))                                      then     ! *************** NavierStokes2D ***************
     
        write (iplt,'(a/a)') 'part 1','tria6'  
      endif





     write (iplt,'(1p,6e12.5)') (esca(1,i),i=1,nume)
     
     close (unit=iplt)
     
     return
     end

! *** output for nodal vetor  

      subroutine ensg_vec (tipo,disp,up,id,ngl,numnp,neq,neqp,istep)
	  
	  use modtapes,     only: filer, iblank, iplt
      use modmesh,      only: invp
      use modvar ,      only: ireordering
      use modtim ,      only: TFUNC
      
      implicit real*8 (a-h,o-z)

      character*70 filename
      character*1  tipo


      dimension   disp (0:neq), id   (ngl,numnp), up (0:neqp) 
      real*8, allocatable ::  ap(:)


      filename = filer
      
!     arquivo Ensight para vetores
      write(*,*)'stop',istep
      write (filename(1:4),'(i4.4)') istep

      write(filename,'(a)') './OUT2/'//filer(1:iblank-1)//'.'//tipo//filename(1:4)

      write(*,'(a)') ' Gravando '//filename(1:50)

      open (unit=iplt,file=filename,form='formatted')
      
      write (iplt,'(a,i5,1x,a)') 'Ensight Vector passo ',istep


! ***  

   
if (neqp.gt.0) then
      
    allocate (ap(-neqp:neq)) ; ap = 0.d0
      
    do i = 1, numnp 
       do j = 1, ngl
       ieq = id(j,i)
	   if (ieq.ge.0) then
	     ap(ieq)=disp(ieq)
       else
         ap(ieq)=up(-ieq)*TFUNC
	   endif
     enddo
    enddo
         
    IF (ireordering .eq. 0) then
              IF (ngl .eq. 1) THEN
                write (iplt,'(1p,6e12.4E3)') (ap(id(1,i)), &
                                              0.d0, &
                                              0.d0, &
        	     			                  i=1,numnp)
              ELSEIF (ngl .eq. 2) THEN
                write (iplt,'(1p,6e12.4E3)') (ap(id(1,i)), &
                                              ap(id(2,i)), &
                                              0.d0, &
        	     			                  i=1,numnp)
                write (iplt,'(1p,6e12.4E3)') (ap(id(1,i)), &
                                              ap(id(2,i)), &
                                              ap(id(3,i)), &
        	     			                  i=1,numnp)
                write (iplt,'(1p,6e12.4E3)') (ap(id(1,i)), &
                                              ap(id(2,i)), &
                                              ap(id(3,i)), &
        	     			                  i=1,numnp)
              ELSE
                write (*,*) 'Error in: Subroutine ensg_vec'
                read (*,*)
                STOP
              ENDIF
    else
              IF (ngl .eq. 1) THEN
                write (iplt,'(1p,6e12.4E3)') (ap(id(1,invp(i))), &
                                              0.d0, &
                                              0.d0, &
        	     			                  i=1,numnp)
              ELSEIF (ngl .eq. 2) THEN
                write (iplt,'(1p,6e12.4E3)') (ap(id(1,invp(i))), &
                                              ap(id(2,invp(i))), &
                                              0.d0, &
        	     			                  i=1,numnp)
              ELSEIF (ngl .eq. 3) THEN
                write (iplt,'(1p,6e12.4E3)') (ap(id(1,invp(i))), &
                                              ap(id(2,invp(i))), &
                                              ap(id(3,invp(i))), &
        	     			                  i=1,numnp)
              ELSEIF (ngl .eq. 6) THEN
                write (iplt,'(1p,6e12.4E3)') (ap(id(1,invp(i))), &
                                              ap(id(2,invp(i))), &
                                              ap(id(3,invp(i))), &
        	     			                  i=1,numnp)
              ELSE
                write (*,*) 'Error in: Subroutine ensg_vec'
                read (*,*)
                STOP
              ENDIF
    endif
                 
          
    deallocate (ap)
      
else
         
    IF (ireordering .eq. 0) then
          IF (ngl .eq. 2) THEN
            write (iplt,'(1p,6e12.4E3)') (disp(id(1,i)), &
                                          disp(id(2,i)), &
                                          0.d0, &
         				                  i=1,numnp)
          ELSEIF (ngl .eq. 3) THEN
            write (iplt,'(1p,6e12.4E3)') (disp(id(1,i)), &
                                          disp(id(2,i)), &
                                          disp(id(3,i)), &
         				                  i=1,numnp)
          ELSEIF (ngl .eq. 6) THEN
            write (iplt,'(1p,6e12.4E3)') (disp(id(1,i)), &
                                          disp(id(2,i)), &
                                          disp(id(3,i)), &
         				                  i=1,numnp)
          endif
    ELSE
          IF (ngl .eq. 2) THEN
            write (iplt,'(1p,6e12.4E3)') (disp(id(1,invp(i))), &
                                          disp(id(2,invp(i))), &
                                          0.d0, &
         				                  i=1,numnp)
          ELSEIF (ngl .eq. 3) THEN
            write (iplt,'(1p,6e12.4E3)') (disp(id(1,invp(i))), &
                                          disp(id(2,invp(i))), &
                                          disp(id(3,invp(i))), &
         				                  i=1,numnp)
          ELSEIF (ngl .eq. 6) THEN
            write (iplt,'(1p,6e12.4E3)') (disp(id(1,invp(i))), &
                                          disp(id(2,invp(i))), &
                                          disp(id(3,invp(i))), &
         				                  i=1,numnp)
          endif
     ENDIF
              
ENDIF
      
   
close (unit=iplt)
   
return
end






!
!      subroutine ensg_sca3dp(tipo,disp0,disp,id,istep)
!      use mvar
!      use mtapes
!      implicit real*8 (a-h,o-z)
!      
!      
!
!      dimension   disp0 (0:neq)	  , &
!                  disp  (0:neqp)  , &
!     		    id   (ngl,numnp),ap(-neqp:neq) 
!      
!      character*70 filename
!      character*1  tipo
!
!!c     arquivo Ensight para escalar 3D
!
!      write (filename(1:3),'(i3.3)') istep
!
!      write(filename,'(a)') filer(1:iblank-1)//'.'//tipo//filename(1:3)
!
!      write(*,'(a)') ' Gravando '//filename(1:50)
!
!      
!      open (unit=iplt,file=filename,form='formatted')
!      
!      write (iplt,'(a,i5,1x,a)') 'Ensight Scalar passo ',istep
!
!      
!      do i = 1, numnp 
!      
!	ieq = id(1,i)
!      
!	if (ieq.ge.0) then
!	ap(ieq)=disp0(ieq)
!      
!	else
!
!      ap(ieq)=disp(-ieq)
!
!	endif
!   
!      enddo
!      
!	write (iplt,'(1p,6e12.5)') (ap(id(1,i)),i=1,numnp)
!      
!      close (unit=iplt)
!      
!      return
!      end
!
!
!!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!
!!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
!      subroutine ensg_esca (tipo,esca,istep,nelx)
! 
!      use mtapes     
!      implicit real*8 (a-h,o-z)
!      
!      dimension   esca     (nelx) 
!      
!      character*70 filename
!      character*2  tipo
!
!      write (filename(1:3),'(i3.3)') istep
!
!      write(filename,'(a)') filer(1:iblank-1)//'.'//tipo//filename(1:3)
!
!  
!      open (unit=iplt,file=filename,form='formatted')
!      
!      write(*,'(a)') ' Gravando '//filename(1:50)
!
!      
!      open (unit=iplt,file=filename,form='formatted')
!      
!      write (iplt,'(a,i5,1x,a)') 'Ensight Escalar passo ',istep
!	  write (iplt,'(a/a)') 'part 1','tria3' 
!      write (iplt,'(1p,6e12.5)') (esca(i),i=1,nelx)
!      
!      close (unit=iplt)
!      
!      return
!      end
!      
!      
!!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
!!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!      subroutine ensg_vec2dp(tipo,dispx0,dispy0,apx,apy,id,istep,neq, &
!                            neqp,ngl,numnp)
!	  use mtapes      
!
!      implicit real*8 (a-h,o-z)
!      
!      
!
!      dimension   dispx0 (0:neq)	, dispy0 (0:neq)		, &
!     		    id   (ngl,numnp),apx(-neqp:neq),apy(-neqp:neq)
!      
!      character*70 filename
!      character*1  tipo
!
!!c     arquivo Ensight para escalar 3D
!
!      write (filename(1:3),'(i3.3)') istep
!
!      write(filename,'(a)') filer(1:iblank-1)//'.'//tipo//filename(1:3)
!
!      write(*,'(a)') ' Gravando '//filename(1:50)
!
!      
!      open (unit=iplt,file=filename,form='formatted')
!      
!      write (iplt,'(a,i5,1x,a)') 'Ensight Vector passo ',istep
!
!      
!      do i = 1, numnp 
!      
!	ieq = id(1,i)
!      
!	if (ieq.gt.0) then
!	
!	apx(ieq)=dispx0(ieq)
!
!	apy(ieq)=dispy0(ieq)
!
!      
!	else
!
!      apx(ieq)=0.d0
!
!	apy(ieq)=0.d0
!
!	endif
!   
!      enddo
!
!
!
!      
!	write (iplt,'(1p,6e12.5)') (apx(id(1,i)), &
!                                 apy(id(1,i)), &
!                                 0.d0,i=1,numnp)
!      
!      close (unit=iplt)
!      
!      return
!      end
!
!
!!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!	subroutine saidanodalp(id,vec,vecp,xtime)
!	
!	
!	use mtapes
!    use mvar
!	
!	implicit real*8 (a-h,o-z)
!
!      COMMON / plot /   neplot(1000) , ngeplot(1000),nnhis
!
!	dimension         vec (0:neq)    , &
!                       vecp (0:neqp)   , &
!          	          id   (ngl,numnp)   
!             
!
!      do ino = 1, nnhis
!
!      
!      if(ID(NGEPLOT(ino),NEPLOT(ino)).ge.0) then
!	
!	write (iplot,*) xtime,VEC(ID(NGEPLOT(ino),NEPLOT(ino)))
!
!    
!      else
!      ieq = - ID(NGEPLOT(ino),NEPLOT(ino))
!	write (iplot,*) xtime,VECp(ieq)
!
!	endif
!      enddo
!
!	return
!	end
!
!
!!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
!
!      


!      subroutine ensg_vec2 (tipo,vec,istep,idvec,nvec)
!      
!      use mtapes
!      use mvar
!      implicit real*8 (a-h,o-z)
!      
!
!
!      dimension   vec (idvec,3)	
!        
!      character*70 filename
!      character*1  tipo
!
!!ccc     arquivo Ensight para vetores
!
!      write (filename(1:4),'(i4.4)') istep
!
!      write(filename,'(a)') filer(1:iblank-1)//'.'//tipo//filename(1:4)
!
!      write(*,'(a)') ' Gravando '//filename(1:50)
!
!      
!      open (unit=iplt,file=filename,form='formatted')
!      
!      write (iplt,'(a,i5,1x,a)') 'Ensight Vector passo ',istep
!      write (iplt,'(1p,6e12.5)') (vec(i,1), &
!                                 vec(i,2), &
!                                vec(i,3), &
!     				                i=1,nvec)
!      
!      close (unit=iplt)
!      
!      return
!      end
!      
!   
!!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!
!      subroutine ensg_par (xyz_par,istep,idpar,numpar)
!	  
!	  use mtapes      
!      implicit real*8 (a-h,o-z)
!      
!      
!      dimension   xyz_par       (idpar,3)         
!      
!      character*70 filename
!                       
!!     arquivo Ensight para particulas
!      
!      write (filename(1:4),'(i4.4)') istep
!      
!      write(filename,'(a)') filer(1:iblank-1)//'.p'//filename(1:4)
!      
!      write(*,'(a)') ' Gravando '//filename(1:50)
!                                  
!      open (unit=iplt,file=filename,form='formatted')                 
!                       
!      write (iplt,'(a)') 'Ensight Particles (BETH) COPPE/UFRJ'
!      write (iplt,'(a,/,i8)') 'particle coordinates',numpar
!      
!      write (iplt,'(i8,1p,3e12.5)') (i,(xyz_par(i,j),j=1,3),i=1,numpar)
!      
!        
!      close (unit=iplt)
!
!      return
!      end
!
!
!!CHICAO 14/01/2002
!      SUBROUTINE ENSG_SCA1 (TIPO,SCA,ISTEP,NUMNP, npinballs)
! 	  
! 	  use mtapes     
!      implicit real*8 (a-h,o-z)
!      
!      
!      dimension   sca     (numnp) 
!      
!      character*70 filename
!      character*1  tipo
!
!!     arquivo Ensight para escalares nodais
!
!      write (filename(1:4),'(i4.4)') istep
!
!      write(filename,'(a)') filer(1:iblank-1)//'.'//tipo//filename(1:4)
!
!      write(*,'(a)') ' Gravando '//filename(1:50)
!
!      
!      open (unit=iplt,file=filename,form='formatted')
!      
!      write (iplt,'(a,i5)') 'Ensight Scalar passo ',istep
!      write (iplt,'(1p,6e12.5)') (sca(i),i=1,npinballs)
!      
!      close (unit=iplt)
!      
!      RETURN
!      END
!!CHICAO
!
!      SUBROUTINE ENSG_ST (tipo,sca,istep,nlin,ncol,ncol1,nres,ilam)
!   	  use mtapes   
!      use mvar 
!      implicit real*8 (a-h,o-z)
!      
!
!      dimension   sca     (nlin,ncol) 
!      
!      character*70 filename
!      character*3  tipo
!
!!     arquivo Ensight para escalares nodais
!
!      write (filename(1:4),'(i4.4)') istep
!
!      write(filename,'(a)') filer(1:iblank-1)//'.'//tipo//filename(1:4)
!
!      write(*,'(a)') ' Gravando '//filename(1:50)
!
!      
!      open (unit=iplt,file=filename,form='formatted')
!      
!      write (iplt,'(a,i5,1x)') 'Ensight Scalar passo ',istep
!     
!      if (ntria.gt.0) then 
!      write (iplt,'(a,i8)') 'part',1
!     
!      write (iplt,'(a)') 'tria3'
!      do k = 1, ilam
!      write (iplt,'(1p,6e12.5)') (sca((i-1)*ilam+k,ncol1),i=1,ntria)
!      enddo
!      endif
!      
!      if (nquad.gt.0) then 
!      write (iplt,'(a,i8)') 'part',2
!      
!      write (iplt,'(a)') 'quad4'
!
!      do k = 1, ilam
!      write (iplt,'(1p,6e12.5)') (sca((i-1)*ilam+k,ncol1),i=ntria+1, &
!                                 nquad+ntria)
!      enddo
!      endif
!      close (unit=iplt)
!      
!      return
!      end
!  


!    --------------------
      SUBROUTINE PLOTTRIA3 ( X,Y,Z,INCID,NUMNP,NUME,NNOEL )
!    --------------------

	  use modtapes      
      IMPLICIT REAL*8 (A-H,O-Z)

 
      DIMENSION      X ( NUMNP  ) ,  &
                     Y ( NUMNP  ) ,  &
                     Z ( NUMNP  ) ,  &
                     INCID ( NUME,NNOEL )


!      CALL BOT ( 'PLOTGE' )


      !WRITE (FILER(IBLANK:IBLANK+4),'(A)') '.TRIA'
       filer = 'tria3.geo'


       OPEN (UNIT=40,FILE=FILER,FORM='FORMATTED')
      
       filer = fileinput     
! ... gravacao no formato ENSIGHT

      WRITE (40,'(A)') 'DADOS PARA FORMATO ENSIGHT5'
      WRITE (40,'(A,/,A,/,A,/,A,/,I8)')                 &
                               'titulo2'          ,    &
                               'node id given'    ,    &
                               'element id given' ,    &
                               'coordinates'      ,    &
                                NUMNP

      WRITE (40,'(I8,3E12.5)') (I,X(I),Y(I),Z(I), I=1, NUMNP )

      WRITE (40,'(A,/,A,/,A,/,I8)')                     &
                                'part 1'           ,    &
                                'malha'            ,    &
                                'tria3'           ,    &
                                 NUME

      WRITE (40,'(4I8)') (J,INCID(J,1),INCID(J,2),INCID(J,3), J=1, NUME ) 
                
      CLOSE (UNIT=40)


!      CALL EOT ( 'PLOTGE' )

      RETURN
      END





   subroutine ensg_vec_plate (tipo,disp,up,id,ngl,numnp,neq,neqp,istep,invp)
	  
	  use modtapes,     only: filer, iblank, iplt
      
      implicit real*8 (a-h,o-z)
      

      dimension   disp (0:neq), id   (ngl,numnp), up (0:neqp) , invp (numnp)
      real*8, allocatable ::  ap(:)
      character*70 filename
      character*1  tipo

!     arquivo Ensight para vetores

      write (filename(1:4),'(i4.4)') istep

      write(filename,'(a)') filer(1:iblank-1)//'.'//tipo//filename(1:4)

      write(*,'(a)') ' Gravando '//filename(1:50)

      
      open (unit=iplt,file=filename,form='formatted')
      
      write (iplt,'(a,i5,1x,a)') 'Ensight Vector passo ',istep


! ***  

   

      
        write (iplt,'(1p,6e12.4E3)') (disp(id(1,invp(i))), &
                                     disp(id(2,invp(i))), &
                                     disp(id(3,invp(i))), &
     				                  i=1,numnp)
   close (unit=iplt)
      
   return
   end



    subroutine ensg_vec_beam (tipo,disp,up,id,ngl,numnp,neq,neqp,istep)
	  
	  use modtapes,     only: filer, iblank, iplt
      
      implicit real*8 (a-h,o-z)
      

      dimension   disp (0:neq), id   (ngl,numnp), up (0:neqp) 
      real*8, allocatable ::  ap(:)
      character*70 filename
      character*1  tipo

!     arquivo Ensight para vetores

      write (filename(1:4),'(i4.4)') istep

      write(filename,'(a)') filer(1:iblank-1)//'.'//tipo//filename(1:4)

      write(*,'(a)') ' Gravando '//filename(1:50)

      
      open (unit=iplt,file=filename,form='formatted')
      
      write (iplt,'(a,i5,1x,a)') 'Ensight Vector passo ',istep


! ***  

   

      
        write (iplt,'(1p,6e12.4E3)') (disp(id(1,i)), &
                                     disp(id(2,i)), &
                                     disp(id(3,i)), &
     				                  i=1,numnp)
   close (unit=iplt)
      
   return
   end
   subroutine ensg_vec_flux (tipo,fluxo,nume,istep)
	  
	  use modtapes,     only: filer, iblank, iplt
      
      implicit real*8 (a-h,o-z)
      

      REAL*8 , INTENT(in) :: fluxo(nume,3)
      real*8, allocatable ::  ap(:)
      character*70 filename
      character*1  tipo

!     arquivo Ensight para vetores de elementos

      write (filename(1:4),'(i4.4)') istep

      write(filename,'(a)') filer(1:iblank-1)//'.'//tipo//filename(1:4)

      write(*,'(a)') ' Gravando '//filename(1:50)

      
      open (unit=iplt,file=filename,form='formatted')
      
      write (iplt,'(a,i5,1x,a)') 'Ensight Vector passo ',istep

      write (iplt,'(a)') 'part 1'
      write (iplt,'(a)') 'tetra4'

! ***  

  

      
        write (iplt,'(1p,6e12.4E3)') (fluxo(iel,1), &
                                     fluxo(iel,2), &
                                     fluxo(iel,3), &
     				                  iel=1,nume)
   close (unit=iplt)
      
   return

   end 

!   subroutine  ensg_scalar_node (tipo,esca,numnp,etype,istep)
!	 use modtapes,     only: filer, iblank, iplt
!     use modmesh,      only: invp
!     use modvar ,      only: ireordering
!
!
!          
!     implicit real*8 (a-h,o-z)
!     
!     character*14 	etype
!
!     
!     dimension   esca     (1,numnp) 
!     
!     character*70 filename
!     character*3  tipo
!
!     write (filename(1:4),'(i4.4)') istep
!
!     write(filename,'(a)') './OUT2/'//filer(1:iblank-1)//'.'//tipo//filename(1:4)
!
! 
!     open (unit=iplt,file=filename,form='formatted')
!     
!     write(*,'(a)') ' Gravando '//filename(1:50)
!
!     
!     open (unit=iplt,file=filename,form='formatted')
!     
!
!
!     write (iplt,'(a,i5,1x,a)') 'Ensight Escalar passo ',istep
!     
!     ! if     (  etype.eq.'calor2d       ')                                      then     ! ***************  CALOR 2D  ******************** 
!     !
!     !   write (iplt,'(a/a)') 'part 1','tria3' 
!     !
!     ! elseif (  etype.eq.'elasticidade3d')                                      then     ! *************** ELASTICIDADE 3D *************** 
!     !
!     !   write (iplt,'(a/a)') 'part 1','tetra4'
!     !
!     ! elseif (  etype.eq.'tetra4_calor')                                      then     ! *************** CALOR 3D *************** 
!     !
!     !   write (iplt,'(a/a)') 'part 1','tetra4'
!     !
!     ! endif
!
!     IF (ireordering .eq. 0) then
!           write (iplt,'(1p,6e12.5)') (esca(1,i),i=1,numnp)
!     else
!           write (iplt,'(1p,6e12.5)') (esca(1,invp(i)),i=1,numnp)
!     endif
!           
!     
!     close (unit=iplt)
!     
!     return
!     end

     subroutine ensg_vec_flux_nodal (tipo,fluxo_nodal,numnp,istep)
	  
	  use modtapes,     only: filer, iblank, iplt
      
      implicit real*8 (a-h,o-z)
      

      REAL*8 , INTENT(in) :: fluxo_nodal(numnp,3)
      real*8, allocatable ::  ap(:)
      character*70 filename
      character*2  tipo

!     arquivo Ensight para vetores de elementos

      write (filename(1:4),'(i4.4)') istep

      write(filename,'(a)') filer(1:iblank-1)//'.'//tipo//filename(1:4)

      write(*,'(a)') ' Gravando '//filename(1:50)

      
      open (unit=iplt,file=filename,form='formatted')
      
      write (iplt,'(a,i5,1x,a)') 'Ensight Vector passo ',istep

      !write (iplt,'(a)') 'part 1'
      !write (iplt,'(a)') 'tetra4'

! ***  

  

      
        write (iplt,'(1p,6e12.4E3)') (fluxo_nodal(ind,1), &
                                     fluxo_nodal(ind,2), &
                                     fluxo_nodal(ind,3), &
     				                  ind=1,numnp)
   close (unit=iplt)
      
   return

   end 

   SUBROUTINE ensg_vec_shell (tipo,disp,up,id,ngl,numnp,neq,neqp,istep)
	  
use modtapes,     only: filer, iblank, iplt
      
implicit real*8 (a-h,o-z)


dimension   disp (0:neq), id   (ngl,numnp), up (0:neqp) 
real*8, allocatable ::  ap(:)
character*70 filename
character*3  tipo

write (filename(1:4),'(i4.4)') istep

write(filename,'(a)') filer(1:iblank-1)//'.'//tipo//filename(1:4)

!write(*,'(a)') ' Gravando '//filename(1:50)


open (unit=iplt,file=filename,form='formatted')

write (iplt,'(a,i5,1x,a)') 'Ensight Vector passo ',istep


 


if (tipo .eq. 'des') THEN
    write (iplt,'(1p,6es12.5)') (disp(id(1,i)), &
                                  disp(id(2,i)), &
                                  disp(id(3,i)), &
				       i=1,numnp)
elseif (tipo .eq. 'rot') THEN
    write (iplt,'(1p,6es12.5)') (disp(id(4,i)), &
                                  disp(id(5,i)), &
                                  disp(id(6,i)), &
				       i=1,numnp)
endif


close (unit=iplt)

open  (unit=iplt,file=filer(1:iblank-1)//'.case')

write (iplt,'(a)') 'FORMAT'
write (iplt,'(a)') 'type:    ensight'
write (iplt,'(a)')
write (iplt,'(a)') 'GEOMETRY'
write (iplt,'(a)') 'model:   '//filer(1:iblank-1)//'.GEO'
write (iplt,'(a)')
write (iplt,'(a)') 'VARIABLE'
write (iplt,'(a)') 'vector per node:     deslocamento    '//filer(1:iblank-1)//'.des0000'
write (iplt,'(a)') 'vector per node:     rotacao         '//filer(1:iblank-1)//'.rot0000'
write (iplt,*)

close (unit=iplt)
   
return

    ENDSUBROUTINE
    
       subroutine  ensg_scalar_node (tipo,u,up,id,numnp,neq,neqp,istep)
	 use modtapes,     only: filer, iblank, iplt
     use modmesh,      only: invp
     use modvar ,      only: ireordering

     implicit real*8 (a-h,o-z)
     
     dimension   u(0:neq), id(1,numnp), up(0:neqp) 
     real*8, allocatable ::  ap(:)
     
     character*70 filename
     character*3  tipo

     write (filename(1:4),'(i4.4)') istep
     write(filename,'(a)') './OUT2/'//filer(1:iblank-1)//'.'//tipo//filename(1:4)
     open (unit=iplt,file=filename,form='formatted')
     write(*,'(a)') ' Gravando '//filename(1:50)
     open (unit=iplt,file=filename,form='formatted')
     write (iplt,'(a,i5,1x,a)') 'Ensight Escalar passo ',istep
     
     allocate (ap(-neqp:neq)) ; ap = 0.d0
     
     if (neqp .gt. 0) then
        do i = 1, numnp 
          do j = 1, 1
             ieq = id(j,i)
	         if (ieq.ge.0) then
	           ap(ieq)=u(ieq)
             else
               ap(ieq)=up(-ieq)
	         endif
          enddo
        enddo
        IF (ireordering .eq. 0) then
              write (iplt,'(1p,6e12.5)') (ap(id(1,i)),i=1,numnp)
        else
              write (iplt,'(1p,6e12.5)') (ap(id(1,invp(i))),i=1,numnp)
        endif
     else
        IF (ireordering .eq. 0) then
              write (iplt,'(1p,6e12.5)') (u(i),i=1,numnp)
        else
              write (iplt,'(1p,6e12.5)') (u(invp(i)),i=1,numnp)
        endif
    endif
           
     
     close (unit=iplt)
     
     return
     end
