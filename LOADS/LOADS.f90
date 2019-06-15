      SUBROUTINE LOADS()
	
!     LEITURA DOS DADOS DA MALHA E DEFINICAO DE VARIOS ARRANJOS

      use    modvar,    only: numnp,nume,nd,nnoel,ngl,neq,neqp, ireordering,nummat,ncprop,grav,etype,neq_bar,neq_PHI,neqp_PHI
      use    modmesh,   only: lm,incid,mtype,id,invp,x,y,z,prop,fno,found,lm_bar,lm_tetra,id_phi,lm_phi
      use    modtapes,  only: iin3, iout, ierror
      use    modloads,  only: fbody,f,up
	
      IMPLICIT NONE            

      !if (etype.eq.'GravNavStok2D') then
      !      CALL RHSV_READ (NUMNP,NGL+1,IIN3,IOUT,IERROR,invp,ireordering,fno,found)
      !else
      CALL RHSV_READ (NUMNP,NGL,IIN3,IOUT,IERROR,invp,ireordering,fno,found)
      !endif
      
      !if (etype .ne. 'tdpflex') then
      CALL RHSV_MAKE (LM,INCID,ID,NNOEL,ND,NUME,NUMNP,NEQ,NEQP,NGL,IIN3,IOUT,IERROR,fno,found,invp)
      CALL RHSV_LOAD (ID,NUMNP,NEQ,NEQP,NGL,IIN3,IOUT,IERROR,fno)
      
      !if (etype.eq.'GravNavStok2D') then
      !      CALL RHSV_MAKE_phi (lm_PHI,INCID,id_PHI,NNOEL,NUME,NUMNP,NEQ_PHI,NEQP_PHI,ngl,IIN3,IOUT,IERROR,fno,found,invp)
      !      CALL RHSV_LOAD_phi (id_PHI,NUMNP,NEQ_PHI,NEQP_PHI,ngl,IIN3,IOUT,IERROR,fno)
      !endif

      !else  
      !CALL RHSV_MAKE2 (LM_tetra,LM_bar,INCID,ID,NNOEL,ND,NUME,NUMNP,NEQ,NEQP,NGL,IIN3,IOUT,IERROR,fno,found,mtype,neq_bar,invp)
      !CALL RHSV_LOAD (ID,NUMNP,NEQ,NEQP,NGL,IIN3,IOUT,IERROR,fno)   
      !endif
      
      !DEALLOCATE (Fno)
      DEALLOCATE (Found)      

      CALL ALLOC(2)

      !if (etype.eq.'epd       ' .or. etype.eq.'ept       ' .or. etype.eq.'epdmohrc       ') then	
      !CALL tria3_elast_body(fbody,lm,x,y,incid,prop,mtype,numnp,nume,nummat,ncprop,nd,nnoel,neq,grav)   
      !
      !elseif     (etype.eq.'tria3_plate        ')            then    
      !CALL tria3_plate_body (fbody,lm,x,y,incid,mtype,prop,numnp,nume,nummat,nd,nnoel,ncprop,neq,grav)
      !
      !elseif     (etype.eq.'truss3d ' )            then    
      !CALL truss3d_body (fbody,lm,x,y,z,incid,mtype,prop,numnp,nume,nummat,nd,nnoel,ncprop,neq,grav)
      !
      !elseif     (etype.eq.'elasticidade3d' .or. etype .eq. "cap3d" )            then    
      !call tetra4_elast_body (fbody,lm,x,y,z,incid,mtype,prop,numnp,nume,nummat,nd,nnoel,ncprop,neq,etype,grav)
      !
      !elseif     (etype.eq.'tdpflex ' )            then    
      !CALL truss3d_body (fbody,lm_bar,x,y,z,incid,mtype,prop,numnp,nume,nummat,6,nnoel,ncprop,neq,grav) 
      !call tetra4_elast_body (fbody,lm,x,y,z,incid,mtype,prop,numnp,nume,nummat,nd,nnoel,ncprop,neq,etype,grav)         
      !
      !else
      fbody(:) = 0.d0

      !endif  

    
      RETURN
      END
      



   
