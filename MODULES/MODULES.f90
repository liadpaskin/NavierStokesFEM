!  **************** BEGIN MODVAR ****************
module modMPIvar

      implicit none
      include 'mpif.h'

        integer ::  numprocs, myid, numnp_proc, numnpi, numnpf, nume_proc, numei, numef, iroot, idebug,ierr
        integer, allocatable :: procsnp(:), procsel(:)


      logical ::    l_root, l_debug

        integer, dimension (MPI_STATUS_SIZE) :: status

 end module modMPIvar

!  **************** BEGIN MODVAR ****************
module modvar

implicit none
        
        integer ::  nd           ,  nnoel        ,  numnp        , &
                    nume         ,  nummat       ,  ngl          , & 
                    neq          ,  neqp = 0     ,  ncprop = 22  , &
                    nwk          ,  noplot       ,  nglplot      , &
                    isolver_type ,  itera        ,  ia           , &
                    ib           ,  ndaux        ,  ireordering  , &
                    nnomax       ,  nad          ,  lmax         , &
                    kmax         ,  nelblk       ,  kmax1        , &
                    niter        ,                  cohes        , &
                    nsubi        ,  kplast       ,  ntria        , &
                    nbar         ,  neq_bar      ,  npoint       , &
                    nsize        ,&
                    neq_PHI          ,  neqp_PHI = 0    , nbody, nlinemax, elrefp

        
        real*8  ::  value_mem    ,  energy ,  etol , grav , areat , ud , Eyoung , fac_dp , xxx , ccc, refPressure
        character*80 	hed
        character*14 	etype
           
 end module modvar 
!  **************** END MODVAR ******************


!  **************** BEGIN modtim ****************
module modtim

implicit none
        
        integer ::  IFUNC , NPTF  , KPFOUN  , NIMP, NFLAG
        
        integer*8 :: it1, it2, clockcount1, clockcount2, nsteps
        
        real*8  ::  TIMEF , DT , TFUNC , DAMP1, DAMP2 , ALPHA , BETA , GAMMA ,             & 
                    A_FUNC,  B_FUNC,  OMEGA_FUNC , PHASE_FUNC, TFTIM(200) , TFVAL(200) , CFLmax, CFLmax_ant

        real*8, allocatable ::   upred (:),  u_n (:),  v_n (:), fint_pred(:) , vpred(:) , a_n(:) , du(:)
        real*8, allocatable ::   upred_PHI (:),  u_n_PHI (:)
           
 end module modtim 
!  **************** END modtim ******************


!  **************** BEGIN MODTAPES **************
module modtapes

implicit none
      
        integer ::   IIN1  ,  IIN2  ,  IIN3   ,   IIN4, &
                     IERROR,  IOUT  ,  iblank , &
                     iplot ,  IPLT 
                     
         character*70 filer, fileinput
                    
end module modtapes
!  **************** END MODTAPES ****************


!  **************** BEGIN MODMESH  **************

module modmesh

            implicit none
            
            logical*1, allocatable :: &
                     found(:,:)

            integer, allocatable :: &
                     id  (:,:),  lm(:,:),  incid (:,:), incidaux2 (:,:), mtype (:) ,  jdiag(:) , ip(:) , ips(:) , &
                     noviz(:) , mask(:) , xls(:) , invp(:) , invp2(:), lm_tetra(:,:) , lm_bar(:,:) , list_tria3(:,:) , &
                     id_phi (:,:), lm_phi (:,:), bodys(:,:)
                     
                     
            real*8, allocatable ::  &
                     x (:),  y (:),  z (:), xi (:),  yi (:),  zi (:), prop (:,:) , Fno(:,:,:)
         
end module modmesh 
!  **************** END MODMESH  ****************


!  **************** BEGIN MODLOADS **************
module modloads

   implicit none

   real*8, allocatable ::  &
                        f(:),    up(:) , fbody(:) , &
                        f_PHI(:),    up_PHI(:)
   
       
end module modloads 
!  **************** END MODLOADS ****************


!  **************** BEGIN MODSOLVER *************
module modsolver
   implicit none
        
   INTEGER*8              :: nALHS!, nALHS_PHI
            
   real*8, allocatable :: &
                    stiff(:,:),      fp(:),    u(:)   , &
                 v(:) , diag(:) , ve1(:,:) , ve2(:,:) , &
                 Stiff_Convec(:,:),Stiff_difus_Lambda(:,:),Stiff_Difus_xMu(:,:), &
                 Forca(:,:), Fint(:,:), sig11(:), sig12(:), sig13(:), sig22(:), sig23(:), sig33(:) , &
                 stiff_PHI(:,:),      fp_PHI(:),    PHI(:), &
                 ALHS(:), CLHS(:), ALHS_PHI(:), CLHS_PHI(:), SedPHI(:,:), fpAcel(:), elPe(:)
                                
   integer, allocatable :: &                               
              maxa(:),  MHT(:) , lm_aux(:,:) , nsub(:) , iplast(:) , mht_aux(:) , ielblk(:)

                                
end module modsolver

module testecrout
    implicit none
    
    real*8, allocatable :: Alhs(:), Clhs(:)
    
    endmodule
    
module modwave
    implicit none
    
    integer :: numeqIn, nsup
    real*8 :: CompOnda, Amp, Depht, numOnda, angFreq
    integer, allocatable :: idsup(:), IDinlet(:,:)
    real*8, allocatable :: xsup(:), ysup(:)
    
endmodule
