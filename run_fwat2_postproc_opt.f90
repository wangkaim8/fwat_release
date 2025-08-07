!=====================================================================
!
!               Full Waveform Adjoint Tomography -v1.1
!               ---------------------------------------
!
!     Main historical authors: Kai Wang
!                              Macquarie Uni, Australia
!                            & University of Toronto, Canada
!                           (c) Martch 2020
!
!=====================================================================
!

 subroutine run_fwat2_postproc_opt(model,evtsetb,evtsete,is_smooth)

  use fullwave_adjoint_tomo_par
  use fwat_input
  use tomography_par, only: USE_ALPHA_BETA_RHO,USE_ISO_KERNELS, step_fac, iter_start,iter_current, &
           INPUT_MODEL_DIR,OUTPUT_MODEL_DIR,KERNEL_OLD_DIR,INPUT_KERNELS_DIR,PRINT_STATISTICS_FILES
  !!! for reading Database
  use specfem_par
  use specfem_par_elastic, only: ELASTIC_SIMULATION,ispec_is_elastic,rho_vp,rho_vs,min_resolved_period
  use specfem_par_acoustic, only: ACOUSTIC_SIMULATION,ispec_is_acoustic
  use specfem_par_poroelastic, only: POROELASTIC_SIMULATION,ispec_is_poroelastic,rho_vpI,rho_vpII,rho_vsI, &
    phistore,tortstore,rhoarraystore
  !!! for model update
  use tomography_model_iso
  use tomography_kernels_iso
  use taper3d
  use constants, only: IMAIN
  implicit none

  real(kind=CUSTOM_REAL) :: distance_min_glob,distance_max_glob
  real(kind=CUSTOM_REAL) :: elemsize_min_glob,elemsize_max_glob
  real(kind=CUSTOM_REAL) :: x_min_glob,x_max_glob
  real(kind=CUSTOM_REAL) :: y_min_glob,y_max_glob
  real(kind=CUSTOM_REAL) :: z_min_glob,z_max_glob


  character(len=MAX_STRING_LEN)                   :: model,model_prev,model_next 
  character(len=MAX_STRING_LEN)                   :: evtsetb,evtsete,is_smooth 
  character(len=MAX_STRING_LEN)                   :: strset,strstep
  character(len=MAX_STRING_LEN) :: input_dir,output_dir,ekernel_dir_list
  character(len=MAX_STRING_LEN) :: kernel_names_comma_delimited
  integer iset, setb, sete,ier
  !real :: t1, t2
  logical :: USE_GPU
  integer :: imod_current,imod_up, imod_down
  integer :: istep
  type(taper_cls) :: taper_3d
    
  read(model(2:),'(I2.2)') imod_current
  imod_up=imod_current+1
  imod_down=imod_current-1
  write(model_next,'(A1,I2.2)') 'M',imod_up
  if (imod_down <0) then
     model_prev='none'
  else
     write(model_prev,'(A1,I2.2)') 'M',imod_down
  endif
  read(evtsetb(4:),'(I3)') setb 
  read(evtsete(4:),'(I3)') sete
  ekernel_dir_list='optimize/ekernel_dir.lst'
  output_dir='optimize/SUM_KERNELS_'//trim(model)
  !**************** Build directories for the storage ****************************
  call read_fwat_par_file(myrank)
  if (myrank == 0) then
     open(unit=OUT_FWAT_LOG,file='output_fwat2_log_'//trim(model)//'.'//trim(evtsetb)&
                                 //'-'//trim(evtsete)//'.txt')
     write(OUT_FWAT_LOG,*) 'Running XFWAT2_POSTPROC_OPT !!!'
     write(OUT_FWAT_LOG,*) 'model,evtsetb,evtsete: ',trim(model),' ', trim(evtsetb),'-',trim(evtsete)
     write(OUT_FWAT_LOG,*) 'SAVE_OUTPUT_EACH_EVENT: ',SAVE_OUTPUT_EACH_EVENT
     write(OUT_FWAT_LOG,*) 'USE_SPH_SMOOTH: ',USE_SPH_SMOOTH
     write(OUT_FWAT_LOG,*) 'SIGMA_H, SIGMA_V: ',SIGMA_H, SIGMA_V
     write(OUT_FWAT_LOG,*) 'model_prev,model_next: ',trim(model_prev),' ',trim(model_next)
     write(OUT_FWAT_LOG,*) 'OPT_METHOD: ',trim(OPT_METHOD) 
     write(OUT_FWAT_LOG,*) 'DO_LS: ',DO_LS 
     write(OUT_FWAT_LOG,*) 'MAX_SLEN: ',MAX_SLEN
     write(OUT_FWAT_LOG,*) 'STEP_LENS: ',STEP_LENS 
     call system('mkdir -p optimize')
     if(imod_current==0) then
       write(*,*) 'cp -r model_initial optimize/MODEL_M00'
       call system('mkdir -p optimize/MODEL_M00')
       call system('cp model_initial/* optimize/MODEL_M00')
     endif
     if (DO_LS) then
       do istep=1,NUM_STEP
         step_fac=STEP_LENS(istep)
         write(strstep,'(F5.3)') step_fac
         call system('mkdir -p optimize/MODEL_'//trim(model)//'_step'//trim(strstep))
       enddo
     else
         call system('mkdir -p optimize/MODEL_'//trim(model_next))
     endif
     call system('mkdir -p '//trim(output_dir))
    
     open(unit=400,file=trim(ekernel_dir_list),iostat=ier)
     do iset=setb,sete
        write(strset,'(I0)') iset
        write(400,'(a)')'solver/'//trim(model)//'.set'//trim(strset)//'/GRADIENT'
     enddo
     close(400)
     write(OUT_FWAT_LOG,*) '*******************************************************'
    
  endif 
  call synchronize_all()
! -----------------------------------------------------------------
! precond and sum 
  if(myrank==0) write(OUT_FWAT_LOG,*) 'This is sum_precond ...'
  call sum_preconditioned_kernels_fwat(ekernel_dir_list,output_dir)
! -----------------------------------------------------------------
! smooth
  ! Kai: Maybe we should move reding Database back to the smoothing subroutine.
  !!!  Read Database before smoothing !!! 
  ! read the value of NSPEC_AB and NGLOB_AB because we need it to define some array sizes below
  call read_mesh_for_init()

  allocate(ibool(NGLLX,NGLLY,NGLLZ,NSPEC_AB),irregular_element_number(NSPEC_AB),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 980')

  if (NSPEC_IRREGULAR > 0) then
    ! allocate arrays for storing the databases
    allocate(xix(NGLLX,NGLLY,NGLLZ,NSPEC_IRREGULAR),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 981')
    allocate(xiy(NGLLX,NGLLY,NGLLZ,NSPEC_IRREGULAR),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 982')
    allocate(xiz(NGLLX,NGLLY,NGLLZ,NSPEC_IRREGULAR),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 983')
    allocate(etax(NGLLX,NGLLY,NGLLZ,NSPEC_IRREGULAR),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 984')
    allocate(etay(NGLLX,NGLLY,NGLLZ,NSPEC_IRREGULAR),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 985')
    allocate(etaz(NGLLX,NGLLY,NGLLZ,NSPEC_IRREGULAR),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 986')
    allocate(gammax(NGLLX,NGLLY,NGLLZ,NSPEC_IRREGULAR),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 987')
    allocate(gammay(NGLLX,NGLLY,NGLLZ,NSPEC_IRREGULAR),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 988')
    allocate(gammaz(NGLLX,NGLLY,NGLLZ,NSPEC_IRREGULAR),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 989')
    allocate(jacobian(NGLLX,NGLLY,NGLLZ,NSPEC_IRREGULAR),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 990')
   else
       ! allocate arrays for storing the databases
    allocate(xix(1,1,1,1),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 991')
    allocate(xiy(1,1,1,1),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 992')
    allocate(xiz(1,1,1,1),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 993')
    allocate(etax(1,1,1,1),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 994')
    allocate(etay(1,1,1,1),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 995')
    allocate(etaz(1,1,1,1),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 996')
    allocate(gammax(1,1,1,1),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 997')
    allocate(gammay(1,1,1,1),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 998')
    allocate(gammaz(1,1,1,1),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 999')
    allocate(jacobian(1,1,1,1),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1000')
  endif

  ! mesh node locations
  allocate(xstore(NGLOB_AB), &
           ystore(NGLOB_AB), &
           zstore(NGLOB_AB),stat=ier)
  if (ier /= 0) stop 'Error allocating arrays for mesh nodes'

  ! material properties
  allocate(kappastore(NGLLX,NGLLY,NGLLZ,NSPEC_AB), &
           mustore(NGLLX,NGLLY,NGLLZ,NSPEC_AB),stat=ier)
  if (ier /= 0) stop 'Error allocating arrays for material properties'

  ! material flags
  allocate(ispec_is_acoustic(NSPEC_AB), &
           ispec_is_elastic(NSPEC_AB), &
           ispec_is_poroelastic(NSPEC_AB),stat=ier)
  if (ier /= 0) stop 'Error allocating arrays for material flags'
  ispec_is_acoustic(:) = .false.
  ispec_is_elastic(:) = .false.
  ispec_is_poroelastic(:) = .false.

  ! reads in external mesh
  call read_mesh_databases()

  ! gets mesh dimensions
  call check_mesh_distances(myrank,NSPEC_AB,NGLOB_AB,ibool,xstore,ystore,zstore, &
                            x_min_glob,x_max_glob,y_min_glob,y_max_glob,z_min_glob,z_max_glob, &
                            elemsize_min_glob,elemsize_max_glob, &
                            distance_min_glob,distance_max_glob)

  ! outputs infos
  if (myrank == 0) then
    print *,'mesh dimensions:'
    print *,'  Xmin and Xmax of the model = ',x_min_glob,x_max_glob
    print *,'  Ymin and Ymax of the model = ',y_min_glob,y_max_glob
    print *,'  Zmin and Zmax of the model = ',z_min_glob,z_max_glob
    print *
    print *,'  Max GLL point distance = ',distance_max_glob
    print *,'  Min GLL point distance = ',distance_min_glob
    print *,'  Max/min ratio = ',distance_max_glob/distance_min_glob
    print *
    print *,'  Max element size = ',elemsize_max_glob
    print *,'  Min element size = ',elemsize_min_glob
    print *,'  Max/min ratio = ',elemsize_max_glob/elemsize_min_glob
    print *
  endif
  if (myrank == 0 .and. IMAIN /= 6) &
    open(unit=IMAIN,file=trim(OUTPUT_FILES)//'/output_mesh_resolution.txt',status='unknown')

  if (ELASTIC_SIMULATION) then
    call check_mesh_resolution(myrank,NSPEC_AB,NGLOB_AB, &
                               ibool,xstore,ystore,zstore, &
                               kappastore,mustore,rho_vp,rho_vs, &
                               DT,model_speed_max,min_resolved_period, &
                               LOCAL_PATH,SAVE_MESH_FILES)

  else if (POROELASTIC_SIMULATION) then
    allocate(rho_vp(NGLLX,NGLLY,NGLLZ,NSPEC_AB))
    allocate(rho_vs(NGLLX,NGLLY,NGLLZ,NSPEC_AB))
    rho_vp = 0.0_CUSTOM_REAL
    rho_vs = 0.0_CUSTOM_REAL
    call check_mesh_resolution_poro(myrank,NSPEC_AB,NGLOB_AB,ibool,xstore,ystore,zstore, &
                                    DT,model_speed_max,min_resolved_period, &
                                    phistore,tortstore,rhoarraystore,rho_vpI,rho_vpII,rho_vsI, &
                                    LOCAL_PATH,SAVE_MESH_FILES)
    deallocate(rho_vp,rho_vs)
  else if (ACOUSTIC_SIMULATION) then
    allocate(rho_vp(NGLLX,NGLLY,NGLLZ,NSPEC_AB),stat=ier)
    if (ier /= 0) stop 'Error allocating array rho_vp'
    allocate(rho_vs(NGLLX,NGLLY,NGLLZ,NSPEC_AB),stat=ier)
    if (ier /= 0) stop 'Error allocating array rho_vs'
    rho_vp = sqrt( kappastore / rhostore ) * rhostore
    rho_vs = 0.0_CUSTOM_REAL
    call check_mesh_resolution(myrank,NSPEC_AB,NGLOB_AB, &
                               ibool,xstore,ystore,zstore, &
                               kappastore,mustore,rho_vp,rho_vs, &
                               DT,model_speed_max,min_resolved_period, &
                               LOCAL_PATH,SAVE_MESH_FILES)
    deallocate(rho_vp,rho_vs)
  endif
  !!!        end reading Database          !!!
  if(trim(is_smooth)=='true')then 
  input_dir=trim(output_dir)//'/'
  USE_GPU=.false.
  if(myrank==0) write(OUT_FWAT_LOG,*) 'This is smoothing ...'

  ! Mijian add for taperring misfit kernel
  
  call bcast_all_singlecr(x_min_glob)
  call bcast_all_singlecr(x_max_glob)
  call bcast_all_singlecr(y_min_glob)
  call bcast_all_singlecr(y_max_glob)
  call bcast_all_singlecr(z_min_glob)
  call bcast_all_singlecr(z_max_glob)

  if (KERNEL_TAPER) then
    call taper_3d%create(x_min_glob,x_max_glob,y_min_glob,y_max_glob,z_min_glob,z_max_glob, &
                         elemsize_min_glob/2., elemsize_min_glob/2, elemsize_min_glob/2, &
                         TAPER_H_SUPPRESS,TAPER_H_BUFFER,TAPER_V)
  endif

  if (USE_ISO_KERNELS) then
    if (USE_PDE_SMOOTH) then
      kernel_names_comma_delimited='bulk_c_kernel'
      call smooth_sem_pde_fwat(sigma_h,sigma_v,kernel_names_comma_delimited,input_dir,&
                     output_dir,USE_GPU)
      kernel_names_comma_delimited='bulk_beta_kernel'
      call smooth_sem_pde_fwat(sigma_h,sigma_v,kernel_names_comma_delimited,input_dir,&
                     output_dir,USE_GPU)
      kernel_names_comma_delimited='rho_kernel'
      call smooth_sem_pde_fwat(sigma_h,sigma_v,kernel_names_comma_delimited,input_dir,&
                     output_dir,USE_GPU)
    else
      kernel_names_comma_delimited='bulk_c_kernel'
      call smooth_sem_fwat(sigma_h,sigma_v,kernel_names_comma_delimited,input_dir,&
                      output_dir,USE_GPU, taper_3d)
      kernel_names_comma_delimited='bulk_beta_kernel'
      call smooth_sem_fwat(sigma_h,sigma_v,kernel_names_comma_delimited,input_dir,&
                      output_dir,USE_GPU, taper_3d)
      kernel_names_comma_delimited='rho_kernel'
      call smooth_sem_fwat(sigma_h,sigma_v,kernel_names_comma_delimited,input_dir,&
                      output_dir,USE_GPU, taper_3d)
    endif ! end if use_pde_smooth
  else if (USE_ALPHA_BETA_RHO) then
    if (USE_SPH_SMOOTH) then ! maybe we don't need to keep this
      kernel_names_comma_delimited='alpha_kernel'
      call smooth_sem_sph_pde_fwat(sigma_h,sigma_v,kernel_names_comma_delimited,input_dir,&
                      output_dir,USE_GPU)
      kernel_names_comma_delimited='beta_kernel'
      call smooth_sem_sph_pde_fwat(sigma_h,sigma_v,kernel_names_comma_delimited,input_dir,&
                      output_dir,USE_GPU)
      kernel_names_comma_delimited='rhop_kernel'
      call smooth_sem_sph_pde_fwat(sigma_h,sigma_v,kernel_names_comma_delimited,input_dir,&
                      output_dir,USE_GPU)
    else if (USE_PDE_SMOOTH) then
      kernel_names_comma_delimited='alpha_kernel'
      call smooth_sem_pde_fwat(sigma_h,sigma_v,kernel_names_comma_delimited,input_dir,&
                     output_dir,USE_GPU)
      kernel_names_comma_delimited='beta_kernel'
      call smooth_sem_pde_fwat(sigma_h,sigma_v,kernel_names_comma_delimited,input_dir,&
                     output_dir,USE_GPU)
      kernel_names_comma_delimited='rhop_kernel'
      call smooth_sem_pde_fwat(sigma_h,sigma_v,kernel_names_comma_delimited,input_dir,&
                     output_dir,USE_GPU)
    else
      kernel_names_comma_delimited='alpha_kernel'
      call smooth_sem_fwat(sigma_h,sigma_v,kernel_names_comma_delimited,input_dir,&
                     output_dir,USE_GPU, taper_3d)
      kernel_names_comma_delimited='beta_kernel'
      call smooth_sem_fwat(sigma_h,sigma_v,kernel_names_comma_delimited,input_dir,&
                     output_dir,USE_GPU, taper_3d)
      kernel_names_comma_delimited='rhop_kernel'
      call smooth_sem_fwat(sigma_h,sigma_v,kernel_names_comma_delimited,input_dir,&
                     output_dir,USE_GPU, taper_3d)
    endif ! endif use (pde_sph/pde/conv) smooth
  else
    if (USE_PDE_SMOOTH) then
      kernel_names_comma_delimited='bulk_c_kernel'
      call smooth_sem_pde_fwat(sigma_h,sigma_v,kernel_names_comma_delimited,input_dir,&
                     output_dir,USE_GPU)
      kernel_names_comma_delimited='bulk_betav_kernel'
      call smooth_sem_pde_fwat(sigma_h,sigma_v,kernel_names_comma_delimited,input_dir,&
                     output_dir,USE_GPU)
      kernel_names_comma_delimited='bulk_betah_kernel'
      call smooth_sem_pde_fwat(sigma_h,sigma_v,kernel_names_comma_delimited,input_dir,&
                     output_dir,USE_GPU)
      kernel_names_comma_delimited='eta_kernel'
      call smooth_sem_pde_fwat(sigma_h,sigma_v,kernel_names_comma_delimited,input_dir,&
                     output_dir,USE_GPU)
    else
      kernel_names_comma_delimited='bulk_c_kernel'
      call smooth_sem_fwat(sigma_h,sigma_v,kernel_names_comma_delimited,input_dir,&
                      output_dir,USE_GPU, taper_3d)
      kernel_names_comma_delimited='bulk_betav_kernel'
      call smooth_sem_fwat(sigma_h,sigma_v,kernel_names_comma_delimited,input_dir,&
                      output_dir,USE_GPU, taper_3d)
      kernel_names_comma_delimited='bulk_betah_kernel'
      call smooth_sem_fwat(sigma_h,sigma_v,kernel_names_comma_delimited,input_dir,&
                      output_dir,USE_GPU, taper_3d)
      kernel_names_comma_delimited='eta_kernel'
      call smooth_sem_fwat(sigma_h,sigma_v,kernel_names_comma_delimited,input_dir,&
                      output_dir,USE_GPU, taper_3d)
    endif !endif use_pde_smooth
  endif
  endif ! end if is_smooth
! -----------------------------------------------------------------
! optimize (SD, CG or LBFGS)
  if(myrank==0) write(OUT_FWAT_LOG,*) 'This is optimize ...'
  INPUT_MODEL_DIR='optimize/MODEL_'//trim(model)//'/'
  if (imod_down<0) then
     KERNEL_OLD_DIR='none'
  else
     KERNEL_OLD_DIR='optimize/SUM_KERNELS_'//trim(model_prev)//'/'
  endif
  INPUT_KERNELS_DIR='optimize/SUM_KERNELS_'//trim(model)//'/'
  PRINT_STATISTICS_FILES=.false.
  INPUT_DATABASES_DIR=trim(LOCAL_PATH)//'/'
  if (DO_LS) then
    do istep=1,NUM_STEP
      step_fac=STEP_LENS(istep)
      iter_start=0
      iter_current=imod_current
      write(strstep,'(F5.3)') step_fac
      OUTPUT_MODEL_DIR='optimize/MODEL_'//trim(model)//'_step'//trim(strstep)//'/'
      call model_update_opt()
    enddo
  else
    step_fac=MAX_SLEN
    iter_start=0
    iter_current=imod_current
    OUTPUT_MODEL_DIR='optimize/MODEL_'//trim(model_next)//'/'
    call model_update_opt() ! run first time to save model 
    OUTPUT_MODEL_DIR=trim(LOCAL_PATH)//'/' ! run second time to be called for next iteration
    call model_update_opt()
  endif 
! -----------------------------------------------------------------
  
  if(myrank==0)  write(OUT_FWAT_LOG,*) '*******************************************************'
  if(myrank==0)  write(OUT_FWAT_LOG,*) 'Finished FWAT stage2 here!!!'
  if(myrank==0) close(OUT_FWAT_LOG)
  if(myrank==0) close(IMAIN)


end subroutine run_fwat2_postproc_opt

