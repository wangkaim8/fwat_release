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
subroutine fullwave_adjoint_tomo_main()

  use fullwave_adjoint_tomo_par

  implicit none

  
  character(len=MAX_STRING_LEN)                   :: model 
  character(len=MAX_STRING_LEN)                   :: evtset 
  character(len=MAX_STRING_LEN)                   :: simu_type 
  character(len=MAX_STRING_LEN)                   :: run_opt
  integer :: myrank


  !!!##############################################################################################################################
  !!! ---------------------------------------------- INITIALIZE RUNTIME ----------------------------------------------------------
  call world_rank(myrank)  
  !!!##############################################################################################################################
  if (command_argument_count() /= 4) then
     if (myrank == 0) then
       print *,'USAGE:  mpirun -np NPROC bin/xfullwave_adjoint_tomo model set simu_type run_opt'
       print *,'model     --- name of current model, such as M00, M01, ...'
       print *,'set       --- name of event set, such as set1, set2, ...'
       print *,'simu_type --- simulation type: noise, tele'
       print *,'run_opt   --- run_opt= 1 (fwd), 2 (fwd_meas), 3 (fwd_meas_adj)'
       stop 'Please check command line arguments'
     endif
   endif

  call get_command_argument(1, model)     
  call get_command_argument(2, evtset)     
  call get_command_argument(3, simu_type)     
  call get_command_argument(4, run_opt)     


  !!!##############################################################################################################################
  !!! -------------------------------  different running mode : forward or FWI ----------------------------------------------------
  !!!##############################################################################################################################
  call run_fwat1_fwd_measure_adj(model,evtset,simu_type,run_opt)
  !call run_fwat2_postproc_opt(model,evtset,evtset)

end subroutine fullwave_adjoint_tomo_main

