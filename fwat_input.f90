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
module fwat_input
use constants, only :MAX_STRING_LEN 

contains

  subroutine read_fwat_par_file(myrank)
  use fullwave_adjoint_tomo_par
  use my_mpi             !! module from specfem

  implicit none

  include "precision.h"  !! from specfem
  integer, intent(in) :: myrank
  ! locals
  character(len=MAX_STRING_LEN)                              :: line, keyw
  integer                                                    :: ipos0,ipos1,ier,icomp



  if (myrank == 0) then

     ! open par file
     VERBOSE_MODE = .false.
     TW_BEFORE=0.
     TW_AFTER=0.
     USE_NEAR_OFFSET=.true.
     USE_SPH_SMOOTH=.false.
     USE_PDE_SMOOTH=.false.
     KERNEL_TAPER = .false.
     open(666,file='fwat_params/FWAT.PAR') 
     do
          read(666,'(a)',end=99) line
          if (is_blank_line(line) .or. line(1:1) == '#') cycle

          !! INDICES TO READ line -----------------------------------------------
          ipos0=index(line,':')+1
          ipos1=index(line,'#')-1
          if (ipos1 < 0 ) ipos1=len_trim(line)

          !! STORE KEYWORD ITEM -------------------------------------------------
          keyw=trim(adjustl(line(1:ipos0-2)))

          !! DIFFERENT ITEM TO READ ---------------------------------------------
          select case (trim(keyw))
          case('VERBOSE_MODE')
            read(line(ipos0:ipos1),*) VERBOSE_MODE
          case('NSCOMP')
            read(line(ipos0:ipos1),*) NSCOMP
            allocate(SCOMPS(NSCOMP))
          case('SCOMPS')
            read(line(ipos0:ipos1),*) SCOMPS(1:NSCOMP)
          case('NRCOMP')
            read(line(ipos0:ipos1),*) NRCOMP
            allocate(RCOMPS(NRCOMP))
          case('RCOMPS')
            read(line(ipos0:ipos1),*) RCOMPS(1:NRCOMP)
            ! Mijian move dat_coord here as a public variable
            do icomp=1,NRCOMP
               if (trim(RCOMPS(icomp))=='R' .or. trim(RCOMPS(icomp))=='T' ) then
                  dat_coord='ZRT'
               else 
                  dat_coord='ZNE'
               endif
            enddo
          case('CH_CODE')
            read(line(ipos0:ipos1),*) CH_CODE
          case('NUM_FILTER')
             read(line(ipos0:ipos1),*) NUM_FILTER
             allocate(SHORT_P(NUM_FILTER))
             allocate(LONG_P(NUM_FILTER))
             allocate(GROUPVEL_MIN(NUM_FILTER))!BinHe added
             allocate(GROUPVEL_MAX(NUM_FILTER))!BinHe added
          case('SHORT_P')
             read(line(ipos0:ipos1),*) SHORT_P(:)
          case('LONG_P')
             read(line(ipos0:ipos1),*) LONG_P(:)
          case('GROUPVEL_MIN')!BinHe added
             read(line(ipos0:ipos1),*) GROUPVEL_MIN(:)
          case('GROUPVEL_MAX')!BinHe added
             read(line(ipos0:ipos1),*) GROUPVEL_MAX(:)
          case('TW_BEFORE')
             read(line(ipos0:ipos1),*) TW_BEFORE
          case('TW_AFTER')
             read(line(ipos0:ipos1),*) TW_AFTER
          case('ADJ_SRC_NORM')
             read(line(ipos0:ipos1),*) ADJ_SRC_NORM
          case('USE_NEAR_OFFSET')!BinHe added
             read(line(ipos0:ipos1),*) USE_NEAR_OFFSET
          case('SUPPRESS_EGF')
             read(line(ipos0:ipos1),*) SUPPRESS_EGF
          case('SAVE_OUTPUT_EACH_EVENT')
             read(line(ipos0:ipos1),*) SAVE_OUTPUT_EACH_EVENT
          case('USE_SPH_SMOOTH')
             read(line(ipos0:ipos1),*) USE_SPH_SMOOTH
          case('USE_PDE_SMOOTH')
             read(line(ipos0:ipos1),*) USE_PDE_SMOOTH
          case('SIGMA_H')
             read(line(ipos0:ipos1),*) SIGMA_H
          case('SIGMA_V')
             read(line(ipos0:ipos1),*) SIGMA_V
          case('OPT_METHOD')
             read(line(ipos0:ipos1),*) OPT_METHOD
          case('PRECOND_TYPE')
             read(line(ipos0:ipos1),*) PRECOND_TYPE
          case('DO_LS')
             read(line(ipos0:ipos1),*) DO_LS
          case('MAX_SLEN')
             read(line(ipos0:ipos1),*) MAX_SLEN
          case('NUM_STEP')
             read(line(ipos0:ipos1),*) NUM_STEP
             allocate(STEP_LENS(NUM_STEP))
          case('STEP_LENS')
             read(line(ipos0:ipos1),*) STEP_LENS(:)
          ! Mijian changed for RF
         !  case('NGAUSS')
         !     read(line(ipos0:ipos1),*) NGAUSS
         !     allocate(F0(NGAUSS))
         !  case('F0')
         !     read(line(ipos0:ipos1),*) F0(:)
         !  case('MINDERR')
         !     read(line(ipos0:ipos1),*) MINDERR
         !  case('RF_TSHIFT')
         !     read(line(ipos0:ipos1),*) RF_TSHIFT
         !  case('MAXIT')
         !     read(line(ipos0:ipos1),*) MAXIT
          case('KERNEL_TAPER')
             read(line(ipos0:ipos1),*) KERNEL_TAPER
          case('TAPER_H_SUPPRESS')
             read(line(ipos0:ipos1),*) TAPER_H_SUPPRESS
          case('TAPER_H_BUFFER')
             read(line(ipos0:ipos1),*) TAPER_H_BUFFER
          case('TAPER_V')
             read(line(ipos0:ipos1),*) TAPER_V
          case default
             write(*,*) 'ERROR KEY WORD NOT MATCH : ', trim(keyw), ' in FWAT.PAR '
             exit
          end select
       enddo

       ! close par file
99       close(666)
  
  endif ! of if (myrank == 0) then
  call MPI_BCAST(VERBOSE_MODE,1,MPI_LOGICAL,0,my_local_mpi_comm_world,ier)
  call MPI_BCAST(NSCOMP, 1,MPI_INTEGER,0,my_local_mpi_comm_world,ier)
  call MPI_BCAST(NRCOMP, 1,MPI_INTEGER,0,my_local_mpi_comm_world,ier)
  call MPI_BCAST(NUM_FILTER, 1,MPI_INTEGER,0,my_local_mpi_comm_world,ier)
  call MPI_BCAST(NUM_STEP, 1,MPI_INTEGER,0,my_local_mpi_comm_world,ier)
  if (myrank>0) then
     allocate(SCOMPS(NSCOMP)) 
     allocate(RCOMPS(NRCOMP)) 
     allocate(SHORT_P(NUM_FILTER))
     allocate(LONG_P(NUM_FILTER))
     allocate(GROUPVEL_MIN(NUM_FILTER))!BinHe added
     allocate(GROUPVEL_MAX(NUM_FILTER))!BinHe added
     allocate(STEP_LENS(NUM_STEP))
  endif
  call MPI_BCAST(SCOMPS,NSCOMP*MAX_STRING_LEN,MPI_CHARACTER,0,my_local_mpi_comm_world,ier)
  call MPI_BCAST(RCOMPS,NRCOMP*MAX_STRING_LEN,MPI_CHARACTER,0,my_local_mpi_comm_world,ier)
  call MPI_BCAST(dat_coord,MAX_STRING_LEN,MPI_CHARACTER,0,my_local_mpi_comm_world,ier)
  call MPI_BCAST(CH_CODE,MAX_STRING_LEN,MPI_CHARACTER,0,my_local_mpi_comm_world,ier)
  call MPI_BCAST(SHORT_P,NUM_FILTER,MPI_REAL,0,my_local_mpi_comm_world,ier)
  call MPI_BCAST(LONG_P,NUM_FILTER,MPI_REAL,0,my_local_mpi_comm_world,ier)
  call MPI_BCAST(GROUPVEL_MIN,NUM_FILTER,MPI_REAL,0,my_local_mpi_comm_world,ier)!BinHe added Kai changed int to real
  call MPI_BCAST(GROUPVEL_MAX,NUM_FILTER,MPI_REAL,0,my_local_mpi_comm_world,ier)!BinHe added Kai changed int to real
  call MPI_BCAST(TW_BEFORE,1,MPI_REAL,0,my_local_mpi_comm_world,ier)
  call MPI_BCAST(TW_AFTER,1,MPI_REAL,0,my_local_mpi_comm_world,ier)
  call MPI_BCAST(ADJ_SRC_NORM,1,MPI_LOGICAL,0,my_local_mpi_comm_world,ier)
  call MPI_BCAST(USE_NEAR_OFFSET,1,MPI_LOGICAL,0,my_local_mpi_comm_world,ier)
  call MPI_BCAST(SUPPRESS_EGF,1,MPI_LOGICAL,0,my_local_mpi_comm_world,ier)
  call MPI_BCAST(SAVE_OUTPUT_EACH_EVENT,1,MPI_LOGICAL,0,my_local_mpi_comm_world,ier)
  call MPI_BCAST(USE_SPH_SMOOTH,1,MPI_LOGICAL,0,my_local_mpi_comm_world,ier)
  call MPI_BCAST(SIGMA_H,1,MPI_REAL,0,my_local_mpi_comm_world,ier)
  call MPI_BCAST(SIGMA_V,1,MPI_REAL,0,my_local_mpi_comm_world,ier)
  call MPI_BCAST(OPT_METHOD,MAX_STRING_LEN,MPI_CHARACTER,0,my_local_mpi_comm_world,ier)
  call MPI_BCAST(PRECOND_TYPE,MAX_STRING_LEN,MPI_CHARACTER,0,my_local_mpi_comm_world,ier)
  call MPI_BCAST(DO_LS,1,MPI_LOGICAL,0,my_local_mpi_comm_world,ier)
  call MPI_BCAST(MAX_SLEN,1,MPI_REAL,0,my_local_mpi_comm_world,ier)
  call MPI_BCAST(STEP_LENS,NUM_STEP,MPI_REAL,0,my_local_mpi_comm_world,ier)
  ! Mijian add for RF
!   call MPI_BCAST(NGAUSS, 1,MPI_INTEGER,0,my_local_mpi_comm_world,ier)
!   if (myrank>0) then
!       allocate(F0(NGAUSS))
!   endif
!   call MPI_BCAST(F0,NGAUSS,MPI_REAL,0,my_local_mpi_comm_world,ier)
!   call MPI_BCAST(rf_tshift,1,MPI_REAL,0,my_local_mpi_comm_world,ier)
!   call MPI_BCAST(maxit,1,MPI_INTEGER,0,my_local_mpi_comm_world,ier)
!   call MPI_BCAST(minderr,1,MPI_REAL,0,my_local_mpi_comm_world,ier)
  call MPI_BCAST(KERNEL_TAPER,1,MPI_LOGICAL,0,my_local_mpi_comm_world,ier)
  call MPI_BCAST(TAPER_H_SUPPRESS,1,MPI_REAL,0,my_local_mpi_comm_world,ier)
  call MPI_BCAST(TAPER_H_BUFFER,1,MPI_REAL,0,my_local_mpi_comm_world,ier)
  call MPI_BCAST(TAPER_V,1,MPI_REAL,0,my_local_mpi_comm_world,ier)

end subroutine read_fwat_par_file

!####################################################################
!----------------------------------------------------------------
! true if blank line or commented line (ie that begins with #)
!----------------------------------------------------------------
  logical function is_blank_line(line)
    
    character(len=MAX_STRING_LEN), intent(in) :: line
    is_blank_line=.false.
    if (len(trim(adjustl(line))) == 0) is_blank_line = .true.
    if (INDEX(trim(adjustl(line)),'#') == 1) is_blank_line = .true.
  end function is_blank_line
!####################################################################yy

end module fwat_input
