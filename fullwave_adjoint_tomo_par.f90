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

module fullwave_adjoint_tomo_par

  !! IMPORT VARIABLES ------------------------------------------------------------------------------------------------
  use specfem_par, only: CUSTOM_REAL, MAX_STRING_LEN, MAX_LENGTH_STATION_NAME, MAX_LENGTH_NETWORK_NAME
  !-------------------------------------------------------------------------------------------------------------------

  implicit none

  !! ------------------------------ compialtion config parameters -----------------------------------------------------------------

  !! log file for inversion
  integer,                       public, parameter  :: OUT_FWAT_LOG=888
  real(kind=CUSTOM_REAL),        public, parameter  :: deg2rad = 0.017453292519943
!################################################# WORKFLOW ######################################################################
 character(len= MAX_STRING_LEN), public, dimension(:),       allocatable  :: station_file
 character(len= MAX_STRING_LEN), public, dimension(:),       allocatable  :: src_solution_file
 character(len= MAX_STRING_LEN), public, dimension(:),       allocatable  :: evtid_names
 character(len= MAX_STRING_LEN), public, dimension(:),       allocatable  :: out_fwd_path
 character(len= MAX_STRING_LEN), public, dimension(:),       allocatable  :: in_dat_path
 character(len= MAX_STRING_LEN), public              :: source_fname
 character(len= MAX_STRING_LEN), public              :: station_fname

!################################################# FWAT PAR ######################################################################

 !!! verbose mode for with outputs for checking the FWI
 logical,                        public                                  :: VERBOSE_MODE

 integer, public   :: NSCOMP
 character(len= MAX_STRING_LEN), public, dimension(:),       allocatable :: SCOMPS
 integer, public   :: NRCOMP
 character(len= MAX_STRING_LEN), public, dimension(:),       allocatable :: RCOMPS
 integer, public   :: NUM_FILTER
 character(len=MAX_STRING_LEN),  public :: dat_coord
 character(len= MAX_STRING_LEN), public :: CH_CODE
 real(kind=CUSTOM_REAL), public, dimension(:),       allocatable     :: SHORT_P
 real(kind=CUSTOM_REAL), public, dimension(:),       allocatable     :: LONG_P
 real(kind=CUSTOM_REAL), public, dimension(:),       allocatable     :: GROUPVEL_MIN !BinHe added advised by Kai
 real(kind=CUSTOM_REAL), public, dimension(:),       allocatable     :: GROUPVEL_MAX  !BinHe added advised by Kai
 logical,                       public :: USE_NEAR_OFFSET !BinHe added
 real(kind=CUSTOM_REAL),        public :: TW_BEFORE, TW_AFTER ! time window before and after ttp
 logical,                       public :: ADJ_SRC_NORM
 logical,                       public :: SUPPRESS_EGF
 logical,                       public :: SAVE_OUTPUT_EACH_EVENT
 logical,                       public :: USE_SPH_SMOOTH ! If smoothing spherical mesh
 logical,                       public :: USE_PDE_SMOOTH ! If use pde-based smoothing mesh 
 real(kind=CUSTOM_REAL),        public :: SIGMA_H, SIGMA_V ! for smoothing
 character(len= MAX_STRING_LEN), public :: OPT_METHOD
 character(len= MAX_STRING_LEN), public :: PRECOND_TYPE
 integer,                       public :: NUM_STEP
 real(kind=CUSTOM_REAL),        public :: MAX_SLEN ! maximum step length for model update
 logical,                       public :: DO_LS
 real(kind=CUSTOM_REAL), public, dimension(:),   allocatable :: STEP_LENS ! for model update and line search
 logical,                       public :: CMT3D_INV=.false. ! for cmt3d_fwd
 character(len=3),              public :: CMT_PAR(9)
 integer,                       public :: FLEX_NWIN ! number of windows made on each measurement
!  real(kind=CUSTOM_REAL),        public :: MINDERR, RF_TSHIFT
!  integer, public                       :: NGAUSS
!  real(kind=CUSTOM_REAL), public, dimension(:), allocatable :: F0
!  integer, public                       :: MAXIT
 logical,                       public :: KERNEL_TAPER
 real(kind=CUSTOM_REAL),        public :: TAPER_H_SUPPRESS, TAPER_H_BUFFER, TAPER_V
end module fullwave_adjoint_tomo_par

!################################################# measure_adj ####################################################################

module ma_constants

  ! number of entries in window_chi output file
  integer, parameter :: N_MEASUREMENT = 5
  integer, parameter :: NCHI = 3*(N_MEASUREMENT-1) + 8

  ! constants
  double precision, parameter :: PI = 3.141592653589793d+00
  double precision, parameter :: TWOPI = 2.0 * PI
  complex(kind=8), parameter :: CCI = cmplx(0.,1.)
  double precision, parameter :: LARGE_VAL = 1.0d8

  ! FFT parameters
  integer, parameter :: LNPT = 16, NPT = 2**LNPT, NDIM = 80000
  double precision, parameter :: FORWARD_FFT = 1.0
  double precision, parameter :: REVERSE_FFT = -1.0

  ! phase correction control parameters, set this between (PI, 2PI),
  ! use a higher value for conservative phase wrapping
  double precision, parameter :: PHASE_STEP = 1.5 * PI

  ! filter parameters for xapiir bandpass subroutine (filter type is BP)
  ! (These should match the filter used in pre-processing.)
  double precision, parameter :: TRBDNDW = 0.3
  double precision, parameter :: APARM = 30.
  integer, parameter :: IORD = 4
  integer, parameter :: PASSES = 2

  ! takes waveform of first trace dat_dtw, without taking the difference waveform to the second trace syn_dtw
  ! this is useful to cissor out later reflections which appear in data (no synthetics needed)
  logical, parameter :: NO_WAVEFORM_DIFFERENCE = .false.

  ! constructs adjoint sources for a "ray density" kernel, where all misfits are equal to one
  logical, parameter :: DO_RAY_DENSITY_SOURCE = .false.

end module ma_constants

module ma_variables

  use ma_constants
!
! multi-taper measurements
!
! Ying Zhou: The fit between the recovered data and the data can be improved
! by either increasing the window width (HWIN above) or by decreasing NPI.
! In her experience, NPI = 2.5 is good for noisy data.
! For synthetic data, we can use a lower NPI.
! number of tapers should be fixed as twice NPI -- see Latex notes
!
! See write_par_file.pl and measure_adj.f90

  character(len=512) :: OUT_DIR

  double precision :: TLONG, TSHORT
  double precision :: WTR, NPI, DT_FAC, ERR_FAC, DT_MAX_SCALE, NCYCLE_IN_WINDOW
  !double precision :: BEFORE_QUALITY, AFTER_QUALITY, BEFORE_TSHIFT, AFTER_TSHIFT
  double precision :: TSHIFT_MIN, TSHIFT_MAX, DLNA_MIN, DLNA_MAX, CC_MIN
  double precision :: DT_SIGMA_MIN, DLNA_SIGMA_MIN

  integer :: ntaper, ipwr_t, ipwr_w, ERROR_TYPE
  integer :: imeas0, imeas, itaper, is_mtm0, is_mtm

  logical :: DISPLAY_DETAILS,OUTPUT_MEASUREMENT_FILES,RUN_BANDPASS,COMPUTE_ADJOINT_SOURCE,USE_PHYSICAL_DISPERSION

end module ma_variables


module ma_weighting

! module for weighting/normalizing measurements

  logical,parameter :: DO_WEIGHTING = .false.

  ! transverse, radial and vertical weights
  double precision :: weight_T, weight_R, weight_Z
  ! body waves: number of picks on vertical, radial and transverse component
  double precision :: num_P_SV_V,num_P_SV_R,num_SH_T
  ! surface waves: number of pick on vertical, radial and transverse
  double precision :: num_Rayleigh_V,num_Rayleigh_R,num_Love_T

  ! typical surface wave speed in km/s, to calculate surface wave arrival times
  ! Love waves faster than Rayleigh
  double precision, parameter :: surface_vel = 4.0

  ! wave type pick
  integer, parameter :: P_SV_V = 1
  integer, parameter :: P_SV_R = 2
  integer, parameter :: SH_T = 3
  integer, parameter :: Rayleigh_V = 4
  integer, parameter :: Rayleigh_R = 5
  integer, parameter :: Love_T = 6

end module ma_weighting

!################################################# flexwin ####################################################################

  module user_parameters
 
  use ma_constants, only: TRBDNDW,APARM,IORD,PASSES,NDIM
  !===================================================================
  ! filter parameters for xapiir subroutine (filter type is BP)
  !double precision, parameter :: TRBDNDW = 0.3
  !double precision, parameter :: APARM = 30.0
  !integer, parameter :: IORD = 5
  !integer, parameter :: PASSES = 2

  ! -------------------------------------------------------------
  ! array dimensions
  ! note that some integer arrays (iM,iL,iR) are NWINDOWS * NWINDOWS
  ! THESE SHOULD PROBABLY BE USER PARAMETERS, SINCE THEY WILL AFFECT
  ! THE SPEED OF THE PROGRAM (ALTHOUGH NOT THE OUTPUT).
  !integer, parameter :: NDIM = 10000
  integer, parameter :: NWINDOWS = 2500

  ! -------------------------------------------------------------
  ! miscellaneous - do not modify!
  ! -------------------------------------------------------------

  ! mathematical constants
  double precision, parameter :: PI = 3.1415926535897
  double precision, parameter :: E  = 2.7182818284590

  ! filter types
  integer, parameter :: HANNING = 1
  integer, parameter :: HAMMING = 2
  integer, parameter :: COSINE  = 3

  ! -------------------------------------------------------------

  end module user_parameters

