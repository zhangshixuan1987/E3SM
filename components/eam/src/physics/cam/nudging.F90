
module nudging
!=====================================================================
!
! See Patrick's original implementation in the section below
!
! Revised by Jian Sun, Kai Zhang, Shixuan Zhang, PNNL, 05/11/2020
!
! Revisions:
!   - Include the linear interpolation of nudge data to the model time step. 
!   - Add the "Nudge_Method" option in the namelist for different methods.
!   - Bug fix for intermittent nudged simulation:
!         * Reset the nudging tendency to zero if Update_Model is false.
!         * Nudge the model data to the same time slice of nudging data. 
!   - Update the nudging code for FV dycore.
!   - Add the "Nudge_Tau" option to control the relaxation 
!     timescale independently
!   - Add the "Nudge_Loc_PhysOut" option to calculate of nudging tendency
!     at the same location where the model state variables are written out.
!   - Add "Nudge_CurrentStep" option to linearly interpolate the nudging data 
!     to current or future time step. It only works appropriately when 
!     the nudging data starts with 00Z.
!   - Add "Nudge_File_Ntime" option to specify how many time slices per file
!   - To restore the functionality of the original nudging code, set the
!     new options with the values below:
!         Nudge_Tau         = -999
!         Nudge_Loc_PhysOut = .False.
!         Nudge_CurrentStep = .False.
!         Nudge_File_Ntime  = 1
!         Nudge_Method      = ‘Step’
!
! The revised nudging code has been evaluated by Sun et al. (2019).
! 
! Reference:
!
!    Sun, J., Zhang, K., Wan, H., Ma, P.‐L., Tang, Q., & Zhang, S. ( 2019).
!    Impact of nudging strategy on the climate representativeness and hindcast
!    skill of constrained EAMv1 simulations. Journal of Advances in Modeling
!    Earth Systems, 11, 3911– 3933. https://doi.org/10.1029/2019MS001831
!
!=====================================================================
!
! Purpose: Implement Nudging of the model state of U,V,T,Q, and/or PS
!          toward specified values from analyses.
!
! Author: Patrick Callaghan
!
! Description:
!    This module assumes that the user has {U,V,T,Q,PS} analyses which
!    have been preprocessed onto the current model grid and are stored
!    in individual files which are indexed with respect to year, month,
!    day, and second of the day. When the model is inbetween the given
!    begining and ending times, forcing is added to nudge the model toward
!    the appropriate analyses values. After the model passes the ending
!    analyses time, the forcing discontinues.
!
! Revisions:
!    01/14/13 - Modified to manage 'GAPS' in analyses data. For now the
!               approach is to coast through the gaps...  If a given
!               analyses file is missing, nudging is turned off for
!               that interval of time. Once an analyses file is found,
!               the Nudging is switched back on.
!    02/22/13 - Modified to add functionality for FV and EUL dynamical
!               cores.
!    03/03/13 - For ne120 runs, the automatic arrays used for reading in
!               U,V,T,Q,PS values were putting too much of a burden on the
!               stack memory. Until Parallel I/O is implemented, the impact
!               on the stack was reduced by using only one automatic array
!               to read in and scatter the data.
!    04/01/13 - Added Heaviside window function for localized nudging
!    04/10/13 - Modified call to physics_ptend_init() to accomodate the
!               new interface (in CESM1_2_BETA05).
!    05/06/13 - 'WRAP_NF' was modified from a generic interface so that
!               now it can only read in 1D arrays from netCDF files.
!               To eliminate errors from future meddling of this sort, all
!               refenences to the 'wrap_nf' module were removed and replaced
!               with direct nf90 calls.
!
! Input/Output Values:
!    Forcing contributions are available for history file output by
!    the names:    {'Nudge_U','Nudge_V','Nudge_T',and 'Nudge_Q'}
!
!    The nudging of the model toward the analyses data is controlled by
!    the 'nudging_nl' namelist in 'user_nl_cam'; whose variables control the
!               with direct nf90 calls.
!
! Input/Output Values:
!    Forcing contributions are available for history file output by
!    the names:    {'Nudge_U','Nudge_V','Nudge_T',and 'Nudge_Q'}
!
!    The nudging of the model toward the analyses data is controlled by
!    the 'nudging_nl' namelist in 'user_nl_cam'; whose variables control the
!    time interval over which nudging is applied, the strength of the nudging
!    tendencies, and its spatial distribution. The strength of the nudging is
!    specified as a fractional coeffcient between [0,1]. The spatial distribution 
!    is specified with a profile index:
!
!        (U,V,T,Q) Profiles:      0 == OFF      (No Nudging of this variable)
!        -------------------      1 == CONSTANT (Spatially Uniform Nudging)
!                                 2 == HEAVISIDE WINDOW FUNCTION
!
!        (PS) Profiles:           0 == OFF (Not Implemented)
!        -------------------      1 == N/A (Not Implemented)
!
!    The Heaviside window function is the product of separate horizonal and vertical 
!    windows that are controled via 14 parameters:
!        Nudge_Hwin_lat0:     Provide the horizontal center of the window in degrees. 
!        Nudge_Hwin_lon0:     The longitude must be in the range [0,360] and the 
!                             latitude should be [-90,+90].
!
!        Nudge_Hwin_latWidth: Specify the lat and lon widths of the window as positive 
!        Nudge_Hwin_lonWidth: values in degrees.Setting a width to a large value (e.g. 999) 
!                             renders the window a constant in that direction.
!
!        Nudge_Hwin_latDelta: Controls the sharpness of the window transition with a 
!        Nudge_Hwin_lonDelta: length in degrees. Small non-zero values yeild a step 
!                             function while a large value leads to a smoother transition.
!        Nudge_Hwin_Invert  : A logical flag used to invert the horizontal window function 
!                             to get its compliment.(e.g. to nudge outside a given window).
!
!        Nudge_Vwin_Lindex:   In the vertical, the window is specified in terms of model 
!        Nudge_Vwin_Ldelta:   level indcies. The High and Low transition levels should 
!        Nudge_Vwin_Hindex:   range from [0,(NCOL+1)]. The transition lengths are also 
!        Nudge_Vwin_Hdelta:   specified in terms of model indices. For a window function 
!                             constant in the vertical, the Low index should be set to 0,
!                             the High index should be set to (NCOL+1), and the transition 
!                             lengths should be set to 0.1
!        Nudge_Vwin_Invert  : A logical flag used to invert the vertical window function 
!                             to get its compliment.
!
!        Nudge_Hwin_lo:       For a given set of spatial parameters, the raw window 
!        Nudge_Hwin_hi:       function may not span the range [0,1], so those values are 
!        Nudge_Vwin_lo:       mapped to the range of values specified in by the user. 
!        Nudge_Vwin_hi:       The 'hi' values are mapped to the maximum of the raw window 
!                             function and 'lo' values are mapped to its minimum. 
!                             Typically the 'hi' values will be set equal to 1, and the 
!                             'lo' values set equal 0 or the desired window minimum. 
!                             Specifying the 'lo' value as 1 and the 'hi' value as 0 acts 
!                             to invert the window function. For a properly specified
!                             window its maximum should be equal to 1: MAX('lo','hi')==1
!
!        EXAMPLE: For a channel window function centered at the equator and independent 
!                 of the vertical (30 levels):
!                        Nudge_Hwin_lo = 0.               Nudge_Vwin_lo = 0.
!                        Nudge_Hwin_hi = 1.               Nudge_Vwin_hi = 1.
!                        Nudge_Hwin_lat0     = 0.         Nudge_Vwin_Lindex = 0.
!                        Nudge_Hwin_latWidth = 30.        Nudge_Vwin_Ldelta = 0.1
!                        Nudge_Hwin_latDelta = 5.0        Nudge_Vwin_Hindex = 31.
!                        Nudge_Hwin_lon0     = 180.       Nudge_Vwin_Hdelta = 0.1
!                        Nudge_Hwin_lonWidth = 999.
!                        Nudge_Hwin_lonDelta = 1.0
!
!                 If on the other hand one desired to apply nudging at the poles and
!                 not at the equator, the settings would be similar but with:
!                        Nudge_Hwin_lo = 1.
!                        Nudge_Hwin_hi = 0.
!
!    The following nudging options are added for more flexible controls of nudged simulations
!
!        Nudge_Method:      The method to perform Nudging analyses. It has to be
!                           one of the following three options:
!                           1. Step: the model meteorology is nudged toward the same (future) time
!                              slice of the constraining data within each nudging window (e.g., 6-hr)
!                           2. Linear: the model meteorology is nudged toward the state at the next
!                              model time step, which is linearly interpolated between two neighboring
!                              time slices of the constraining data
!                           3. IMT: the model meteorology is nudged toward the constraining data at
!                              the same model time step and nudging is only applied when the 
!                              constrainging  data is available
! 
!                           The first two options are continuous nudging while the
!                           third option is intermittent nudging. 
!
!                           Example: Nudge_Method = 'Linear'
!
!        Nudge_Tau:         User-defined relaxation time scale (unit: hour).
!                           If Nudge_Tau > 0, the relaxation time scale is determined by Nudge_Tau;
!                           If not, the relaxation time scale is determined by Nudge_Times_Per_Day 
!                           and nudging strength specified in the namelist (e.g. Nudge_Ucoef);
!
!                           Example: Nudge_Tau = -999
!
!        Nudge_Loc_PhysOut: If TRUE, change the location of calculation of nudging tendency to the
!                           same location where the nudging data is written out
!                           If FALSE, calculate the nudging tendency at the beginning of tphysbc (the
!                           same as the default model)
!
!                           Example: Nudge_Loc_PhysOut = TRUE
!
!        Nudge_CurrentStep: If TRUE, linearly interpolate nudging data to current model time step
!                           If FALSE, linearly interpolate nudging data to future model time step (the same as default model)
!
!                           Example: Nudge_CurrentStep = FALSE 
!
!        Nudge_File_Ntime:  Number of time slices per nudging data file
!                           The current nudging code only works for the nudging data file with one-day data.
!                           It does not work correctly when the nudging data file contains multiple-day data.
!                           Thus, Nudge_File_Ntime has to equal to 1 or Nudge_Times_Per_Day.
!                           If there are multiple time slices per file, the revised nudging code assumes 
!                           that the time slice in a file always starts from 00Z.
!
!                           Example: Nudge_File_Ntime = 1
!
!    &nudging_nl
!      Nudge_Model         - LOGICAL toggle to activate nudging.
!      Nudge_Path          - CHAR path to the analyses files.
!      Nudge_File_Template - CHAR Analyses filename with year, month, day, and second
!                                 values replaced by %y, %m, %d, and %s
!                                 respectively.         
!      Nudge_Times_Per_Day - INT Number of analyses files available per day.           
!      Model_Times_Per_Day - INT Number of times to update the model state
!      (used for nudging)              
!                                each day. The value is restricted to be longer
!                                than the              
!                                current model timestep and shorter than the
!                                analyses              
!                                timestep. As this number is increased, the
!                                nudging               
!                                force has the form of newtonian cooling.              
!      Nudge_Uprof         - INT index of profile structure to use for U.  [0,1,2]
!      Nudge_Vprof         - INT index of profile structure to use for V.  [0,1,2]
!      Nudge_Tprof         - INT index of profile structure to use for T.  [0,1,2]
!      Nudge_Qprof         - INT index of profile structure to use for Q.  [0,1,2]
!      Nudge_PSprof        - INT index of profile structure to use for PS. [0,N/A]
!      Nudge_Ucoef         - REAL fractional nudging coeffcient for U.
!                                    Utau=(Nudge_Ucoef/analyses_timestep)
!      Nudge_Vcoef         - REAL fractional nudging coeffcient for V.
!                                    Vtau=(Nudge_Vcoef/analyses_timestep)
!      Nudge_Tcoef         - REAL fractional nudging coeffcient for T.
!                                    Ttau=(Nudge_Tcoef/analyses_timestep)
!      Nudge_Qcoef         - REAL fractional nudging coeffcient for Q.
!                                    Qtau=(Nudge_Qcoef/analyses_timestep)
!      Nudge_PScoef        - REAL fractional nudging coeffcient for PS.
!                                    PStau=(Nudge_PScoef/analyses_timestep)
!      Nudge_Beg_Year      - INT nudging begining year.
!      Nudge_Beg_Month     - INT nudging begining month.
!      Nudge_Beg_Day       - INT nudging begining day.
!      Nudge_End_Year      - INT nudging ending year.
!      Nudge_End_Month     - INT nudging ending month.
!      Nudge_End_Day       - INT nudging ending day.
!      Nudge_Hwin_lo       - REAL value mapped to RAW horizontal window minimum. [0]
!      Nudge_Hwin_hi       - REAL value mapped to RAW horizontal window maximum. [1]
!      Nudge_Vwin_lo       - REAL value mapped to RAW vertical window minimum.   [0]
!      Nudge_Vwin_hi       - REAL value mapped to RAW vertical window maximum.   [1]
!      Nudge_Hwin_lat0     - REAL latitudinal center of window in degrees.
!      Nudge_Hwin_lon0     - REAL longitudinal center of window in degrees.
!      Nudge_Hwin_latWidth - REAL latitudinal width of window in degrees.
!      Nudge_Hwin_lonWidth - REAL longitudinal width of window in degrees.
!      Nudge_Hwin_latDelta - REAL latitudinal transition length of window in degrees.
!      Nudge_Hwin_lonDelta - REAL longitudinal transition length of window in degrees.
!      Nudge_Vwin_Lindex   - REAL LO model index of transition
!      Nudge_Vwin_Hindex   - REAL HI model index of transition
!      Nudge_Vwin_Ldelta   - REAL LO transition length
!      Nudge_Vwin_Hdelta   - REAL HI transition length
!      Nudge_Method        - CHAR method to perform Nudging analyses. [Step,Linear,IMT]
!      Nudge_Tau           - REAL relaxation time scale if Nudge_Tau > 0
!      Nudge_Loc_PhysOut   - LOGICAL change the location of calculation of nudging tendency 
!      Nudge_CurrentStep   - LOGICAL linearly interpolate nudging data to current
!                                    or future model time step
!      Nudge_File_Ntime    - INT  Number of time slices per nudging data file
!    /
!
! A typical namelist setting of nudged global simulation (UV-only, tau = 6h) 
! with this revised nudging code is provided below:
!
!    &nudging_nl
!      Nudge_Model         = .True.
!      Nudge_Path          = 'PATH_TO_DATA'
!      Nudge_File_Template = 'E3SM.cam.h1.%y-%m-%d-00000.nc'
!      Nudge_Times_Per_Day = 4            ! nudging input data frequency
!      Model_Times_Per_Day = 48           ! should not be larger than 48 if dtime = 1800s
!      Nudge_Uprof         = 1
!      Nudge_Ucoef         = 1.
!      Nudge_Vprof         = 1
!      Nudge_Vcoef         = 1.
!      Nudge_Tprof         = 0
!      Nudge_Tcoef         = 0.
!      Nudge_Qprof         = 0
!      Nudge_Qcoef         = 0.
!      Nudge_PSprof        = 0
!      Nudge_PScoef        = 0.
!      Nudge_Beg_Year      = 0000         ! begin year of nudged simulation
!      Nudge_Beg_Month     = 1            ! begin month of nudged simulation
!      Nudge_Beg_Day       = 1            ! begin day of nudged simulation
!      Nudge_End_Year      = 9999         ! end year of nudged simulation
!      Nudge_End_Month     = 1            ! end month of nudged simulation
!      Nudge_End_Day       = 1            ! end day of nudged simulation
!      Nudge_Method        = 'Linear'     ! use linear-interpolation nudging
!      Nudge_File_Ntime    = 4            ! if there are four time slices per nudging input data file 
!      Nudge_Loc_PhysOut   = .TRUE.
!    /
!
! If a user prefers one time slice per nudging input data file, simply change
! two settings above to the following settings:
!
!      Nudge_File_Template = 'E3SM.cam.h1.%y-%m-%d-%s.nc'
!      Nudge_File_Ntime    = 1
!
! More advanced example settings for the CONUS RRM nudged simulation (i.e., Table 3 of Tang et al. (2019), 
! using UV-only, tau=6h, four time slices per nudging input data file, nudging spatial window) 
! with this revised nudging code are:
!
!    &nudging_nl
!      Nudge_Model         = .true.
!      Nudge_Path          = 'PATH_TO_DATA'
!      Nudge_File_Template = ‘interim_se_%y%m%d00_%y%m%d18_TQUV.nc’
!      Nudge_Times_Per_Day = 4
!      Model_Times_Per_Day = 96 
!      Nudge_Uprof         = 2 
!      Nudge_Ucoef         = 1.00 
!      Nudge_Vprof         = 2 
!      Nudge_Vcoef         = 1.00 
!      Nudge_Tprof         = 0 
!      Nudge_Tcoef         = 0.00 
!      Nudge_Qprof         = 0 
!      Nudge_Qcoef         = 0.00 
!      Nudge_PSprof        = 0 
!      Nudge_PScoef        = 0.00 
!      Nudge_Beg_Year      = 2011 
!      Nudge_Beg_Month     = 1 
!      Nudge_Beg_Day       = 1
!      Nudge_End_Year      = 2011 
!      Nudge_End_Month     = 12 
!      Nudge_End_Day       = 31 
!      Nudge_Hwin_lo       = 1.0
!      Nudge_Hwin_hi       = 0.0 
!      Nudge_Hwin_lat0     = 38.0 
!      Nudge_Hwin_latWidth = 34.0
!      Nudge_Hwin_latDelta = 3.8 
!      Nudge_Hwin_lon0     = 254.0 
!      Nudge_Hwin_lonWidth = 44.0
!      Nudge_Hwin_lonDelta = 3.8 
!      Nudge_Vwin_lo       = 0.0 
!      Nudge_Vwin_hi       = 1.0 
!      Nudge_Vwin_Hindex   = 73.0 
!      Nudge_Vwin_Hdelta   = 0.1 
!      Nudge_Vwin_Lindex   = 0.0 
!      Nudge_Vwin_Ldelta   = 0.1
!      Nudge_Method        = 'Linear'
!      Nudge_File_Ntime    = 4
!      Nudge_Loc_PhysOut   = .TRUE.
!    /
!
! Reference:
!
!    Tang, Q., Klein, S. A., Xie, S., Lin, W., Golaz, J.-C., Roesler, E. L.,
!    Taylor, M. A., Rasch, P. J., Bader, D. C., Berg, L. K., Caldwell, P.,
!    Giangrande, S. E., Neale, R. B., Qian, Y., Riihimaki, L. D., Zender, C. S.,
!    Zhang, Y., and Zheng, X.: Regionally refined test bed in E3SM atmosphere
!    model version 1 (EAMv1) and applications for high-resolution modeling,
!    Geosci. Model Dev., 12, 2679–2706,
!    https://doi.org/10.5194/gmd-12-2679-2019, 2019.
!
!================
!  DIAG NOTE:
!================
!   The interface for reading and using analyses data is not complete for the FV
!   dynamical core. Wind values stored in the available data set are the values
!   on the staggered grid US,VS rather than U,V. To test the implementation of
!   the nudging for this case, the US,VS values were read in a loaded as if they
!   were U,V. The implementation of this hack is tagged with 'DIAG' where code
!   changed are needed to undo and fix what I have done.
!================
!
! TO DO:
! -----------
!    ** Currently the surface pressure is read in, but there is no forcing
!       meachnism implemented.
!    ** Analyses data is read in and then distributed to processing elements
!       via 'scatted_field_to_chunk' calls. The SE's want this to be changed
!       to parallel I/O calls.
!    ** Possibly implement time variation to nudging coeffcients, so that
!       rather than just bashing the model with a sledge hammer, the user has the
!       option to ramp up the nudging coefs over a startup time frame via a
!       heavyside step function.
!
!=====================================================================
  ! Useful modules
  !------------------
  use shr_kind_mod,only:r8=>SHR_KIND_R8,r4=>SHR_KIND_R4,cs=>SHR_KIND_CS,cl=>SHR_KIND_CL
  use time_manager,only:timemgr_time_ge,timemgr_time_inc,get_curr_date,dtime
  use phys_grid   ,only:scatter_field_to_chunk
  use cam_abortutils  ,only:endrun
  use spmd_utils  ,only:masterproc
  use cam_logfile ,only:iulog
  use shr_log_mod, only:errMsg => shr_log_errMsg
  use ioFileMod, only: getfil
  use cam_pio_utils, only: cam_pio_openfile
  use perf_mod
#ifdef SPMD
  use mpishorthand
#endif
  use torch_ftn

  ! Set all Global values and routines to private by default
  ! and then explicitly set their exposure.
  !----------------------------------------------------------
  implicit none
  private
  public:: Nudge_Model,Nudge_ON
  public:: Nudge_Allow_Missing_File
  public:: Nudge_Pdep_Weight_On
  public:: Nudge_Lin_Relax_On
  public:: Nudge_PS_Adjust_On
  public:: Nudge_PS_On
  public:: Nudge_Q_Adjust_On
  public:: Nudge_Land,Nudge_SRF_On
  public:: Nudge_SRF_Flux_On
  public:: Nudge_SRF_Q_On
  public:: Nudge_SRF_PSWgt_On
  public:: Nudge_SRF_Prec_On
  public:: Nudge_SRF_RadFlux_On
  public:: Nudge_SRF_State_On
  public:: nudging_readnl
  public:: nudging_init
  public:: nudging_timestep_init
  public:: nudging_timestep_tend
  public:: nudging_calc_tend
  public:: nudging_update_land_surface
  public:: nudging_update_srf_flux
  public:: Nudge_Loc_PhysOut
  private::se2latlon_interp_init
  private::latlon2se_interp_init
  private::dist_latlon_2pts 
  private::coord_ind_weight
  private::update_nudging_tend
  private::update_nudge_prof
  private::ps_nudging
  private::nudging_update_analyses_se
  private::nudging_update_srf_analyses_se
  private::nudging_update_analyses_eul
  private::nudging_update_analyses_fv
  private::nudging_set_PSprofile
  private::nudging_set_SRFprofile
  private::nudging_set_profile
  private::linear_interpolation
  private::linear_interpolation_2d
  private::read_and_scatter_se
  private::read_and_scatter_se_2d
  private::read_and_scatter_fv
  private::open_netcdf

  ! Machine Learning Bias Correction  
  public:: mltbc_nudge
  public:: mltbc_regional_on
  public:: mltbc_predict_option
  public:: mltbc_interp_test
  public:: mltbc_patch_nxy
  public:: mltbc_patch_nlon
  public:: mltbc_patch_nlat
  public:: mltbc_patch_dx
  public:: mltbc_patch_dy
  public:: mltbc_patch_biln
  public:: mltbc_timestep_init

  private::mltbc_advance
  private::mltbc_gather_data
  private::mltbc_gather_patch
  private::mltbc_reg_tendadv
  private::mltbc_glb_tendadv
  private::mltbc_don_statadv
  private::mltbc_don_encoder
  private::mltbc_don_tendadv
  private::mltbc_don_decoder

  ! Nudging Parameters
  !--------------------
  logical::         Nudge_Model       =.false.
  logical::         Nudge_ON          =.false.
  logical::         Nudge_File_Present=.false.
  logical::         Nudge_Initialized =.false.
  logical::         Nudge_Allow_Missing_File = .false. 
  logical::         Nudge_Pdep_Weight_On = .false.
  logical::         Nudge_SRF_File_Present=.false.
  logical::         Nudge_SRF_PSWgt_On   = .false.
  logical::         Nudge_SRF_Prec_On    = .false.
  logical::         Nudge_SRF_RadFlux_On = .false.
  logical::         Nudge_SRF_State_On   = .false. 
  logical::         Nudge_Lin_Relax_On   = .false.
  logical::         Nudge_PS_Adjust_On   = .false.
  logical::         Nudge_PS_On          = .false.
  logical::         Nudge_Q_Adjust_On    = .false.
  logical::         Nudge_Land           = .false.
  logical::         Nudge_SRF_On         = .false.
  logical::         Nudge_SRF_Flux_On    = .false. 
  logical::         Nudge_SRF_Q_On       = .false.
  logical::         Nudge_Hwin_Invert    = .false.
  logical::         Nudge_Vwin_Invert    = .false.
  character(len=cl) Nudge_Path
  character(len=cl) Nudge_File,Nudge_File_Template
  character(len=cl) Nudge_SRF_File,Nudge_SRF_File_Template
  integer           Nudge_Times_Per_Day
  integer           Model_Times_Per_Day
  real(r8)          Nudge_Ucoef,Nudge_Vcoef
  integer           Nudge_Uprof,Nudge_Vprof
  real(r8)          Nudge_Qcoef,Nudge_Tcoef
  integer           Nudge_Qprof,Nudge_Tprof
  real(r8)          Nudge_PScoef
  integer           Nudge_PSprof
  real(r8)          Nudge_SRFcoef
  integer           Nudge_SRFprof
  integer           Nudge_Beg_Year ,Nudge_Beg_Month
  integer           Nudge_Beg_Day  ,Nudge_Beg_Sec
  integer           Nudge_End_Year ,Nudge_End_Month
  integer           Nudge_End_Day  ,Nudge_End_Sec
  integer           Nudge_Curr_Year,Nudge_Curr_Month
  integer           Nudge_Curr_Day ,Nudge_Curr_Sec
  integer           Nudge_Next_Year,Nudge_Next_Month
  integer           Nudge_Next_Day ,Nudge_Next_Sec
  integer           Nudge_Step
  integer           Model_Curr_Year,Model_Curr_Month
  integer           Model_Curr_Day ,Model_Curr_Sec
  integer           Model_Next_Year,Model_Next_Month
  integer           Model_Next_Day ,Model_Next_Sec
  integer           Model_Step
  real(r8)          Nudge_Hwin_lo
  real(r8)          Nudge_Hwin_hi
  real(r8)          Nudge_Hwin_lat0
  real(r8)          Nudge_Hwin_latWidth
  real(r8)          Nudge_Hwin_latDelta
  real(r8)          Nudge_Hwin_lon0
  real(r8)          Nudge_Hwin_lonWidth
  real(r8)          Nudge_Hwin_lonDelta
  real(r8)          Nudge_Vwin_lo
  real(r8)          Nudge_Vwin_hi
  real(r8)          Nudge_Vwin_Hindex
  real(r8)          Nudge_Vwin_Hdelta
  real(r8)          Nudge_Vwin_Lindex
  real(r8)          Nudge_Vwin_Ldelta
  real(r8)          Nudge_Hwin_latWidthH
  real(r8)          Nudge_Hwin_lonWidthH
  real(r8)          Nudge_Hwin_max
  real(r8)          Nudge_Hwin_min

  ! Nudging State Arrays
  !-----------------------
  integer Nudge_nlon,Nudge_nlat,Nudge_ncol,Nudge_nlev
!DIAG
  integer Nudge_slat
!DIAG
  real(r8),allocatable::Target_U(:,:,:)     !(pcols,pver,begchunk:endchunk)
  real(r8),allocatable::Target_V(:,:,:)     !(pcols,pver,begchunk:endchunk)
  real(r8),allocatable::Target_T(:,:,:)     !(pcols,pver,begchunk:endchunk)
  real(r8),allocatable::Target_Q(:,:,:)     !(pcols,pver,begchunk:endchunk)
  real(r8),allocatable::Target_PS(:,:)      !(pcols,begchunk:endchunk)
  real(r8),allocatable::Target_PHIS(:,:)    !(pcols,begchunk:endchunk)
  real(r8),allocatable::Model_U(:,:,:)      !(pcols,pver,begchunk:endchunk)
  real(r8),allocatable::Model_V(:,:,:)      !(pcols,pver,begchunk:endchunk)
  real(r8),allocatable::Model_T(:,:,:)      !(pcols,pver,begchunk:endchunk)
  real(r8),allocatable::Model_Q(:,:,:)      !(pcols,pver,begchunk:endchunk)
  real(r8),allocatable::Model_PS(:,:)       !(pcols,begchunk:endchunk)
  real(r8),allocatable::Model_PHIS(:,:)     !(pcols,begchunk:endchunk)
  real(r8),allocatable::Nudge_Utau(:,:,:)   !(pcols,pver,begchunk:endchunk)
  real(r8),allocatable::Nudge_Vtau(:,:,:)   !(pcols,pver,begchunk:endchunk)
  real(r8),allocatable::Nudge_Ttau(:,:,:)   !(pcols,pver,begchunk:endchunk)
  real(r8),allocatable::Nudge_Qtau(:,:,:)   !(pcols,pver,begchunk:endchunk)
  real(r8),allocatable::Nudge_PStau(:,:)    !(pcols,begchunk:endchunk)
  real(r8),allocatable::Nudge_SRFtau(:,:)   !(pcols,begchunk:endchunk)
  real(r8),allocatable::Nudge_Ustep(:,:,:)  !(pcols,pver,begchunk:endchunk)
  real(r8),allocatable::Nudge_Vstep(:,:,:)  !(pcols,pver,begchunk:endchunk)
  real(r8),allocatable::Nudge_Tstep(:,:,:)  !(pcols,pver,begchunk:endchunk)
  real(r8),allocatable::Nudge_Qstep(:,:,:)  !(pcols,pver,begchunk:endchunk)
  real(r8),allocatable::Nudge_PSstep(:,:)   !(pcols,begchunk:endchunk)
  logical            :: l_Update_Model, & 
                        l_Update_Nudge, &
                        l_After_Beg,    &
                        l_Before_End
  real(r8),parameter :: sec_per_hour = 3600._r8
  real(r8),parameter :: fillvalue = 1.e+20_r8               ! fill value for netcdf fields
  integer,parameter  :: imissing = -99999                   ! special value for integer missing data 
  real(r8),parameter :: rmissing = -1.e36_r8                ! special value for missing data 

  character(len=10)  :: Nudge_Method                        ! nudge method 
  logical            :: Nudge_Loc_PhysOut                   ! whether nudging tendency is calculated at the same 
                                                            ! location where the model state variables are written out
  real(r8)           :: Nudge_Tau                           ! nudge relaxation timescale
  real(r8)           :: Nudge_SRF_Tau                       ! nudge relaxation timescale for land surface 
  logical            :: Nudge_CurrentStep                   ! .true. if linearly interpolated to current model time step
  integer            :: Nudge_File_Ntime                    ! number of time slices per nudging data file 
  integer            :: Nudge_SRF_File_Ntime                ! number of time slices per nudging data file for surface 
  logical :: first_file                                     ! the flag for first nudge data
  logical :: first_srf_file                                 ! the flag for first surface nudge data 
  real(r8), allocatable, dimension(:,:,:,:) :: INTP_U       ! (pcols,pver,begchunk:endchunk,:)
  real(r8), allocatable, dimension(:,:,:,:) :: INTP_V       ! (pcols,pver,begchunk:endchunk,:)
  real(r8), allocatable, dimension(:,:,:,:) :: INTP_T       ! (pcols,pver,begchunk:endchunk,:)
  real(r8), allocatable, dimension(:,:,:,:) :: INTP_Q       ! (pcols,pver,begchunk:endchunk,:)
  real(r8), allocatable, dimension(:,:,:)   :: INTP_PS      ! (pcols,begchunk:endchunk,:)
  real(r8), allocatable, dimension(:,:,:)   :: INTP_PHIS    ! (pcols,begchunk:endchunk,:)

!Variables for machine learning bias correction 
  real(r8), allocatable :: Model_rlat(:)            !(Nudge_ncol)
  real(r8), allocatable :: Model_rlon(:)            !(Nudge_ncol)
  real(r8), allocatable :: Model_wgth(:)            !(Nudge_ncol)
  real(r8), allocatable :: Model_Var(:,:)           !(Nudge_ncol,Nudge_nlev)

  real(r8), allocatable :: mltbc_patch_lon(:)      !(mltbc_patch_nxy)
  real(r8), allocatable :: mltbc_patch_lat(:)      !(mltbc_patch_nxy)
  integer,  allocatable :: mltbc_se2latlon_ind(:,:) !(mltbc_patch_nxy,5)
  integer,  allocatable :: mltbc_latlon2se_ind(:,:) !(Nudge_ncol,4)
  real(r8), allocatable :: mltbc_se2latlon_wgt(:,:) !(mltbc_patch_nxy,5)
  real(r8), allocatable :: mltbc_latlon2se_wgt(:,:) !(Nudge_ncol,4)

!Variables for surface nudging 
  real(r8), allocatable :: Target_U10(:,:)     !(pcols,begchunk:endchunk)
  real(r8), allocatable :: Target_V10(:,:)     !(pcols,begchunk:endchunk)
  real(r8), allocatable :: Target_T2(:,:)      !(pcols,begchunk:endchunk)
  real(r8), allocatable :: Target_TD2(:,:)     !(pcols,begchunk:endchunk)
  real(r8), allocatable :: Target_TS(:,:)      !(pcols,begchunk:endchunk)
  real(r8), allocatable :: Target_Q2(:,:)      !(pcols,begchunk:endchunk)
  real(r8), allocatable :: Target_PRECC(:,:)   !(pcols,begchunk:endchunk)
  real(r8), allocatable :: Target_PRECL(:,:)   !(pcols,begchunk:endchunk)
  real(r8), allocatable :: Target_PRECSC(:,:)  !(pcols,begchunk:endchunk)
  real(r8), allocatable :: Target_PRECSL(:,:)  !(pcols,begchunk:endchunk)
  real(r8), allocatable :: Target_EVAP(:,:)    !(pcols,begchunk:endchunk)
  real(r8), allocatable :: Target_SHFLX(:,:)   !(pcols,begchunk:endchunk)
  real(r8), allocatable :: Target_LHFLX(:,:)   !(pcols,begchunk:endchunk)
  real(r8), allocatable :: Target_NETSW(:,:)   !(pcols,begchunk:endchunk)
  real(r8), allocatable :: Target_FLWDS(:,:)   !(pcols,begchunk:endchunk)
  real(r8), allocatable :: Target_FSDS(:,:)    !(pcols,begchunk:endchunk)
  real(r8), allocatable :: Target_FSDSD(:,:)   !(pcols,begchunk:endchunk)
  real(r8), allocatable :: Target_FSDSUV(:,:)  !(pcols,begchunk:endchunk)

  real(r8), allocatable, dimension(:,:,:) :: INTP_U10   ! (pcols,begchunk:endchunk,:)
  real(r8), allocatable, dimension(:,:,:) :: INTP_V10   ! (pcols,begchunk:endchunk,:)
  real(r8), allocatable, dimension(:,:,:) :: INTP_T2    ! (pcols,begchunk:endchunk,:)
  real(r8), allocatable, dimension(:,:,:) :: INTP_TD2   ! (pcols,begchunk:endchunk,:)
  real(r8), allocatable, dimension(:,:,:) :: INTP_TS    ! (pcols,begchunk:endchunk,:)
  real(r8), allocatable, dimension(:,:,:) :: INTP_PRECC ! (pcols,begchunk:endchunk,:)
  real(r8), allocatable, dimension(:,:,:) :: INTP_PRECL ! (pcols,begchunk:endchunk,:)
  real(r8), allocatable, dimension(:,:,:) :: INTP_PRECSC! (pcols,begchunk:endchunk,:)
  real(r8), allocatable, dimension(:,:,:) :: INTP_PRECSL! (pcols,begchunk:endchunk,:)
  real(r8), allocatable, dimension(:,:,:) :: INTP_EVAP  ! (pcols,begchunk:endchunk,:)
  real(r8), allocatable, dimension(:,:,:) :: INTP_LHFLX ! (pcols,begchunk:endchunk,:)
  real(r8), allocatable, dimension(:,:,:) :: INTP_SHFLX ! (pcols,begchunk:endchunk,:)
  real(r8), allocatable, dimension(:,:,:) :: INTP_FSNS  ! (pcols,begchunk:endchunk,:)
  real(r8), allocatable, dimension(:,:,:) :: INTP_FLDS  ! (pcols,begchunk:endchunk,:)
  real(r8), allocatable, dimension(:,:,:) :: INTP_FSDS  ! (pcols,begchunk:endchunk,:)
  real(r8), allocatable, dimension(:,:,:) :: INTP_FSDSD ! (pcols,begchunk:endchunk,:)
  real(r8), allocatable, dimension(:,:,:) :: INTP_FSDSUV! (pcols,begchunk:endchunk,:)

  integer  :: Nudge_PS_OPT    ! option for the surface pressure nudging 
  integer  :: Nudge_UV_OPT    ! option for the zonal wind nudging 
  integer  :: Nudge_T_OPT     ! option for the temperature nudging 
  integer  :: Nudge_Q_OPT     ! option for the humidity nudging 

  integer  :: Nudge_NO_PBL_UV ! option for excluding UV nudging within PBL 
  integer  :: Nudge_NO_PBL_T  ! option for excluding T nudging within PBL 
  integer  :: Nudge_NO_PBL_Q  ! option for excluding Q nudging within PBL 

  !For surface flux nudging 
  integer  :: prec_dp_idx  = 0
  integer  :: snow_dp_idx  = 0
  integer  :: prec_sh_idx  = 0
  integer  :: snow_sh_idx  = 0
  integer  :: prec_sed_idx = 0
  integer  :: snow_sed_idx = 0
  integer  :: prec_pcw_idx = 0
  integer  :: snow_pcw_idx = 0
  integer  :: vmag_gust_idx= 0

  !Parameters determined with experiment 
  !From p_relax upwards, nudging is reduced linearly 
  real(r8), parameter :: p_uv_relax = 30.E2_r8  ! p_relax for u/v wind 
  real(r8), parameter :: p_T_relax  = 10.E2_r8  ! p_relax for temperature  
  real(r8), parameter :: p_q_relax  = 100.E2_r8 ! p_relax for humidity 
  real(r8), parameter :: p_norelax  = 1.0_r8    ! from p_norelax upwards, no nudging
  real(r8), parameter :: z_min      = 150._r8   ! height levels below which nudging is turned off  

  !Parameters for machine learning bias correction 
  logical  :: mltbc_nudge       = .false.
  logical  :: mltbc_regional_on = .false.
  logical  :: mltbc_patch_biln  = .false.
  logical  :: mltbc_interp_test = .false.

  integer  :: mltbc_predict_option
  integer  :: mltbc_patch_nlon
  integer  :: mltbc_patch_nlat
  integer  :: mltbc_patch_nxy
  real(r8) :: mltbc_patch_dx
  real(r8) :: mltbc_patch_dy

  !logical variables
  logical  :: l_mltbc_encoder
  logical  :: l_mltbc_decoder
  logical  :: l_mltbc_predictor

  character(len=cl) :: mltbc_model_path
  character(len=cl) :: file_encoder     ! Machine Learning Decoder pt file
  character(len=cl) :: file_decoder     ! Machine Learning Encoder pt file
  character(len=cl) :: file_predictor   ! Machine Learning Model   pt file

contains
  !================================================================
  subroutine nudging_readnl(nlfile)
   !
   ! NUDGING_READNL: Initialize default values controlling the Nudging
   !                 process. Then read namelist values to override
   !                 them.
   !===============================================================
   use ppgrid        ,only: pver
   use namelist_utils,only:find_group_name
   use units         ,only:getunit,freeunit
   !
   ! Arguments
   !-------------
   character(len=*),intent(in)::nlfile
   !
   ! Local Values
   !---------------
   integer ierr,unitn

   namelist /nudging_nl/ Nudge_Model,Nudge_Path,                       &
                         Nudge_Allow_Missing_File,                     & 
                         Nudge_File_Template,Nudge_Times_Per_Day,      &
                         Model_Times_Per_Day,                          &
                         Nudge_Ucoef ,Nudge_Uprof,                     &
                         Nudge_Vcoef ,Nudge_Vprof,                     &
                         Nudge_Qcoef ,Nudge_Qprof,                     &
                         Nudge_Tcoef ,Nudge_Tprof,                     &
                         Nudge_PScoef,Nudge_PSprof,                    &
                         Nudge_SRFcoef,Nudge_SRFprof,                  &
                         Nudge_Beg_Year,Nudge_Beg_Month,Nudge_Beg_Day, &
                         Nudge_End_Year,Nudge_End_Month,Nudge_End_Day, &
                         Nudge_Hwin_lo,Nudge_Hwin_hi,                  &
                         Nudge_Vwin_lo,Nudge_Vwin_hi,                  &
                         Nudge_Hwin_lat0,Nudge_Hwin_lon0,              &
                         Nudge_Hwin_latWidth,Nudge_Hwin_lonWidth,      &
                         Nudge_Hwin_latDelta,Nudge_Hwin_lonDelta,      &
                         Nudge_Vwin_Lindex,Nudge_Vwin_Hindex,          &
                         Nudge_Vwin_Ldelta,Nudge_Vwin_Hdelta,          &
                         Nudge_Hwin_Invert,Nudge_Vwin_Invert,          &
                         Nudge_Method, Nudge_Tau, Nudge_Loc_PhysOut,   &
                         Nudge_CurrentStep, Nudge_File_Ntime,          & 
                         Nudge_NO_PBL_UV, Nudge_NO_PBL_T,              &
                         Nudge_NO_PBL_Q, Nudge_Pdep_Weight_On,         &
                         Nudge_Lin_Relax_On, Nudge_PS_OPT,             & 
                         Nudge_UV_OPT, Nudge_T_OPT, Nudge_Q_OPT,       &
                         Nudge_PS_Adjust_On, Nudge_Q_Adjust_On,        & 
                         Nudge_SRF_File_Template,Nudge_SRF_File_Ntime, &
                         Nudge_Land, Nudge_SRF_Flux_On,                & 
                         Nudge_SRF_PSWgt_On, Nudge_SRF_Prec_On,        & 
                         Nudge_SRF_RadFlux_On, Nudge_SRF_State_On,     & 
                         Nudge_SRF_Q_On, Nudge_SRF_Tau,                &
                         mltbc_model_path, mltbc_regional_on,          & 
                         mltbc_predict_option, mltbc_interp_test,      & 
                         mltbc_patch_biln

  
   ! Nudging is NOT initialized yet, For now
   ! Nudging will always begin/end at midnight.
   !--------------------------------------------
   Nudge_Initialized =.false.
   Nudge_ON          =.false.
   Nudge_PS_On       =.false.
   Nudge_SRF_On      =.false.
   Nudge_SRF_Flux_On =.false.
   Nudge_SRF_Q_On    =.false. 
   Nudge_File_Present=.false.
   Nudge_SRF_File_Present=.false.
   Nudge_Beg_Sec=0
   Nudge_End_Sec=0

   ! Set Default Namelist values
   !-----------------------------
   Nudge_Model        =.false.
   Nudge_Allow_Missing_File = .false.
   Nudge_Pdep_Weight_On = .false.
   Nudge_Lin_Relax_On   = .false.
   Nudge_PS_Adjust_On   = .false.
   Nudge_Q_Adjust_On    = .false.
   Nudge_Land           = .false.
   Nudge_SRF_PSWgt_On   = .false. 
   Nudge_SRF_Prec_On    = .false.
   Nudge_SRF_RadFlux_On = .false.
   Nudge_SRF_State_On   = .false.
   Nudge_Path         ='./Data/YOTC_ne30np4_001/'
   Nudge_File_Template='YOTC_ne30np4_L30.cam2.i.%y-%m-%d-%s.nc'
   Nudge_Times_Per_Day=4
   Model_Times_Per_Day=4
   Nudge_Ucoef   =0._r8
   Nudge_Vcoef   =0._r8
   Nudge_Qcoef   =0._r8
   Nudge_Tcoef   =0._r8
   Nudge_PScoef  =0._r8
   Nudge_SRFcoef =0._r8
   Nudge_Uprof   =0
   Nudge_Vprof   =0
   Nudge_Qprof   =0
   Nudge_Tprof   =0
   Nudge_PSprof  =0
   Nudge_SRFprof =0
   Nudge_Beg_Year =2008
   Nudge_Beg_Month=5
   Nudge_Beg_Day  =1
   Nudge_End_Year =2008
   Nudge_End_Month=9
   Nudge_End_Day  =1
   Nudge_Hwin_lo      =0._r8
   Nudge_Hwin_hi      =1.0_r8
   Nudge_Hwin_lat0    =0._r8
   Nudge_Hwin_latWidth=9999._r8
   Nudge_Hwin_latDelta=1.0_r8
   Nudge_Hwin_lon0    =180._r8
   Nudge_Hwin_lonWidth=9999._r8
   Nudge_Hwin_lonDelta=1.0_r8
   Nudge_Hwin_Invert  = .false.
   Nudge_Vwin_lo      =0._r8
   Nudge_Vwin_hi      =1.0_r8
   Nudge_Vwin_Hindex  =float(pver+1)
   Nudge_Vwin_Hdelta  =0.1_r8
   Nudge_Vwin_Lindex  =0._r8
   Nudge_Vwin_Ldelta  =0.1_r8
   Nudge_Vwin_Invert  = .false.
   Nudge_Method       = 'Linear'
   Nudge_Loc_PhysOut  = .false.
   Nudge_Tau          = -999._r8
   Nudge_CurrentStep  = .false.
   Nudge_File_Ntime   = 0

   Nudge_SRF_File_Template = 'YOTC_sfc_ne30np4_L30.cam2.i.%y-%m-%d-%s.nc'
   Nudge_SRF_File_Ntime    = 0 
   Nudge_SRF_Tau           = -999._r8 

   Nudge_UV_OPT       = 0 
   Nudge_T_OPT        = 0
   Nudge_Q_OPT        = 0

   Nudge_NO_PBL_UV    = 0 
   Nudge_NO_PBL_T     = 0 
   Nudge_NO_PBL_Q     = 0
   
   ! Set Default values for machine learing 
   !-----------------------------
   mltbc_patch_biln     = .false.
   mltbc_interp_test    = .false.
   mltbc_nudge          = .false.
   mltbc_regional_on    = .false.
   mltbc_predict_option = 0
   mltbc_patch_nlon     = 1
   mltbc_patch_nlat     = 1
   mltbc_patch_nxy      = 1
   mltbc_patch_dx       = 1.0_r8
   mltbc_patch_dy       = 1.0_r8
   mltbc_model_path     = './mltbc_PT_File/'

   ! Read in namelist values
   !------------------------
   if(masterproc) then
     unitn = getunit()
     open(unitn,file=trim(nlfile),status='old')
     call find_group_name(unitn,'nudging_nl',status=ierr)
     if(ierr.eq.0) then
       read(unitn,nudging_nl,iostat=ierr)
       if(ierr.ne.0) then
         call endrun('nudging_readnl:: ERROR reading namelist')
       endif
     endif
     close(unitn)
     call freeunit(unitn)
   endif

   if (trim(Nudge_Method).eq."MLTBC") then 
     mltbc_nudge = .true. 
   else
     mltbc_nudge = .false.
   end if 

   ! Set hi/lo values according to the given '_Invert' parameters
   !--------------------------------------------------------------
   if(Nudge_Hwin_Invert) then
     Nudge_Hwin_lo = 1.0_r8
     Nudge_Hwin_hi = 0.0_r8
   else
     Nudge_Hwin_lo = 0.0_r8
     Nudge_Hwin_hi = 1.0_r8
   endif

   if(Nudge_Vwin_Invert) then
     Nudge_Vwin_lo = 1.0_r8
     Nudge_Vwin_hi = 0.0_r8
   else
     Nudge_Vwin_lo = 0.0_r8
     Nudge_Vwin_hi = 1.0_r8
   endif

   ! Check for valid namelist values
   !----------------------------------
   if((max(Nudge_Hwin_lo,Nudge_Hwin_hi).ne.1.0).or. &
      (max(Nudge_Vwin_lo,Nudge_Vwin_hi).ne.1.0)   ) then
     write(iulog,*) 'NUDGING: The window function must have a maximum value of 1'
     write(iulog,*) 'NUDGING:  Nudge_Hwin_lo=',Nudge_Hwin_lo
     write(iulog,*) 'NUDGING:  Nudge_Hwin_hi=',Nudge_Hwin_hi
     write(iulog,*) 'NUDGING:  Nudge_Vwin_lo=',Nudge_Vwin_lo
     write(iulog,*) 'NUDGING:  Nudge_Vwin_hi=',Nudge_Vwin_hi
     call endrun('nudging_readnl:: ERROR in namelist')
   endif

   if((Nudge_Hwin_lat0.lt.-90.).or.(Nudge_Hwin_lat0.gt.+90.)) then
     write(iulog,*) 'NUDGING: Window lat0 must be in [-90,+90]'
     write(iulog,*) 'NUDGING:  Nudge_Hwin_lat0=',Nudge_Hwin_lat0
     call endrun('nudging_readnl:: ERROR in namelist')
   endif

   if((Nudge_Hwin_lon0.lt.0.).or.(Nudge_Hwin_lon0.ge.360.)) then
     write(iulog,*) 'NUDGING: Window lon0 must be in [0,+360)'
     write(iulog,*) 'NUDGING:  Nudge_Hwin_lon0=',Nudge_Hwin_lon0
     call endrun('nudging_readnl:: ERROR in namelist')
   endif

   if((Nudge_Vwin_Lindex.gt.Nudge_Vwin_Hindex)                         .or. &
      (Nudge_Vwin_Hindex.gt.float(pver+1)).or.(Nudge_Vwin_Hindex.lt.0.).or. &
      (Nudge_Vwin_Lindex.gt.float(pver+1)).or.(Nudge_Vwin_Lindex.lt.0.)   ) then
     write(iulog,*) 'NUDGING: Window Lindex must be in [0,pver+1]'
     write(iulog,*) 'NUDGING: Window Hindex must be in [0,pver+1]'
     write(iulog,*) 'NUDGING: Lindex must be LE than Hindex'
     write(iulog,*) 'NUDGING:  Nudge_Vwin_Lindex=',Nudge_Vwin_Lindex
     write(iulog,*) 'NUDGING:  Nudge_Vwin_Hindex=',Nudge_Vwin_Hindex
     call endrun('nudging_readnl:: ERROR in namelist')
   endif

   if((Nudge_Hwin_latDelta.le.0.).or.(Nudge_Hwin_lonDelta.le.0.).or. &
      (Nudge_Vwin_Hdelta  .le.0.).or.(Nudge_Vwin_Ldelta  .le.0.)    ) then
     write(iulog,*) 'NUDGING: Window Deltas must be positive'
     write(iulog,*) 'NUDGING:  Nudge_Hwin_latDelta=',Nudge_Hwin_latDelta
     write(iulog,*) 'NUDGING:  Nudge_Hwin_lonDelta=',Nudge_Hwin_lonDelta
     write(iulog,*) 'NUDGING:  Nudge_Vwin_Hdelta=',Nudge_Vwin_Hdelta
     write(iulog,*) 'NUDGING:  Nudge_Vwin_Ldelta=',Nudge_Vwin_Ldelta
     call endrun('nudging_readnl:: ERROR in namelist')

   endif

   if((Nudge_Hwin_latWidth.le.0.).or.(Nudge_Hwin_lonWidth.le.0.)) then
     write(iulog,*) 'NUDGING: Window widths must be positive'
     write(iulog,*) 'NUDGING:  Nudge_Hwin_latWidth=',Nudge_Hwin_latWidth
     write(iulog,*) 'NUDGING:  Nudge_Hwin_lonWidth=',Nudge_Hwin_lonWidth
     call endrun('nudging_readnl:: ERROR in namelist')
   endif

   !if ( Nudge_PSprof .ne. 0 ) then
   !  write(iulog,*) 'NUDGING: PS nudging scheme is not implemented yet'
   !  write(iulog,*) 'NUDGING: Nudge_PSprof must be set to zero at this moment'
   !  call endrun('nudging_readnl:: ERROR in namelist')
   !end if

   ! Broadcast namelist variables
   !------------------------------
#ifdef SPMD
   call mpibcast(Nudge_Path              ,len(Nudge_Path)         ,mpichar,0,mpicom)
   call mpibcast(Nudge_File_Template     ,len(Nudge_File_Template),mpichar,0,mpicom)
   call mpibcast(Nudge_Model             , 1, mpilog, 0, mpicom)
   call mpibcast(Nudge_Initialized       , 1, mpilog, 0, mpicom)
   call mpibcast(Nudge_ON                , 1, mpilog, 0, mpicom)
   call mpibcast(Nudge_File_Present      , 1, mpilog, 0, mpicom)
   call mpibcast(Nudge_Allow_Missing_File, 1, mpilog, 0, mpicom)
   call mpibcast(Nudge_Pdep_Weight_On    , 1, mpilog, 0, mpicom)
   call mpibcast(Nudge_Lin_Relax_On      , 1, mpilog, 0, mpicom)
   call mpibcast(Nudge_PS_Adjust_On      , 1, mpilog, 0, mpicom)
   call mpibcast(Nudge_PS_On             , 1, mpilog, 0, mpicom)
   call mpibcast(Nudge_Q_Adjust_On       , 1, mpilog, 0, mpicom)
   call mpibcast(Nudge_Times_Per_Day     , 1, mpiint, 0, mpicom)
   call mpibcast(Model_Times_Per_Day     , 1, mpiint, 0, mpicom)
   call mpibcast(Nudge_Ucoef    , 1, mpir8 , 0, mpicom)
   call mpibcast(Nudge_Vcoef    , 1, mpir8 , 0, mpicom)
   call mpibcast(Nudge_Tcoef    , 1, mpir8 , 0, mpicom)
   call mpibcast(Nudge_Qcoef    , 1, mpir8 , 0, mpicom)
   call mpibcast(Nudge_PScoef   , 1, mpir8 , 0, mpicom)
   call mpibcast(Nudge_Uprof    , 1, mpiint, 0, mpicom)
   call mpibcast(Nudge_Vprof    , 1, mpiint, 0, mpicom)
   call mpibcast(Nudge_Tprof    , 1, mpiint, 0, mpicom)
   call mpibcast(Nudge_Qprof    , 1, mpiint, 0, mpicom)
   call mpibcast(Nudge_PSprof   , 1, mpiint, 0, mpicom)
   call mpibcast(Nudge_Beg_Year , 1, mpiint, 0, mpicom)
   call mpibcast(Nudge_Beg_Month, 1, mpiint, 0, mpicom)
   call mpibcast(Nudge_Beg_Day  , 1, mpiint, 0, mpicom)
   call mpibcast(Nudge_Beg_Sec  , 1, mpiint, 0, mpicom)
   call mpibcast(Nudge_End_Year , 1, mpiint, 0, mpicom)
   call mpibcast(Nudge_End_Month, 1, mpiint, 0, mpicom)
   call mpibcast(Nudge_End_Day  , 1, mpiint, 0, mpicom)
   call mpibcast(Nudge_End_Sec  , 1, mpiint, 0, mpicom)
   call mpibcast(Nudge_Hwin_lo      , 1, mpir8 , 0, mpicom)
   call mpibcast(Nudge_Hwin_hi      , 1, mpir8 , 0, mpicom)
   call mpibcast(Nudge_Hwin_lat0    , 1, mpir8 , 0, mpicom)
   call mpibcast(Nudge_Hwin_latWidth, 1, mpir8 , 0, mpicom)
   call mpibcast(Nudge_Hwin_latDelta, 1, mpir8 , 0, mpicom)
   call mpibcast(Nudge_Hwin_lon0    , 1, mpir8 , 0, mpicom)
   call mpibcast(Nudge_Hwin_lonWidth, 1, mpir8 , 0, mpicom)
   call mpibcast(Nudge_Hwin_lonDelta, 1, mpir8 , 0, mpicom)
   call mpibcast(Nudge_Hwin_Invert,   1, mpilog, 0, mpicom)
   call mpibcast(Nudge_Vwin_lo      , 1, mpir8 , 0, mpicom)
   call mpibcast(Nudge_Vwin_hi      , 1, mpir8 , 0, mpicom)
   call mpibcast(Nudge_Vwin_Hindex  , 1, mpir8 , 0, mpicom)
   call mpibcast(Nudge_Vwin_Hdelta  , 1, mpir8 , 0, mpicom)
   call mpibcast(Nudge_Vwin_Lindex  , 1, mpir8 , 0, mpicom)
   call mpibcast(Nudge_Vwin_Ldelta  , 1, mpir8 , 0, mpicom)
   call mpibcast(Nudge_Vwin_Invert,   1, mpilog, 0, mpicom)
   call mpibcast(Nudge_Method,len(Nudge_Method),mpichar,0,mpicom)
   call mpibcast(Nudge_Loc_PhysOut,1,mpilog,0,mpicom)
   call mpibcast(Nudge_Tau,1,mpir8,0,mpicom)
   call mpibcast(Nudge_CurrentStep,1,mpilog,0,mpicom)
   call mpibcast(Nudge_File_Ntime,1,mpiint,0,mpicom)

   call mpibcast(Nudge_SRF_File_Template ,len(Nudge_SRF_File_Template),mpichar,0,mpicom)
   call mpibcast(Nudge_SRF_File_Ntime    , 1, mpiint, 0, mpicom)
   call mpibcast(Nudge_SRFcoef           , 1, mpir8 , 0, mpicom)
   call mpibcast(Nudge_SRFprof           , 1, mpiint, 0, mpicom)
   call mpibcast(Nudge_SRF_Tau           , 1, mpir8,  0, mpicom)
   call mpibcast(Nudge_Land              , 1, mpilog, 0, mpicom)
   call mpibcast(Nudge_SRF_On            , 1, mpilog, 0, mpicom)
   call mpibcast(Nudge_SRF_Flux_On       , 1, mpilog, 0, mpicom)
   call mpibcast(Nudge_SRF_Q_On          , 1, mpilog, 0, mpicom)
   call mpibcast(Nudge_SRF_PSWgt_On      , 1, mpilog, 0, mpicom)
   call mpibcast(Nudge_SRF_Prec_On       , 1, mpilog, 0, mpicom)
   call mpibcast(Nudge_SRF_RadFlux_On    , 1, mpilog, 0, mpicom)
   call mpibcast(Nudge_SRF_State_On      , 1, mpilog, 0, mpicom)
   call mpibcast(Nudge_SRF_File_Present  , 1, mpilog, 0, mpicom)
   call mpibcast(Nudge_PS_OPT            , 1, mpiint, 0, mpicom)
   call mpibcast(Nudge_UV_OPT            , 1, mpiint, 0, mpicom)
   call mpibcast(Nudge_T_OPT             , 1, mpiint, 0, mpicom)
   call mpibcast(Nudge_Q_OPT             , 1, mpiint, 0, mpicom)
   call mpibcast(Nudge_NO_PBL_UV         , 1, mpiint, 0, mpicom)
   call mpibcast(Nudge_NO_PBL_T          , 1, mpiint, 0, mpicom)
   call mpibcast(Nudge_NO_PBL_Q          , 1, mpiint, 0, mpicom)
   call mpibcast(mltbc_model_path,len(mltbc_model_path)      ,mpichar,0,mpicom)
   call mpibcast(mltbc_nudge             , 1, mpilog, 0, mpicom)
   call mpibcast(mltbc_regional_on       , 1, mpilog, 0, mpicom)
   call mpibcast(mltbc_predict_option    , 1, mpiint, 0, mpicom)
   call mpibcast(mltbc_patch_biln        , 1, mpilog, 0, mpicom)
   call mpibcast(mltbc_interp_test       , 1, mpilog, 0, mpicom)
#endif

   if ( Nudge_ON .and. (Nudge_File_Ntime .ne. Nudge_Times_Per_Day) .and. (Nudge_File_Ntime .ne. 1) ) then
     write(iulog,*) 'NUDGING: Nudge_File_Ntime=',Nudge_File_Ntime
     write(iulog,*) 'NUDGING: Nudge_File_Ntime must equal to Nudge_Times_Per_Day or 1'
     call endrun('nudging_readnl:: ERROR in namelist')
   end if

   if ( Nudge_SRF_On .and. (Nudge_SRF_File_Ntime .ne. Nudge_Times_Per_Day) .and. (Nudge_SRF_File_Ntime .ne. 1) ) then
     write(iulog,*) 'NUDGING: Nudge_SRF_File_Ntime=',Nudge_SRF_File_Ntime
     write(iulog,*) 'NUDGING: Nudge_SRF_File_Ntime must equal to Nudge_Times_Per_Day or 1'
     call endrun('nudging_readnl:: ERROR in namelist')
   end if


   ! End Routine
   !------------
   return
  end subroutine ! nudging_readnl
  !================================================================


  !================================================================
  subroutine nudging_init()
   !
   ! NUDGING_INIT: Allocate space and initialize Nudging values
   !===============================================================
   use ppgrid        ,only: pver,pcols,begchunk,endchunk
   use error_messages,only: alloc_err
   use dycore        ,only: dycore_is
   use dyn_grid      ,only: get_horiz_grid_dim_d, get_dyn_grid_parm, get_horiz_grid_d
   use phys_grid     ,only: get_rlat_p,get_rlon_p,get_ncols_p,get_lat_p,get_lon_p, & 
                            get_rlat_all_p,get_rlon_all_p,get_lon_all_p,get_lat_all_p           
   use cam_history   ,only: addfld, horiz_only
   use shr_const_mod ,only: SHR_CONST_PI
   use constituents  ,only: pcnst

   ! Local values
   !----------------
   integer  Year,Month,Day,Sec
   integer  YMD1,YMD
   logical  After_Beg,Before_End
   integer  istat,lchnk,ncol,icol,ilev
   integer  hdim1_d,hdim2_d
   real(r8) :: lons(pcols), lats(pcols)
   real(r8) :: rlat,rlon
   real(r8) :: Wprof(pver)
   real(r8) :: lonp,lon0,lonn,latp,lat0,latn
   real(r8) :: Val1_p,Val2_p,Val3_p,Val4_p
   real(r8) :: Val1_0,Val2_0,Val3_0,Val4_0
   real(r8) :: Val1_n,Val2_n,Val3_n,Val4_n
   integer  :: i,j,m,n

   ! Allocate Space for Nudging data arrays
   !-----------------------------------------
   allocate(Target_U(pcols,pver,begchunk:endchunk),stat=istat)
   call alloc_err(istat,'nudging_init','Target_U',pcols*pver*((endchunk-begchunk)+1))
   allocate(Model_U(pcols,pver,begchunk:endchunk),stat=istat)
   call alloc_err(istat,'nudging_init','Model_U',pcols*pver*((endchunk-begchunk)+1))

   allocate(Target_V(pcols,pver,begchunk:endchunk),stat=istat)
   call alloc_err(istat,'nudging_init','Target_V',pcols*pver*((endchunk-begchunk)+1))
   allocate(Model_V(pcols,pver,begchunk:endchunk),stat=istat)
   call alloc_err(istat,'nudging_init','Model_V',pcols*pver*((endchunk-begchunk)+1))

   allocate(Target_T(pcols,pver,begchunk:endchunk),stat=istat)
   call alloc_err(istat,'nudging_init','Target_T',pcols*pver*((endchunk-begchunk)+1))
   allocate(Model_T(pcols,pver,begchunk:endchunk),stat=istat)
   call alloc_err(istat,'nudging_init','Model_T',pcols*pver*((endchunk-begchunk)+1))

   allocate(Target_Q(pcols,pver,begchunk:endchunk),stat=istat)
   call alloc_err(istat,'nudging_init','Target_Q',pcols*pver*((endchunk-begchunk)+1))
   allocate(Model_Q(pcols,pver,begchunk:endchunk),stat=istat)
   call alloc_err(istat,'nudging_init','Model_Q',pcols*pver*((endchunk-begchunk)+1))

   allocate(Target_PS(pcols,begchunk:endchunk),stat=istat)
   call alloc_err(istat,'nudging_init','Target_PS',pcols*((endchunk-begchunk)+1))
   allocate(Model_PS(pcols,begchunk:endchunk),stat=istat)
   call alloc_err(istat,'nudging_init','Model_PS',pcols*((endchunk-begchunk)+1))

   allocate(Target_PHIS(pcols,begchunk:endchunk),stat=istat)
   call alloc_err(istat,'nudging_init','Target_PHIS',pcols*((endchunk-begchunk)+1))
   allocate(Model_PHIS(pcols,begchunk:endchunk),stat=istat)
   call alloc_err(istat,'nudging_init','Model_PHIS',pcols*((endchunk-begchunk)+1))

   !-------------------------------------------
   ! Allocate Space for spatial dependence of
   ! Nudging Coefs and Nudging Forcing.
   !-------------------------------------------
   allocate(Nudge_Utau(pcols,pver,begchunk:endchunk),stat=istat)
   call alloc_err(istat,'nudging_init','Nudge_Utau',pcols*pver*((endchunk-begchunk)+1))
   allocate(Nudge_Ustep(pcols,pver,begchunk:endchunk),stat=istat)
   call alloc_err(istat,'nudging_init','Nudge_Ustep',pcols*pver*((endchunk-begchunk)+1))

   allocate(Nudge_Vtau(pcols,pver,begchunk:endchunk),stat=istat)
   call alloc_err(istat,'nudging_init','Nudge_Vtau',pcols*pver*((endchunk-begchunk)+1))
   allocate(Nudge_Vstep(pcols,pver,begchunk:endchunk),stat=istat)
   call alloc_err(istat,'nudging_init','Nudge_Vstep',pcols*pver*((endchunk-begchunk)+1))

   allocate(Nudge_Ttau(pcols,pver,begchunk:endchunk),stat=istat)
   call alloc_err(istat,'nudging_init','Nudge_Ttau',pcols*pver*((endchunk-begchunk)+1))
   allocate(Nudge_Tstep(pcols,pver,begchunk:endchunk),stat=istat)
   call alloc_err(istat,'nudging_init','Nudge_Tstep',pcols*pver*((endchunk-begchunk)+1))

   allocate(Nudge_Qtau(pcols,pver,begchunk:endchunk),stat=istat)
   call alloc_err(istat,'nudging_init','Nudge_Qtau',pcols*pver*((endchunk-begchunk)+1))
   allocate(Nudge_Qstep(pcols,pver,begchunk:endchunk),stat=istat)
   call alloc_err(istat,'nudging_init','Nudge_Qstep',pcols*pver*((endchunk-begchunk)+1))

   allocate(Nudge_PStau(pcols,begchunk:endchunk),stat=istat)
   call alloc_err(istat,'nudging_init','Nudge_PStau',pcols*((endchunk-begchunk)+1))
   allocate(Nudge_PSstep(pcols,begchunk:endchunk),stat=istat)
   call alloc_err(istat,'nudging_init','Nudge_PSstep',pcols*((endchunk-begchunk)+1))

   !-------------------------------------------
   ! Allocate Space for surface nudging 
   !-------------------------------------------
   if (Nudge_Land) then 
     allocate(Target_U10(pcols,begchunk:endchunk),stat=istat)
     call alloc_err(istat,'nudging_init','Target_U10',pcols*((endchunk-begchunk)+1))
     allocate(Target_V10(pcols,begchunk:endchunk),stat=istat)
     call alloc_err(istat,'nudging_init','Target_V10',pcols*((endchunk-begchunk)+1))
     allocate(Target_T2(pcols,begchunk:endchunk),stat=istat)
     call alloc_err(istat,'nudging_init','Target_T2',pcols*((endchunk-begchunk)+1))
     allocate(Target_TD2(pcols,begchunk:endchunk),stat=istat)
     call alloc_err(istat,'nudging_init','Target_TD2',pcols*((endchunk-begchunk)+1))
     allocate(Target_TS(pcols,begchunk:endchunk),stat=istat)
     call alloc_err(istat,'nudging_init','Target_TS',pcols*((endchunk-begchunk)+1))
     allocate(Target_Q2(pcols,begchunk:endchunk),stat=istat)
     call alloc_err(istat,'nudging_init','Target_Q2',pcols*((endchunk-begchunk)+1))
     allocate(Target_PRECC(pcols,begchunk:endchunk),stat=istat)
     call alloc_err(istat,'nudging_init','Target_PRECC',pcols*((endchunk-begchunk)+1))
     allocate(Target_PRECL(pcols,begchunk:endchunk),stat=istat)
     call alloc_err(istat,'nudging_init','Target_PRECL',pcols*((endchunk-begchunk)+1))
     allocate(Target_PRECSC(pcols,begchunk:endchunk),stat=istat)
     call alloc_err(istat,'nudging_init','Target_PRECSC',pcols*((endchunk-begchunk)+1))
     allocate(Target_PRECSL(pcols,begchunk:endchunk),stat=istat)
     call alloc_err(istat,'nudging_init','Target_PRECSL',pcols*((endchunk-begchunk)+1))
     allocate(Target_EVAP(pcols,begchunk:endchunk),stat=istat)
     call alloc_err(istat,'nudging_init','Target_EVAP',pcols*((endchunk-begchunk)+1))
     allocate(Target_SHFLX(pcols,begchunk:endchunk),stat=istat)
     call alloc_err(istat,'nudging_init','Target_SHFLX',pcols*((endchunk-begchunk)+1))
     allocate(Target_LHFLX(pcols,begchunk:endchunk),stat=istat)
     call alloc_err(istat,'nudging_init','Target_LHFLX',pcols*((endchunk-begchunk)+1))
     allocate(Target_NETSW(pcols,begchunk:endchunk),stat=istat)
     call alloc_err(istat,'nudging_init','Target_NETSW',pcols*((endchunk-begchunk)+1))
     allocate(Target_FLWDS(pcols,begchunk:endchunk),stat=istat)
     call alloc_err(istat,'nudging_init','Target_FLWDS',pcols*((endchunk-begchunk)+1))
     allocate(Target_FSDS(pcols,begchunk:endchunk),stat=istat)
     call alloc_err(istat,'nudging_init','Target_FSDS',pcols*((endchunk-begchunk)+1))
     allocate(Target_FSDSD(pcols,begchunk:endchunk),stat=istat)
     call alloc_err(istat,'nudging_init','Target_FSDSD',pcols*((endchunk-begchunk)+1))
     allocate(Target_FSDSUV(pcols,begchunk:endchunk),stat=istat)
     call alloc_err(istat,'nudging_init','Target_FSDSUV',pcols*((endchunk-begchunk)+1))
     allocate(Nudge_SRFtau(pcols,begchunk:endchunk),stat=istat)
     call alloc_err(istat,'nudging_init','Nudge_SRFtau',pcols*((endchunk-begchunk)+1))
   end if 

   !-----------------------------------------------------
   ! Register output fields with the cam history module
   !-----------------------------------------------------
   call addfld('Nudge_U',(/ 'lev' /),'A','m/s/s'  ,'U Nudging Tendency')
   call addfld('Nudge_V',(/ 'lev' /),'A','m/s/s'  ,'V Nudging Tendency')
   call addfld('Nudge_T',(/ 'lev' /),'A','K/s'    ,'T Nudging Tendency')
   call addfld('Nudge_Q',(/ 'lev' /),'A','kg/kg/s','Q Nudging Tendency')
   call addfld('Nudge_PS', horiz_only, 'A','Pa/s' ,'PS Nudging Tendency')

   !-----------------------------------------------------
   call addfld('Z3_bf_ndg',(/ 'lev' /), 'A','m'      ,'Geopotential Height Before Nudging')
   call addfld('U_bf_ndg',(/ 'lev' /),  'A','m/s'    ,'Zonal Wind Before Nudging')
   call addfld('V_bf_ndg',(/ 'lev' /),  'A','m/s'    ,'Meridional Wind Before Nudging')
   call addfld('T_bf_ndg',(/ 'lev' /),  'A','K'      ,'Temperature Before Nudging')
   call addfld('Q_bf_ndg',(/ 'lev' /),  'A','kg/kg'  ,'Specific Humidity Before Nudging')
   call addfld('PS_bf_ndg', horiz_only, 'A','Pa'     ,'Surface Pressure Before Nudging')

   call addfld('Z3_af_ndg',(/ 'lev' /), 'A','m'      ,'Geopotential Height After Nudging')
   call addfld('U_af_ndg',(/ 'lev' /),  'A','m/s'    ,'Zonal Wind After Nudging')
   call addfld('V_af_ndg',(/ 'lev' /),  'A','m/s'    ,'Meridional Wind After Nudging')
   call addfld('T_af_ndg',(/ 'lev' /),  'A','K'      ,'Temperature After Nudging')
   call addfld('Q_af_ndg',(/ 'lev' /),  'A','kg/kg'  ,'Specific Humidity After Nudging')
   call addfld('PS_af_ndg', horiz_only, 'A','Pa'     ,'Surface Pressure After Nudging')

   call addfld('Z3_ref',(/ 'lev' /), 'A','m'      ,'Reference for Geopotential Height')
   call addfld('U_ref',(/ 'lev' /),  'A','m/s'    ,'Reference for Zonal Wind')
   call addfld('V_ref',(/ 'lev' /),  'A','m/s'    ,'Reference for Meridional Wind')
   call addfld('T_ref',(/ 'lev' /),  'A','K/s'    ,'Reference for Temperature')
   call addfld('Q_ref',(/ 'lev' /),  'A','kg/kg'  ,'Reference for Specific Humidity')
   call addfld('PS_ref', horiz_only, 'A','Pa'     ,'Reference for Surface pressure')
   call addfld('PHIS_ref', horiz_only, 'A','m2/s2','Reference for Surface geopotential')

   call addfld('Nudge_U_vint',horiz_only,'A','kg/m/s2','Vertical integral of U Nudging Tendency')
   call addfld('Nudge_V_vint',horiz_only,'A','kg/m/s2','Vertical integral of V Nudging Tendency')
   call addfld('Nudge_T_vint',horiz_only,'A','W/m2'   ,'Vertical integral of T Nudging Tendency')
   call addfld('Nudge_Q_vint',horiz_only,'A','kg/m2/s','Vertical integral of Q Nudging Tendency')

   call addfld('Nudge_PRECC', horiz_only, 'A','m/s/s'    ,'PRECC Nudging Tendency')
   call addfld('Nudge_PRECL', horiz_only, 'A','m/s/s'    ,'PRECL Nudging Tendency')
   call addfld('Nudge_PRECSC',horiz_only, 'A','m/s/s'    ,'PRECSC Nudging Tendency')
   call addfld('Nudge_PRECSL',horiz_only, 'A','m/s/s'    ,'PRECSL Nudging Tendency')
   call addfld('Nudge_LHFLX', horiz_only, 'A','W/m2/s'   ,'LHFLX Nudging Tendency')
   call addfld('Nudge_SHFLX', horiz_only, 'A','W/m2/s'   ,'SHFLX Nudging Tendency')
   call addfld('Nudge_QFLX',  horiz_only, 'A','kg/m2/s2' ,'QFLX Nudging Tendency')
   call addfld('Nudge_SOLL',  horiz_only, 'A','W/m2/s'   ,'SOLL Nudging Tendency')
   call addfld('Nudge_SOLS',  horiz_only, 'A','W/m2/s'   ,'SOLS Nudging Tendency')
   call addfld('Nudge_SOLLD', horiz_only, 'A','W/m2/s'   ,'SOLLD Nudging Tendency')
   call addfld('Nudge_SOLSD', horiz_only, 'A','W/m2/s'   ,'SOLSD Nudging Tendency')
   call addfld('Nudge_NETSW', horiz_only, 'A','W/m2/s'   ,'NETSW Nudging Tendency')
   call addfld('Nudge_FLWDS', horiz_only, 'A','W/m2/s'   ,'FLWDS Nudging Tendency')

   call addfld('PRECC_bf_ndg',  horiz_only, 'A','m/s'     ,'Convective rain rate  Before Nudging')
   call addfld('PRECL_bf_ndg',  horiz_only, 'A','m/s'     ,'Large-scale rain rate Before Nudging')
   call addfld('PRECSC_bf_ndg', horiz_only, 'A','m/s'     ,'Convective snow rate  Before Nudging')
   call addfld('PRECSL_bf_ndg', horiz_only, 'A','m/s'     ,'Large-scale snow rate Before Nudging')
   call addfld('LHFLX_bf_ndg',  horiz_only, 'A','W/m2'    ,'Surface latent heat flux Before Nudging')
   call addfld('SHFLX_bf_ndg',  horiz_only, 'A','W/m2'    ,'Surface sensible heat flux Before Nudging')
   call addfld('QFLX_bf_ndg',   horiz_only, 'A','kg/m2/s' ,'Surface water flux Before Nudging')
   call addfld('SOLL_bf_ndg',   horiz_only, 'A','W/m2'    ,'Near IR direct solar radiation Before Nudging')
   call addfld('SOLS_bf_ndg',   horiz_only, 'A','W/m2'    ,'Visible direct solar radiation Before Nudging')
   call addfld('SOLLD_bf_ndg',  horiz_only, 'A','W/m2'    ,'Near IR diffuse solar radiation Before Nudging')
   call addfld('SOLSD_bf_ndg',  horiz_only, 'A','W/m2'    ,'Visible diffuse solar radiation Before Nudging')
   call addfld('NETSW_bf_ndg',  horiz_only, 'A','W/m2'    ,'Net solar radiation Before Nudging')
   call addfld('FLWDS_bf_ndg',  horiz_only, 'A','W/m2'    ,'Downward longwave radiation Before Nudging')

   call addfld('PRECC_af_ndg',  horiz_only, 'A','m/s'     ,'Convective rain rate After Nudging')
   call addfld('PRECL_af_ndg',  horiz_only, 'A','m/s'     ,'Large-scale rain rate After Nudging')
   call addfld('PRECSC_af_ndg', horiz_only, 'A','m/s'     ,'Convective snow rate After Nudging')
   call addfld('PRECSL_af_ndg', horiz_only, 'A','m/s'     ,'Large-scale snow rate After Nudging')
   call addfld('LHFLX_af_ndg',  horiz_only, 'A','W/m2'    ,'Surface latent heat flux After Nudging')
   call addfld('SHFLX_af_ndg',  horiz_only, 'A','W/m2'    ,'Surface sensible heat flux After Nudging')
   call addfld('QFLX_af_ndg',   horiz_only, 'A','kg/m2/s' ,'Surface water flux After Nudging')
   call addfld('SOLL_af_ndg',   horiz_only, 'A','W/m2'    ,'Near IR direct solar radiation After Nudging')
   call addfld('SOLS_af_ndg',   horiz_only, 'A','W/m2'    ,'Visible direct solar radiation After Nudging')
   call addfld('SOLLD_af_ndg',  horiz_only, 'A','W/m2'    ,'Near IR diffuse solar radiation After Nudging')
   call addfld('SOLSD_af_ndg',  horiz_only, 'A','W/m2'    ,'Visible diffuse solar radiation After Nudging')
   call addfld('NETSW_af_ndg',  horiz_only, 'A','W/m2'    ,'Net solar radiation After Nudging')
   call addfld('FLWDS_af_ndg',  horiz_only, 'A','W/m2'    ,'Downward longwave radiation After Nudging')

   call addfld('PRECC_ref',  horiz_only, 'A','m/s'     ,'Reference convective rain rate')
   call addfld('PRECL_ref',  horiz_only, 'A','m/s'     ,'Reference large-scale rain rate')
   call addfld('PRECSC_ref', horiz_only, 'A','m/s'     ,'Reference convective snow rate')
   call addfld('PRECSL_ref', horiz_only, 'A','m/s'     ,'Reference large-scale snow rate')
   call addfld('LHFLX_ref',  horiz_only, 'A','W/m2'    ,'Reference surface latent heat flux')
   call addfld('SHFLX_ref',  horiz_only, 'A','W/m2'    ,'Reference surface sensible heat flux')
   call addfld('QFLX_ref',   horiz_only, 'A','kg/m2/s' ,'Reference surface water flux')
   call addfld('SOLL_ref',   horiz_only, 'A','W/m2'    ,'Reference near IR direct solar radiation')
   call addfld('SOLS_ref',   horiz_only, 'A','W/m2'    ,'Reference visible direct solar radiation')
   call addfld('SOLLD_ref',  horiz_only, 'A','W/m2'    ,'Reference near IR diffuse solar radiation')
   call addfld('SOLSD_ref',  horiz_only, 'A','W/m2'    ,'Reference visible diffuse solar radiation')
   call addfld('NETSW_ref',  horiz_only, 'A','W/m2'    ,'Reference net solar radiation')
   call addfld('FLWDS_ref',  horiz_only, 'A','W/m2'    ,'Reference downward longwave radiation')

   !-----------------------------------------------------
   ! Register output fields for Machine Learning Nudging 
   !-----------------------------------------------------
   call addfld('DON_RGD_UERR',  (/ 'lev' /), 'A', 'm/s'   , 'DeeONet Interpolation Error for U' )
   call addfld('DON_RGD_VERR',  (/ 'lev' /), 'A', 'm/s'   , 'DeeONet Interpolation Error for V' )
   call addfld('DON_RGD_TERR',  (/ 'lev' /), 'A', 'K'     , 'DeeONet Interpolation Error for T' )
   call addfld('DON_RGD_QERR',  (/ 'lev' /), 'A', 'kg/kg' , 'DeeONet Interpolation Error for Q' )
   call addfld('DON_RGD_PSERR',  horiz_only, 'A', 'Pa'    , 'DeeONet Interpolation Error for PS' )

   !-----------------------------------------
   ! Values initialized only by masterproc
   !-----------------------------------------
   if(masterproc) then

     ! Set the Stepping intervals for Model and Nudging values
     ! Ensure that the Model_Step is not smaller than one timestep
     !  and not larger then the Nudge_Step.
     !--------------------------------------------------------
     Model_Step=86400/Model_Times_Per_Day
     Nudge_Step=86400/Nudge_Times_Per_Day

     if(Model_Step.lt.dtime) then
       write(iulog,*) ' '
       write(iulog,*) 'NUDGING: Model_Step cannot be less than a model timestep'
       write(iulog,*) 'NUDGING:  Setting Model_Step=dtime , dtime=',dtime
       write(iulog,*) ' '
       Model_Step=dtime
     endif
     if(Model_Step.gt.Nudge_Step) then
       write(iulog,*) ' '
       write(iulog,*) 'NUDGING: Model_Step cannot be more than Nudge_Step'
       write(iulog,*) 'NUDGING:  Setting Model_Step=Nudge_Step, Nudge_Step=',Nudge_Step
       write(iulog,*) ' '
       Model_Step=Nudge_Step
     endif

     ! Initialize column and level dimensions
     !--------------------------------------------------------
     call get_horiz_grid_dim_d(hdim1_d,hdim2_d)
     Nudge_nlon=hdim1_d
     Nudge_nlat=hdim2_d
     Nudge_ncol=hdim1_d*hdim2_d
     Nudge_nlev=pver
!DIAG
     Nudge_slat=Nudge_nlat-1
!DIAG

     ! Check the time relative to the nudging window
     !------------------------------------------------
     call get_curr_date(Year,Month,Day,Sec)

     YMD=(Year*10000) + (Month*100) + Day
     YMD1=(Nudge_Beg_Year*10000) + (Nudge_Beg_Month*100) + Nudge_Beg_Day

     call timemgr_time_ge(YMD1,Nudge_Beg_Sec,         &
                          YMD ,Sec          ,After_Beg)
     YMD1=(Nudge_End_Year*10000) + (Nudge_End_Month*100) + Nudge_End_Day

     call timemgr_time_ge(YMD ,Sec          ,          &
                          YMD1,Nudge_End_Sec,Before_End)

     if((After_Beg).and.(Before_End)) then
       ! Set Time indicies so that the next call to
       ! timestep_init will initialize the data arrays.
       !--------------------------------------------
       Model_Next_Year =Year
       Model_Next_Month=Month
       Model_Next_Day  =Day
       Model_Next_Sec  =(Sec/Model_Step)*Model_Step
       Nudge_Next_Year =Year
       Nudge_Next_Month=Month
       Nudge_Next_Day  =Day
       Nudge_Next_Sec  =(Sec/Nudge_Step)*Nudge_Step

     elseif(.not.After_Beg) then
       ! Set Time indicies to Nudging start,
       ! timestep_init will initialize the data arrays.
       !--------------------------------------------
       Model_Next_Year =Nudge_Beg_Year
       Model_Next_Month=Nudge_Beg_Month
       Model_Next_Day  =Nudge_Beg_Day
       Model_Next_Sec  =Nudge_Beg_Sec
       Nudge_Next_Year =Nudge_Beg_Year
       Nudge_Next_Month=Nudge_Beg_Month
       Nudge_Next_Day  =Nudge_Beg_Day
       Nudge_Next_Sec  =Nudge_Beg_Sec

     elseif(.not.Before_End) then
       ! Nudging will never occur, so switch it off
       !--------------------------------------------
       mltbc_nudge = .false.
       Nudge_Model    = .false.
       Nudge_ON       = .false.
       Nudge_Land     = .false.
       Nudge_SRF_On   = .false. 
       write(iulog,*) ' '
       write(iulog,*) 'NUDGING: WARNING - Nudging has been requested by it will'
       write(iulog,*) 'NUDGING:           never occur for the given time values'
       write(iulog,*) ' '
     endif

     ! Initialize values for window function
     !----------------------------------------
     lonp= 180.
     lon0=   0.
     lonn=-180.
     latp=  90.-Nudge_Hwin_lat0
     lat0=   0.
     latn= -90.-Nudge_Hwin_lat0

     Nudge_Hwin_lonWidthH=Nudge_Hwin_lonWidth/2.
     Nudge_Hwin_latWidthH=Nudge_Hwin_latWidth/2.

     Val1_p=(1.+tanh((Nudge_Hwin_lonWidthH+lonp)/Nudge_Hwin_lonDelta))/2.
     Val2_p=(1.+tanh((Nudge_Hwin_lonWidthH-lonp)/Nudge_Hwin_lonDelta))/2.
     Val3_p=(1.+tanh((Nudge_Hwin_latWidthH+latp)/Nudge_Hwin_latDelta))/2.
     Val4_p=(1.+tanh((Nudge_Hwin_latWidthH-latp)/Nudge_Hwin_latDelta))/2.

     Val1_0=(1.+tanh((Nudge_Hwin_lonWidthH+lon0)/Nudge_Hwin_lonDelta))/2.
     Val2_0=(1.+tanh((Nudge_Hwin_lonWidthH-lon0)/Nudge_Hwin_lonDelta))/2.
     Val3_0=(1.+tanh((Nudge_Hwin_latWidthH+lat0)/Nudge_Hwin_latDelta))/2.
     Val4_0=(1.+tanh((Nudge_Hwin_latWidthH-lat0)/Nudge_Hwin_latDelta))/2.

     Val1_n=(1.+tanh((Nudge_Hwin_lonWidthH+lonn)/Nudge_Hwin_lonDelta))/2.
     Val2_n=(1.+tanh((Nudge_Hwin_lonWidthH-lonn)/Nudge_Hwin_lonDelta))/2.
     Val3_n=(1.+tanh((Nudge_Hwin_latWidthH+latn)/Nudge_Hwin_latDelta))/2.
     Val4_n=(1.+tanh((Nudge_Hwin_latWidthH-latn)/Nudge_Hwin_latDelta))/2.

     Nudge_Hwin_max=     Val1_0*Val2_0*Val3_0*Val4_0
     Nudge_Hwin_min=min((Val1_p*Val2_p*Val3_n*Val4_n), &
                        (Val1_p*Val2_p*Val3_p*Val4_p), &
                        (Val1_n*Val2_n*Val3_n*Val4_n), &
                        (Val1_n*Val2_n*Val3_p*Val4_p))

     ! Initialization is done,
     !--------------------------
     Nudge_Initialized=.true.

     ! Check that this is a valid DYCORE model
     !------------------------------------------
     if((.not.dycore_is('UNSTRUCTURED')).and. &
        (.not.dycore_is('EUL')         ).and. &
        (.not.dycore_is('LR')          )      ) then
       call endrun('NUDGING IS CURRENTLY ONLY CONFIGURED FOR CAM-SE, FV, or EUL')
     endif

     ! Informational Output
     !---------------------------
     write(iulog,*) ' '
     write(iulog,*) '---------------------------------------------------------'
     write(iulog,*) '  MODEL NUDGING INITIALIZED WITH THE FOLLOWING SETTINGS: '
     write(iulog,*) '---------------------------------------------------------'
     write(iulog,*) 'NUDGING: Nudge_Model=',Nudge_Model
     write(iulog,*) 'NUDGING: Nudge_Allow_Missing_File=',Nudge_Allow_Missing_File
     write(iulog,*) 'NUDGING: Nudge_Pdep_Weight_On=',Nudge_Pdep_Weight_On
     write(iulog,*) 'NUDGING: Nudge_Lin_Relax_On=', Nudge_Lin_Relax_On
     write(iulog,*) 'NUDGING: Nudge_PS_On=',Nudge_PS_On
     write(iulog,*) 'NUDGING: Nudge_Q_Adjust_On=',Nudge_Q_Adjust_On
     write(iulog,*) 'NUDGING: Nudge_PS_Adjust_On=',Nudge_PS_Adjust_On
     write(iulog,*) 'NUDGING: Nudge_Land=',Nudge_Land
     write(iulog,*) 'NUDGING: Nudge_SRF_Flux_On=',Nudge_SRF_Flux_On
     write(iulog,*) 'NUDGING: Nudge_SRF_Q_On=',Nudge_SRF_Q_On
     write(iulog,*) 'NUDGING: Nudge_SRF_PSWgt_On=',Nudge_SRF_PSWgt_On
     write(iulog,*) 'NUDGING: Nudge_SRF_Prec_On=',Nudge_SRF_Prec_On
     write(iulog,*) 'NUDGING: Nudge_SRF_RadFlux_On=',Nudge_SRF_RadFlux_On
     write(iulog,*) 'NUDGING: Nudge_SRF_State_On=',Nudge_SRF_State_On
     write(iulog,*) 'NUDGING: Nudge_Path=',Nudge_Path
     write(iulog,*) 'NUDGING: Nudge_File_Template=',Nudge_File_Template
     write(iulog,*) 'NUDGING: Nudge_SRF_File_Template=',Nudge_SRF_File_Template
     write(iulog,*) 'NUDGING: Nudge_Times_Per_Day=',Nudge_Times_Per_Day
     write(iulog,*) 'NUDGING: Model_Times_Per_Day=',Model_Times_Per_Day
     write(iulog,*) 'NUDGING: Nudge_Step=',Nudge_Step
     write(iulog,*) 'NUDGING: Model_Step=',Model_Step
     write(iulog,*) 'NUDGING: Nudge_Ucoef  =',Nudge_Ucoef
     write(iulog,*) 'NUDGING: Nudge_Vcoef  =',Nudge_Vcoef
     write(iulog,*) 'NUDGING: Nudge_Qcoef  =',Nudge_Qcoef
     write(iulog,*) 'NUDGING: Nudge_Tcoef  =',Nudge_Tcoef
     write(iulog,*) 'NUDGING: Nudge_PScoef =',Nudge_PScoef
     write(iulog,*) 'NUDGING: Nudge_SRFcoef=',Nudge_SRFcoef
     write(iulog,*) 'NUDGING: Nudge_Uprof  =',Nudge_Uprof
     write(iulog,*) 'NUDGING: Nudge_Vprof  =',Nudge_Vprof
     write(iulog,*) 'NUDGING: Nudge_Qprof  =',Nudge_Qprof
     write(iulog,*) 'NUDGING: Nudge_Tprof  =',Nudge_Tprof
     write(iulog,*) 'NUDGING: Nudge_PSprof =',Nudge_PSprof
     write(iulog,*) 'NUDGING: Nudge_SRFprof=',Nudge_SRFprof
     write(iulog,*) 'NUDGING: Nudge_Beg_Year =',Nudge_Beg_Year
     write(iulog,*) 'NUDGING: Nudge_Beg_Month=',Nudge_Beg_Month
     write(iulog,*) 'NUDGING: Nudge_Beg_Day  =',Nudge_Beg_Day
     write(iulog,*) 'NUDGING: Nudge_End_Year =',Nudge_End_Year
     write(iulog,*) 'NUDGING: Nudge_End_Month=',Nudge_End_Month
     write(iulog,*) 'NUDGING: Nudge_End_Day  =',Nudge_End_Day
     write(iulog,*) 'NUDGING: Nudge_Hwin_lo       =',Nudge_Hwin_lo
     write(iulog,*) 'NUDGING: Nudge_Hwin_hi       =',Nudge_Hwin_hi
     write(iulog,*) 'NUDGING: Nudge_Hwin_lat0     =',Nudge_Hwin_lat0
     write(iulog,*) 'NUDGING: Nudge_Hwin_latWidth =',Nudge_Hwin_latWidth
     write(iulog,*) 'NUDGING: Nudge_Hwin_latDelta =',Nudge_Hwin_latDelta
     write(iulog,*) 'NUDGING: Nudge_Hwin_lon0     =',Nudge_Hwin_lon0
     write(iulog,*) 'NUDGING: Nudge_Hwin_lonWidth =',Nudge_Hwin_lonWidth
     write(iulog,*) 'NUDGING: Nudge_Hwin_lonDelta =',Nudge_Hwin_lonDelta
     write(iulog,*) 'NUDGING: Nudge_Hwin_Invert   =',Nudge_Hwin_Invert
     write(iulog,*) 'NUDGING: Nudge_Vwin_lo       =',Nudge_Vwin_lo
     write(iulog,*) 'NUDGING: Nudge_Vwin_hi       =',Nudge_Vwin_hi
     write(iulog,*) 'NUDGING: Nudge_Vwin_Hindex   =',Nudge_Vwin_Hindex
     write(iulog,*) 'NUDGING: Nudge_Vwin_Hdelta   =',Nudge_Vwin_Hdelta
     write(iulog,*) 'NUDGING: Nudge_Vwin_Lindex   =',Nudge_Vwin_Lindex
     write(iulog,*) 'NUDGING: Nudge_Vwin_Ldelta   =',Nudge_Vwin_Ldelta
     write(iulog,*) 'NUDGING: Nudge_Vwin_Invert   =',Nudge_Vwin_Invert
     write(iulog,*) 'NUDGING: Nudge_Hwin_latWidthH=',Nudge_Hwin_latWidthH
     write(iulog,*) 'NUDGING: Nudge_Hwin_lonWidthH=',Nudge_Hwin_lonWidthH
     write(iulog,*) 'NUDGING: Nudge_Hwin_max      =',Nudge_Hwin_max
     write(iulog,*) 'NUDGING: Nudge_Hwin_min      =',Nudge_Hwin_min
     write(iulog,*) 'NUDGING: Nudge_Initialized   =',Nudge_Initialized
     write(iulog,*) 'NUDGING: Nudge_Method        =',Nudge_Method
     write(iulog,*) 'NUDGING: Nudge_Loc_PhysOut   =',Nudge_Loc_PhysOut
     write(iulog,*) 'NUDGING: Nudge_Tau           =',Nudge_Tau
     write(iulog,*) 'NUDGING: Nudge_SRF_Tau       =',Nudge_SRF_Tau
     write(iulog,*) 'NUDGING: Nudge_CurrentStep   =',Nudge_CurrentStep
     write(iulog,*) 'NUDGING: Nudge_File_Ntime    =',Nudge_File_Ntime
     write(iulog,*) 'NUDGING: Nudge_PS_OPT        =',Nudge_PS_OPT
     write(iulog,*) 'NUDGING: Nudge_UV_OPT        =',Nudge_UV_OPT
     write(iulog,*) 'NUDGING: Nudge_T_OPT         =',Nudge_T_OPT
     write(iulog,*) 'NUDGING: Nudge_Q_OPT         =',Nudge_Q_OPT
     write(iulog,*) 'NUDGING: Nudge_NO_PBL_UV     =',Nudge_NO_PBL_UV
     write(iulog,*) 'NUDGING: Nudge_NO_PBL_T      =',Nudge_NO_PBL_T
     write(iulog,*) 'NUDGING: Nudge_NO_PBL_Q      =',Nudge_NO_PBL_Q
     write(iulog,*) 'NUDGING: mltbc_nudge         =',mltbc_nudge
     write(iulog,*) 'NUDGING: mltbc_regional_on   =',mltbc_regional_on
     write(iulog,*) 'NUDGING: mltbc_predict_option=',mltbc_predict_option
     write(iulog,*) 'NUDGING: mltbc_interp_test   =',mltbc_interp_test
     write(iulog,*) 'NUDGING: mltbc_patch_biln    =',mltbc_patch_biln
     write(iulog,*) 'NUDGING: mltbc_model_path    =',mltbc_model_path
     write(iulog,*) ' '
     write(iulog,*) ' '

   endif ! (masterproc) then

   ! Broadcast other variables that have changed
   !---------------------------------------------
#ifdef SPMD
   call mpibcast(Model_Step          , 1, mpiint, 0, mpicom)
   call mpibcast(Nudge_Step          , 1, mpiint, 0, mpicom)
   call mpibcast(Model_Next_Year     , 1, mpiint, 0, mpicom)
   call mpibcast(Model_Next_Month    , 1, mpiint, 0, mpicom)
   call mpibcast(Model_Next_Day      , 1, mpiint, 0, mpicom)
   call mpibcast(Model_Next_Sec      , 1, mpiint, 0, mpicom)
   call mpibcast(Nudge_Next_Year     , 1, mpiint, 0, mpicom)
   call mpibcast(Nudge_Next_Month    , 1, mpiint, 0, mpicom)
   call mpibcast(Nudge_Next_Day      , 1, mpiint, 0, mpicom)
   call mpibcast(Nudge_Next_Sec      , 1, mpiint, 0, mpicom)
   call mpibcast(Nudge_Model         , 1, mpilog, 0, mpicom)
   call mpibcast(Nudge_Allow_Missing_File, 1, mpilog, 0, mpicom)
   call mpibcast(Nudge_Pdep_Weight_On, 1, mpilog, 0, mpicom) 
   call mpibcast(Nudge_Lin_Relax_On  , 1, mpilog, 0, mpicom)
   call mpibcast(Nudge_PS_Adjust_On  , 1, mpilog, 0, mpicom)
   call mpibcast(Nudge_PS_On         , 1, mpilog, 0, mpicom)
   call mpibcast(Nudge_Q_Adjust_On   , 1, mpilog, 0, mpicom)
   call mpibcast(Nudge_Land          , 1, mpilog, 0, mpicom)
   call mpibcast(Nudge_SRF_On        , 1, mpilog, 0, mpicom)
   call mpibcast(Nudge_SRF_Flux_On   , 1, mpilog, 0, mpicom)
   call mpibcast(Nudge_SRF_Q_On      , 1, mpilog, 0, mpicom)
   call mpibcast(Nudge_SRF_PSWgt_On  , 1, mpilog, 0, mpicom)
   call mpibcast(Nudge_SRF_Prec_On   , 1, mpilog, 0, mpicom)
   call mpibcast(Nudge_SRF_RadFlux_On, 1, mpilog, 0, mpicom)
   call mpibcast(Nudge_SRF_State_On  , 1, mpilog, 0, mpicom)
   call mpibcast(Nudge_ON            , 1, mpilog, 0, mpicom)
   call mpibcast(Nudge_Initialized   , 1, mpilog, 0, mpicom)
   call mpibcast(Nudge_ncol          , 1, mpiint, 0, mpicom)
   call mpibcast(Nudge_nlev          , 1, mpiint, 0, mpicom)
   call mpibcast(Nudge_nlon          , 1, mpiint, 0, mpicom)
   call mpibcast(Nudge_nlat          , 1, mpiint, 0, mpicom)
   call mpibcast(Nudge_Hwin_max      , 1, mpir8 , 0, mpicom)
   call mpibcast(Nudge_Hwin_min      , 1, mpir8 , 0, mpicom)
   call mpibcast(Nudge_Hwin_lonWidthH, 1, mpir8 , 0, mpicom)
   call mpibcast(Nudge_Hwin_latWidthH, 1, mpir8 , 0, mpicom)
   call mpibcast(Nudge_PS_OPT        , 1, mpiint, 0, mpicom)
   call mpibcast(Nudge_UV_OPT        , 1, mpiint, 0, mpicom)
   call mpibcast(Nudge_T_OPT         , 1, mpiint, 0, mpicom)
   call mpibcast(Nudge_Q_OPT         , 1, mpiint, 0, mpicom)
   call mpibcast(Nudge_NO_PBL_UV     , 1, mpiint, 0, mpicom)
   call mpibcast(Nudge_NO_PBL_T      , 1, mpiint, 0, mpicom)
   call mpibcast(Nudge_NO_PBL_Q      , 1, mpiint, 0, mpicom)
   call mpibcast(mltbc_nudge         , 1, mpilog, 0, mpicom)
   call mpibcast(mltbc_regional_on   , 1, mpilog, 0, mpicom)
   call mpibcast(mltbc_predict_option, 1, mpiint, 0, mpicom)
   call mpibcast(mltbc_patch_biln    , 1, mpilog, 0, mpicom)
   call mpibcast(mltbc_interp_test   , 1, mpilog, 0, mpicom)
!DIAG
   call mpibcast(Nudge_slat       , 1, mpiint, 0, mpicom)
!DIAG
#endif

   ! Initialize Nudging Coeffcient profiles in local arrays
   ! Load zeros into nudging arrays
   !------------------------------------------------------
   do lchnk=begchunk,endchunk
     ncol=get_ncols_p(lchnk)
     do icol=1,ncol
       rlat=get_rlat_p(lchnk,icol)*180._r8/SHR_CONST_PI
       rlon=get_rlon_p(lchnk,icol)*180._r8/SHR_CONST_PI

       call nudging_set_profile(rlat,rlon,Nudge_Uprof, Wprof,pver)
       Nudge_Utau(icol,:,lchnk)=Wprof(:)

       call nudging_set_profile(rlat,rlon,Nudge_Vprof, Wprof,pver)
       Nudge_Vtau(icol,:,lchnk)=Wprof(:)

       call nudging_set_profile(rlat,rlon,Nudge_Tprof, Wprof,pver)
       Nudge_Ttau(icol,:,lchnk)=Wprof(:)

       call nudging_set_profile(rlat,rlon,Nudge_Qprof, Wprof,pver)
       Nudge_Qtau(icol,:,lchnk)=Wprof(:)

       Nudge_PStau(icol,lchnk)=nudging_set_PSprofile(rlat,rlon,Nudge_PSprof)

     end do

     if (mltbc_nudge) then
       Nudge_Utau(:ncol,:pver,lchnk) = Nudge_Utau(:ncol,:pver,lchnk) * Nudge_Ucoef
       Nudge_Vtau(:ncol,:pver,lchnk) = Nudge_Vtau(:ncol,:pver,lchnk) * Nudge_Vcoef
       Nudge_Ttau(:ncol,:pver,lchnk) = Nudge_Ttau(:ncol,:pver,lchnk) * Nudge_Tcoef
       Nudge_Qtau(:ncol,:pver,lchnk) = Nudge_Qtau(:ncol,:pver,lchnk) * Nudge_Qcoef
       Nudge_PStau(:ncol,lchnk)      = Nudge_PStau(:ncol,lchnk)      * Nudge_PScoef
     else 

       if  (Nudge_Tau .le. 0._r8) then

         Nudge_Utau(:ncol,:pver,lchnk) =                             &
         Nudge_Utau(:ncol,:pver,lchnk) * Nudge_Ucoef/float(Nudge_Step)

         Nudge_Vtau(:ncol,:pver,lchnk) =                             &
         Nudge_Vtau(:ncol,:pver,lchnk) * Nudge_Vcoef/float(Nudge_Step)

         Nudge_Ttau(:ncol,:pver,lchnk) =                             &
         Nudge_Ttau(:ncol,:pver,lchnk) * Nudge_Tcoef/float(Nudge_Step)

         Nudge_Qtau(:ncol,:pver,lchnk) =                             &
         Nudge_Qtau(:ncol,:pver,lchnk) * Nudge_Qcoef/float(Nudge_Step)

         Nudge_PStau(:ncol,lchnk)=                              &
         Nudge_PStau(:ncol,lchnk)* Nudge_PScoef/float(Nudge_Step)

       else          ! use Nudge_Tau directy as relaxation timescale

         Nudge_Utau(:ncol,:pver,lchnk) =                        &
         Nudge_Utau(:ncol,:pver,lchnk) * Nudge_Ucoef / Nudge_Tau / sec_per_hour

         Nudge_Vtau(:ncol,:pver,lchnk) =                        &
         Nudge_Vtau(:ncol,:pver,lchnk) * Nudge_Vcoef / Nudge_Tau / sec_per_hour 

         Nudge_Ttau(:ncol,:pver,lchnk) =                        &
         Nudge_Ttau(:ncol,:pver,lchnk) * Nudge_Tcoef / Nudge_Tau / sec_per_hour 

         Nudge_Qtau(:ncol,:pver,lchnk) =                        &
         Nudge_Qtau(:ncol,:pver,lchnk) * Nudge_Qcoef / Nudge_Tau / sec_per_hour 

         Nudge_PStau(:ncol,lchnk) =                        &
         Nudge_PStau(:ncol,lchnk) * Nudge_PScoef / Nudge_Tau / sec_per_hour 

       end if
     end if 

     if (Nudge_PSprof .ne. 0) then
        Nudge_PS_On    = .true.
     else
        Nudge_PS_On    = .false.
     end if

     Nudge_Ustep(:pcols,:pver,lchnk)=0._r8
     Target_U(:pcols,:pver,lchnk)=0._r8

     Nudge_Vstep(:pcols,:pver,lchnk)=0._r8
     Target_V(:pcols,:pver,lchnk)=0._r8

     Nudge_Tstep(:pcols,:pver,lchnk)=0._r8
     Target_T(:pcols,:pver,lchnk)=0._r8

     Nudge_Qstep(:pcols,:pver,lchnk)=0._r8
     Target_Q(:pcols,:pver,lchnk)=0._r8

     Nudge_PSstep(:pcols,lchnk)=0._r8
     Target_PS(:pcols,lchnk)=0._r8
     Target_PHIS(:pcols,lchnk)=0._r8

     if (Nudge_Land) then

       do icol=1,ncol
         rlat=get_rlat_p(lchnk,icol)*180._r8/SHR_CONST_PI
         rlon=get_rlon_p(lchnk,icol)*180._r8/SHR_CONST_PI
         Nudge_SRFtau(icol,lchnk)=nudging_set_SRFprofile(rlat,rlon,Nudge_SRFprof)
       end do

       if (Nudge_SRF_Tau .le. 0._r8) then
         Nudge_SRFtau(:ncol,lchnk) =  &
         Nudge_SRFtau(:ncol,lchnk) * Nudge_SRFcoef /float(Nudge_Step)
       else 
         Nudge_SRFtau(:ncol,lchnk) =  &
         Nudge_SRFtau(:ncol,lchnk) * Nudge_SRFcoef / Nudge_SRF_Tau / sec_per_hour
       end if 

       Target_U10(:pcols,lchnk)=0._r8
       Target_V10(:pcols,lchnk)=0._r8
       Target_T2(:pcols,lchnk)=0._r8
       Target_TD2(:pcols,lchnk)=0._r8
       Target_TS(:pcols,lchnk)=0._r8
       Target_Q2(:pcols,lchnk)=0._r8
       Target_EVAP(:pcols,lchnk)=0._r8
       Target_LHFLX(:pcols,lchnk)=0._r8
       Target_SHFLX(:pcols,lchnk)=0._r8
       Target_PRECC(:pcols,lchnk)=0._r8
       Target_PRECL(:pcols,lchnk)=0._r8
       Target_PRECSC(:pcols,lchnk)=0._r8
       Target_PRECSL(:pcols,lchnk)=0._r8
       Target_NETSW(:pcols,lchnk)=0._r8
       Target_FLWDS(:pcols,lchnk)=0._r8
       Target_FSDS(:pcols,lchnk)=0._r8
       Target_FSDSD(:pcols,lchnk)=0._r8
       Target_FSDSUV(:pcols,lchnk)=0._r8

     end if 

   end do  ! lchnk loop 

   first_file     = .true.
   first_srf_file = .true.
   select case (Nudge_Method)
      case ('Step')
           ! use Xanal directly, no interpolation is needed
           !-----------------------------------------------
      case ('IMT')
           ! use Xanal directly, no interpolation is needed
           !-----------------------------------------------
      case ('Linear')
           allocate(INTP_U(pcols,pver,begchunk:endchunk,2),stat=istat)
           call alloc_err(istat,'nudging_init','INTP_U',2*pcols*pver*((endchunk-begchunk)+1))

           allocate(INTP_V(pcols,pver,begchunk:endchunk,2),stat=istat)
           call alloc_err(istat,'nudging_init','INTP_V',2*pcols*pver*((endchunk-begchunk)+1))

           allocate(INTP_T(pcols,pver,begchunk:endchunk,2),stat=istat)
           call alloc_err(istat,'nudging_init','INTP_T',2*pcols*pver*((endchunk-begchunk)+1))

           allocate(INTP_Q(pcols,pver,begchunk:endchunk,2),stat=istat)
           call alloc_err(istat,'nudging_init','INTP_Q',2*pcols*pver*((endchunk-begchunk)+1))

           allocate(INTP_PS(pcols,begchunk:endchunk,2),stat=istat)
           call alloc_err(istat,'nudging_init','INTP_PS',2*pcols*((endchunk-begchunk)+1))

           allocate(INTP_PHIS(pcols,begchunk:endchunk,2),stat=istat)
           call alloc_err(istat,'nudging_init','INTP_PHIS',2*pcols*((endchunk-begchunk)+1))

           if (Nudge_Land) then
              allocate(INTP_U10(pcols,begchunk:endchunk,2),stat=istat)
              call alloc_err(istat,'nudging_init','INTP_U10',2*pcols*((endchunk-begchunk)+1))
              allocate(INTP_V10(pcols,begchunk:endchunk,2),stat=istat)
              call alloc_err(istat,'nudging_init','INTP_V10',2*pcols*((endchunk-begchunk)+1))
              allocate(INTP_T2(pcols,begchunk:endchunk,2),stat=istat)
              call alloc_err(istat,'nudging_init','INTP_T2',2*pcols*((endchunk-begchunk)+1))
              allocate(INTP_TD2(pcols,begchunk:endchunk,2),stat=istat)
              call alloc_err(istat,'nudging_init','INTP_TD2',2*pcols*((endchunk-begchunk)+1))
              allocate(INTP_TS(pcols,begchunk:endchunk,2),stat=istat)
              call alloc_err(istat,'nudging_init','INTP_TS',2*pcols*((endchunk-begchunk)+1))
              allocate(INTP_PRECC(pcols,begchunk:endchunk,2),stat=istat)
              call alloc_err(istat,'nudging_init','INTP_PRECC',2*pcols*((endchunk-begchunk)+1))
              allocate(INTP_PRECL(pcols,begchunk:endchunk,2),stat=istat)
              call alloc_err(istat,'nudging_init','INTP_PRECL',2*pcols*((endchunk-begchunk)+1))
              allocate(INTP_PRECSC(pcols,begchunk:endchunk,2),stat=istat)
              call alloc_err(istat,'nudging_init','INTP_PRECSC',2*pcols*((endchunk-begchunk)+1))
              allocate(INTP_PRECSL(pcols,begchunk:endchunk,2),stat=istat)
              call alloc_err(istat,'nudging_init','INTP_PRECSL',2*pcols*((endchunk-begchunk)+1))
              allocate(INTP_EVAP(pcols,begchunk:endchunk,2),stat=istat)
              call alloc_err(istat,'nudging_init','INTP_EVAP',2*pcols*((endchunk-begchunk)+1))
              allocate(INTP_SHFLX(pcols,begchunk:endchunk,2),stat=istat)
              call alloc_err(istat,'nudging_init','INTP_SHFLX',2*pcols*((endchunk-begchunk)+1))
              allocate(INTP_LHFLX(pcols,begchunk:endchunk,2),stat=istat)
              call alloc_err(istat,'nudging_init','INTP_LHFLX',2*pcols*((endchunk-begchunk)+1))
              allocate(INTP_FSNS(pcols,begchunk:endchunk,2),stat=istat)
              call alloc_err(istat,'nudging_init','INTP_FSNS',2*pcols*((endchunk-begchunk)+1))
              allocate(INTP_FLDS(pcols,begchunk:endchunk,2),stat=istat)
              call alloc_err(istat,'nudging_init','INTP_FLDS',2*pcols*((endchunk-begchunk)+1))
              allocate(INTP_FSDS(pcols,begchunk:endchunk,2),stat=istat)
              call alloc_err(istat,'nudging_init','INTP_FSDS',2*pcols*((endchunk-begchunk)+1))
              allocate(INTP_FSDSD(pcols,begchunk:endchunk,2),stat=istat)
              call alloc_err(istat,'nudging_init','INTP_FSDSD',2*pcols*((endchunk-begchunk)+1))
              allocate(INTP_FSDSUV(pcols,begchunk:endchunk,2),stat=istat)
              call alloc_err(istat,'nudging_init','INTP_FSDSUV',2*pcols*((endchunk-begchunk)+1))
            end if

      case ('MLTBC')
           ! use Xanal directly, no interpolation is needed
           !-----------------------------------------------
           allocate(Model_rlat(Nudge_ncol),stat=istat)
           call alloc_err(istat,'Machine Learning NUDGING','Model_rlat',Nudge_ncol)
           allocate(Model_rlon(Nudge_ncol),stat=istat)
           call alloc_err(istat,'Machine Learning NUDGING','Model_rlon',Nudge_ncol)
           allocate(Model_Var(Nudge_ncol,Nudge_nlev),stat=istat)
           call alloc_err(istat,'nudging_init','Model_Var',Nudge_ncol*Nudge_nlev)

           !extract the lat/lon/weight info            
           call get_horiz_grid_d(Nudge_ncol,clat_d_out=Model_rlat,clon_d_out=Model_rlon)

#ifdef SPMD
           call mpibcast(Model_rlat, Nudge_ncol, mpir8,   0, mpicom)
           call mpibcast(Model_rlon, Nudge_ncol, mpir8,   0, mpicom)
          !call mpibcast(Model_wgth, Nudge_ncol, mpir8,   0, mpicom)
#endif
          ! Initialize Machine Learning lat/lon info in local arrays
          !------------------------------------------------------
           if (mltbc_regional_on) then
             !info to construct lat-lon for limited region
             mltbc_patch_dx   = 1.0_r8
             mltbc_patch_dy   = 1.0_r8
             mltbc_patch_nlon = int(Nudge_Hwin_lonWidth)
             mltbc_patch_nlat = int(Nudge_Hwin_latWidth)
             mltbc_patch_nxy  = mltbc_patch_nlon * mltbc_patch_nlat
             val1_0 = Nudge_Hwin_lon0 - mltbc_patch_nlon * mltbc_patch_dx / 2.0_r8 + 0.5_r8
             val2_0 = Nudge_Hwin_lon0 + mltbc_patch_nlon * mltbc_patch_dx / 2.0_r8 - 0.5_r8
             val3_0 = Nudge_Hwin_lat0 - mltbc_patch_nlat * mltbc_patch_dy / 2.0_r8 + 0.5_r8
             val4_0 = Nudge_Hwin_lat0 + mltbc_patch_nlat * mltbc_patch_dy / 2.0_r8 - 0.5_r8

             if ( (mltbc_predict_option == 2) .and. (mltbc_patch_nxy /= 9216) ) then 
               write(iulog,*) 'Machine Learning NUDGING: Convolution model 2 only works on 96x96 lat-lon grid'
               write(iulog,*) 'Machine Learning NUDGING: Invalid Nudge_Hwin_latWidth and/or Nudge_Hwin_lonWidth'
               call endrun('Machine Learning NUDGING:: ERROR in namelist for Conv2d setup')
             end if 
             if( (mltbc_predict_option /= 2) .and. mltbc_patch_nxy /= 4900 ) then
               write(iulog,*) 'Machine Learning NUDGING: Convolution model 1 only works on 70x70 lat-lon grid'
               write(iulog,*) 'Machine Learning NUDGING: Invalid Nudge_Hwin_latWidth and/or Nudge_Hwin_lonWidth'
               call endrun('Machine Learning NUDGING:: ERROR in namelist for Conv2d setup')
             end if 

             if((val1_0.lt.0.).or.(val2_0.ge.360.)) then
               write(iulog,*) 'Machine Learning NUDGING: Conv2d Window lon range not in [0,+360]'
               write(iulog,*) 'Machine Learning NUDGING: Conv2d lons,lone=',val1_0,val2_0 
               call endrun('Machine Learning NUDGING:: ERROR in namelist for Conv2d setup')
             endif
             if((val3_0.lt.-90.).or.(val4_0.gt.+90.)) then
               write(iulog,*) 'Machine Learning NUDGING: Conv2d Window lat range not in [-90,+90]'
               write(iulog,*) 'Machine Learning NUDGING: Conv2d lats,late=',val3_0,val4_0
               call endrun('Machine Learning NUDGING:: ERROR in namelist for Conv2d setup')
             endif

             allocate(mltbc_patch_lat(mltbc_patch_nlat),stat=istat)
             call alloc_err(istat,'Machine Learning NUDGING','mltbc_patch_lat',mltbc_patch_nlat)
             allocate(mltbc_patch_lon(mltbc_patch_nlon),stat=istat)
             call alloc_err(istat,'Machine Learning NUDGING','mltbc_patch_lon',mltbc_patch_nlon)
             allocate(mltbc_se2latlon_ind(mltbc_patch_nxy,5),stat=istat)
             call alloc_err(istat,'Machine Learning NUDGING','mltbc_se2latlon_ind',5*mltbc_patch_nxy)
             allocate(mltbc_se2latlon_wgt(mltbc_patch_nxy,5),stat=istat)
             call alloc_err(istat,'Machine Learning NUDGING','mltbc_se2latlon_wgt',5*mltbc_patch_nxy)
             allocate(mltbc_latlon2se_ind(Nudge_ncol,5),stat=istat)
             call alloc_err(istat,'Machine Learning NUDGING','mltbc_latlon2se_ind',5*Nudge_ncol)
             allocate(mltbc_latlon2se_wgt(Nudge_ncol,5),stat=istat)
             call alloc_err(istat,'Machine Learning NUDGING','mltbc_latlon2se_wgt',5*Nudge_ncol)

             if (masterproc) then
               do i = 1, mltbc_patch_nlon
                 mltbc_patch_lon(i) = (val1_0+(i-1)*mltbc_patch_dx)*SHR_CONST_PI/180._r8
               end do
               do j = 1, mltbc_patch_nlat
                 mltbc_patch_lat(j) = (val3_0+(j-1)*mltbc_patch_dy)*SHR_CONST_PI/180._r8
               end do

               call t_startf('mltbc_don_interp_init')
               !initialize the weigthing fuction to interpolate model grid to lat-lon 
               call se2latlon_interp_init(Nudge_ncol, Model_rlon, Model_rlat, & 
                                          mltbc_patch_nlon, mltbc_patch_nlat, & 
                                          mltbc_patch_lon, mltbc_patch_lat,  &
                                          mltbc_se2latlon_ind, mltbc_se2latlon_wgt)  

               !initialize the weigthing fuction to interpolate lat-lon to model grid 
               call latlon2se_interp_init(Nudge_ncol, Model_rlon, Model_rlat, &
                                          mltbc_patch_nlon, mltbc_patch_nlat, &
                                          mltbc_patch_lon, mltbc_patch_lat,  &
                                          mltbc_latlon2se_ind, mltbc_latlon2se_wgt) 

               call t_stopf('mltbc_don_interp_init')

             end if 

#ifdef SPMD
             call mpibcast(mltbc_patch_nlon,   1,                     mpiint,  0, mpicom)
             call mpibcast(mltbc_patch_nlat,   1,                     mpiint,  0, mpicom)
             call mpibcast(mltbc_patch_nxy,    1,                     mpiint,  0, mpicom)
             call mpibcast(mltbc_patch_lon,    mltbc_patch_nlon,  mpir8,   0, mpicom)
             call mpibcast(mltbc_patch_lat,    mltbc_patch_nlat,  mpir8,   0, mpicom)
             call mpibcast(mltbc_se2latlon_ind, 5*mltbc_patch_nxy, mpiint,  0, mpicom)
             call mpibcast(mltbc_se2latlon_wgt, 5*mltbc_patch_nxy, mpir8,   0, mpicom)
             call mpibcast(mltbc_latlon2se_ind, 5*Nudge_ncol,          mpiint,  0, mpicom)
             call mpibcast(mltbc_latlon2se_wgt, 5*Nudge_ncol,          mpir8,   0, mpicom)
#endif
           end if       

           if (masterproc) then
             if (mltbc_regional_on) then
               write(iulog,*) 'Machine Learning NUDGING: regional ML model'
             else 
               write(iulog,*) 'Machine Learning NUDGING: global   ML model'
             end if   
           end if 
      case default
           call endrun('nudging_init error: nudge method should &
                       &be either Step, Linear , IMT or MLTBC...')
   end select        

   ! End Routine
   !------------
   return
  end subroutine ! nudging_init
  !================================================================


  !================================================================
  subroutine nudging_timestep_init(phys_state)
   !
   ! NUDGING_TIMESTEP_INIT:
   !                 Check the current time and update Model/Nudging
   !                 arrays when necessary. Toggle the Nudging flag
   !                 when the time is withing the nudging window.
   !===============================================================
   use physics_types,only: physics_state
   use constituents ,only: cnst_get_ind
   use dycore       ,only: dycore_is
   use ppgrid       ,only: pver,pcols,begchunk,endchunk
   use filenames    ,only: interpret_filename_spec
   use cam_history  ,only: outfld
   use physconst    ,only: cpair, gravit, rga
   use phys_grid    ,only: get_rlat_all_p, get_rlon_all_p

   ! Arguments
   !-----------
   type(physics_state),intent(in):: phys_state(begchunk:endchunk)
  
   ! Local values
   !----------------
   integer Year,Month,Day,Sec
   integer YMD1,YMD2,YMD
   logical Update_Model,Update_Nudge,Sync_Error
   logical After_Beg   ,Before_End
   integer lchnk,ncol,indw, k 
   character(len=2000) err_str

   if (mltbc_nudge) return

   ! Check if Nudging is initialized
   !---------------------------------
   if(.not.Nudge_Initialized) then
     call endrun('nudging_timestep_init:: Nudging NOT Initialized')
   endif

   ! Get Current time
   !--------------------
   call get_curr_date(Year,Month,Day,Sec)
   YMD=(Year*10000) + (Month*100) + Day

   !-------------------------------------------------------
   ! Determine if the current time is AFTER the begining time
   ! and if it is BEFORE the ending time.
   !-------------------------------------------------------
   YMD1=(Nudge_Beg_Year*10000) + (Nudge_Beg_Month*100) + Nudge_Beg_Day
   call timemgr_time_ge(YMD1,Nudge_Beg_Sec,         &
                        YMD ,Sec          ,After_Beg)

   YMD1=(Nudge_End_Year*10000) + (Nudge_End_Month*100) + Nudge_End_Day
   call timemgr_time_ge(YMD ,Sec,                    &
                        YMD1,Nudge_End_Sec,Before_End)

   !--------------------------------------------------------------
   ! When past the NEXT time, Update Model Arrays and time indices
   !--------------------------------------------------------------
   YMD1=(Model_Next_Year*10000) + (Model_Next_Month*100) + Model_Next_Day
   call timemgr_time_ge(YMD1,Model_Next_Sec,            &
                        YMD ,Sec           ,Update_Model)

   if((Before_End).and.(Update_Model)) then
     ! Increment the Model times by the current interval
     !---------------------------------------------------
     Model_Curr_Year =Model_Next_Year
     Model_Curr_Month=Model_Next_Month
     Model_Curr_Day  =Model_Next_Day
     Model_Curr_Sec  =Model_Next_Sec
     YMD1=(Model_Curr_Year*10000) + (Model_Curr_Month*100) + Model_Curr_Day
     call timemgr_time_inc(YMD1,Model_Curr_Sec,              &
                           YMD2,Model_Next_Sec,Model_Step,0,0)

     ! Check for Sync Error where NEXT model time after the update
     ! is before the current time. If so, reset the next model
     ! time to a Model_Step after the current time.
     !--------------------------------------------------------------
     call timemgr_time_ge(YMD2,Model_Next_Sec,            &
                          YMD ,Sec           ,Sync_Error)
     if(Sync_Error) then
       Model_Curr_Year =Year
       Model_Curr_Month=Month
       Model_Curr_Day  =Day
       Model_Curr_Sec  =Sec
       call timemgr_time_inc(YMD ,Model_Curr_Sec,              &
                             YMD2,Model_Next_Sec,Model_Step,0,0)
       write(iulog,*) 'NUDGING: WARNING - Model_Time Sync ERROR... CORRECTED'
     endif
     Model_Next_Year =(YMD2/10000)
     YMD2            = YMD2-(Model_Next_Year*10000)
     Model_Next_Month=(YMD2/100)
     Model_Next_Day  = YMD2-(Model_Next_Month*100)

     if ( .not. Nudge_Loc_PhysOut ) then
        ! Load values at Current into the Model arrays
        !-----------------------------------------------
        call cnst_get_ind('Q',indw)
        do lchnk=begchunk,endchunk
          ncol=phys_state(lchnk)%ncol
          if ( Nudge_Uprof .ne. 0 ) then
             Model_U(:ncol,:pver,lchnk)=phys_state(lchnk)%u(:ncol,:pver)
          end if
          if ( Nudge_Vprof .ne. 0 ) then
             Model_V(:ncol,:pver,lchnk)=phys_state(lchnk)%v(:ncol,:pver)
          end if
          if ( Nudge_Tprof .ne. 0 ) then
             Model_T(:ncol,:pver,lchnk)=phys_state(lchnk)%t(:ncol,:pver)
          end if
          if ( Nudge_Qprof .ne. 0 ) then
             Model_Q(:ncol,:pver,lchnk)=phys_state(lchnk)%q(:ncol,:pver,indw)
          end if
          if ( Nudge_PSprof .ne. 0 ) then
             Model_PS(:ncol,lchnk)=phys_state(lchnk)%ps(:ncol)
          end if
        end do
     end if
   end if

   !----------------------------------------------------------------
   ! When past the NEXT time, Update Nudging Arrays and time indices
   !----------------------------------------------------------------
   YMD1=(Nudge_Next_Year*10000) + (Nudge_Next_Month*100) + Nudge_Next_Day
   call timemgr_time_ge(YMD1,Nudge_Next_Sec,            &
                        YMD ,Sec           ,Update_Nudge)

   ! update the logical variables for calculation
   ! of nudging tendency when Nudge_Loc_PhysOut is true
   if ( Nudge_Loc_PhysOut ) then
      l_Update_Model = Update_Model
      l_Update_Nudge = Update_Nudge
         l_After_Beg = After_Beg
        l_Before_End = Before_End
   end if

   if((Before_End).and.(Update_Nudge)) then
     ! Increment the Nudge times by the current interval
     !---------------------------------------------------
     Nudge_Curr_Year =Nudge_Next_Year
     Nudge_Curr_Month=Nudge_Next_Month
     Nudge_Curr_Day  =Nudge_Next_Day
     Nudge_Curr_Sec  =Nudge_Next_Sec
     YMD1=(Nudge_Curr_Year*10000) + (Nudge_Curr_Month*100) + Nudge_Curr_Day
     call timemgr_time_inc(YMD1,Nudge_Curr_Sec,              &
                           YMD2,Nudge_Next_Sec,Nudge_Step,0,0)
     Nudge_Next_Year =(YMD2/10000)
     YMD2            = YMD2-(Nudge_Next_Year*10000)
     Nudge_Next_Month=(YMD2/100)
     Nudge_Next_Day  = YMD2-(Nudge_Next_Month*100)

     ! Update the Nudge arrays with analysis
     ! data at the NEXT time
     !-----------------------------------------------
     Nudge_File=interpret_filename_spec(Nudge_File_Template      , &
                                         yr_spec=Nudge_Next_Year , &
                                        mon_spec=Nudge_Next_Month, &
                                        day_spec=Nudge_Next_Day  , &
                                        sec_spec=Nudge_Next_Sec    )

     Nudge_SRF_File=interpret_filename_spec(Nudge_SRF_File_Template  , &
                                            yr_spec=Nudge_Next_Year  , &
                                            mon_spec=Nudge_Next_Month, &
                                            day_spec=Nudge_Next_Day  , &
                                            sec_spec=Nudge_Next_Sec    )  

     if(masterproc) then
       write(iulog,*) 'NUDGING: Reading analyses:',trim(Nudge_Path)//trim(Nudge_File)
       if(Nudge_Land) then 
         write(iulog,*) 'NUDGING: Reading surface analyses:',trim(Nudge_Path)//trim(Nudge_SRF_File)
       end if 
     endif
     ! How to manage MISSING values when there are 'Gaps' in the analyses data?
     !  Check for analyses file existence. If it is there, then read data.
     !  If it is not, then issue a warning and switch off nudging to 'coast'
     !  thru the gap.
     !------------------------------------------------------------------------
     if(dycore_is('UNSTRUCTURED')) then 
       ! read nudging data from reanalysis or user specified file 
       call nudging_update_analyses_se(trim(Nudge_Path)//trim(Nudge_File))
  
       if(Nudge_Land) then 
         call nudging_update_srf_analyses_se(trim(Nudge_Path)//trim(Nudge_SRF_File))
       end if
     elseif(dycore_is('EUL')) then
       call nudging_update_analyses_eul(trim(Nudge_Path)//trim(Nudge_File))
     else !if(dycore_is('LR')) then
       call nudging_update_analyses_fv(trim(Nudge_Path)//trim(Nudge_File))
     endif

     if ( .not. dycore_is('EUL') ) then
       select case (Nudge_Method)
          case ('Linear')
             call t_startf ('nudging_interp')
             call linear_interpolation   (INTP_U, Target_U)
             call linear_interpolation   (INTP_V, Target_V)
             call linear_interpolation   (INTP_Q, Target_Q)
             call linear_interpolation   (INTP_T, Target_T)
             call linear_interpolation_2d(INTP_PS, Target_PS)
             call linear_interpolation_2d(INTP_PHIS, Target_PHIS)
             if (Nudge_Land) then
               call linear_interpolation_2d(INTP_U10,    Target_U10)
               call linear_interpolation_2d(INTP_V10,    Target_V10)
               call linear_interpolation_2d(INTP_T2,     Target_T2)
               call linear_interpolation_2d(INTP_TD2,    Target_TD2)
               call linear_interpolation_2d(INTP_TS,     Target_TS)
               call linear_interpolation_2d(INTP_PRECC,  Target_PRECC)
               call linear_interpolation_2d(INTP_PRECL,  Target_PRECL)
               call linear_interpolation_2d(INTP_PRECSC, Target_PRECSC)
               call linear_interpolation_2d(INTP_PRECSL, Target_PRECSL)
               call linear_interpolation_2d(INTP_EVAP,   Target_EVAP)
               call linear_interpolation_2d(INTP_SHFLX,  Target_SHFLX)
               call linear_interpolation_2d(INTP_LHFLX,  Target_LHFLX)
               call linear_interpolation_2d(INTP_FSNS,   Target_NETSW)
               call linear_interpolation_2d(INTP_FLDS,   Target_FLWDS)
               call linear_interpolation_2d(INTP_FSDS,   Target_FSDS)
               call linear_interpolation_2d(INTP_FSDSD,  Target_FSDSD)
               call linear_interpolation_2d(INTP_FSDSUV, Target_FSDSUV)
             end if ! Nudge_Land
             call t_stopf ('nudging_interp')
          case default
             ! No interpolation is needed for Step or IMT 
       end select
     end if
   end if 
   !-------------------------------------------------------
   ! Toggle Nudging flag when the time interval is between
   ! beginning and ending times, and the analyses file exists.
   !-------------------------------------------------------
   if((After_Beg).and.(Before_End)) then
     if(Nudge_File_Present) then
       Nudge_ON=.true.
     else
       Nudge_ON=.false.
       if(Nudge_Allow_Missing_File) then
         if(masterproc) then
           write(iulog,*) 'NUDGING: WARNING - Nudging data file NOT FOUND. Switching '
           write(iulog,*) 'NUDGING:           nudging OFF to coast thru the gap. '
         endif   
       else
         write(err_str,*) 'NUDGING: Nudging data file (',trim(adjustl(Nudge_File)),') NOT FOUND; ', errmsg(__FILE__, __LINE__)
          call endrun(err_str)
       endif 
     endif
     if(Nudge_Land) then
       !For land surfac nudging file 
       if(Nudge_SRF_File_Present) then
         Nudge_SRF_On=.true.
       else
         Nudge_SRF_On=.false.
         if(Nudge_Allow_Missing_File) then
           if(masterproc) then
             write(iulog,*) 'NUDGING: WARNING - Nudging data file NOT FOUND. Switching '
             write(iulog,*) 'NUDGING:           nudging OFF to coast thru the gap. '
           endif
         else
           write(err_str,*) 'NUDGING: Nudging data file (',trim(adjustl(Nudge_SRF_File)),') NOT FOUND; ', errmsg(__FILE__, __LINE__)
           call endrun(err_str)
         endif
       endif
     end if 
   else
     Nudge_ON=.false.
     Nudge_SRF_On=.false.
   endif

   !---------------------------------------------------
   ! If Data arrays have changed update stepping arrays
   !---------------------------------------------------
   if ( .not. Nudge_Loc_PhysOut ) then
      if ((Before_End).and.((Update_Nudge).or.(Update_Model))) then
        do lchnk=begchunk,endchunk
          ncol=phys_state(lchnk)%ncol
          if (Nudge_Uprof .ne. 0) then
            Nudge_Ustep(:ncol,:pver,lchnk)=(  Target_U(:ncol,:pver,lchnk)  &
                                              -Model_U(:ncol,:pver,lchnk)) &
                                           *Nudge_Utau(:ncol,:pver,lchnk)
          end if
          if (Nudge_Vprof .ne. 0) then
            Nudge_Vstep(:ncol,:pver,lchnk)=(  Target_V(:ncol,:pver,lchnk)  &
                                              -Model_V(:ncol,:pver,lchnk)) &
                                           *Nudge_Vtau(:ncol,:pver,lchnk)
          end if
          if (Nudge_Tprof .ne. 0) then
            Nudge_Tstep(:ncol,:pver,lchnk)=(  Target_T(:ncol,:pver,lchnk)  &
                                              -Model_T(:ncol,:pver,lchnk)) &
                                           *Nudge_Ttau(:ncol,:pver,lchnk)
          end if
          if (Nudge_Qprof .ne. 0) then
            Nudge_Qstep(:ncol,:pver,lchnk)=(  Target_Q(:ncol,:pver,lchnk)  &
                                              -Model_Q(:ncol,:pver,lchnk)) &
                                           *Nudge_Qtau(:ncol,:pver,lchnk)
          end if
          if (Nudge_PSprof .ne. 0) then
            Nudge_PSstep(:ncol,lchnk)=(  Target_PS(:ncol,lchnk)  &
                                              -Model_PS(:ncol,lchnk)) &
                                           *Nudge_PStau(:ncol,lchnk)
          end if
        end do

         !******************
         ! DIAG
         !******************
!        if(masterproc) then
!         write(iulog,*) 'PFC: Target_T(1,:pver,begchunk)=',Target_T(1,:pver,begchunk)  
!         write(iulog,*) 'PFC:  Model_T(1,:pver,begchunk)=',Model_T(1,:pver,begchunk)
!         write(iulog,*) 'PFC: Nudge_Tstep(1,:pver,begchunk)=',Nudge_Tstep(1,:pver,begchunk)
!         write(iulog,*) 'PFC: Nudge_Xstep arrays updated:'
!        endif

      else

         do lchnk=begchunk,endchunk
            ncol=phys_state(lchnk)%ncol
            ! The following lines are used to reset the nudging tendency
            ! to zero in order to perform an intermittent simulation
            if (Nudge_Uprof .ne. 0) then
                Nudge_Ustep(:ncol,:pver,lchnk) = 0._r8
            end if
            if (Nudge_Vprof .ne. 0) then
                Nudge_Vstep(:ncol,:pver,lchnk) = 0._r8
            end if
            if (Nudge_Tprof .ne. 0) then
                Nudge_Tstep(:ncol,:pver,lchnk) = 0._r8
            end if
            if (Nudge_Qprof .ne. 0) then
                Nudge_Qstep(:ncol,:pver,lchnk) = 0._r8
            end if
            if (Nudge_PSprof .ne. 0) then
                Nudge_PSstep(:ncol,lchnk)      = 0._r8
            end if
         end do

      end if

   end if   ! if for Nudge_Loc_PhysOut

   ! End Routine
   !------------
   return
  end subroutine ! nudging_timestep_init
  !================================================================

  !================================================================
  subroutine mltbc_timestep_init(state,dtime)
   !
   ! DEEPONET_TIMESTEP_INIT:
   !                 Check the current time and update Model/Nudging
   !                 arrays when necessary. Toggle the Nudging flag
   !                 when the time is withing the nudging window.
   !===============================================================
   use physics_types,only: physics_state
   use constituents, only: cnst_get_ind
   use dycore,       only: dycore_is
   use ppgrid,       only: pver,pverp,pcols,begchunk,endchunk 
   use phys_grid,    only: get_ncols_p, get_rlat_all_p, get_rlon_all_p, get_lon_all_p, &
                           get_lat_all_p, gather_chunk_to_field,scatter_field_to_chunk
   use dyn_grid,     only: get_dyn_grid_parm,get_horiz_grid_dim_d, &
                           get_horiz_grid_d, get_dyn_grid_parm_real1d
   use filenames,    only: interpret_filename_spec
   use time_manager, only: get_nstep
   use cam_history,  only: outfld
   use physconst,    only: cpair, gravit, rga
   use phys_grid,    only: get_rlat_all_p, get_rlon_all_p
   use parallel_mod, only: par

   ! Arguments
   !-----------
   type(physics_state),  intent(in) :: state(begchunk:endchunk)
   real(r8),             intent(in) :: dtime

   ! Local values
   !----------------
   integer Year,Month,Day,Sec
   integer YMD1,YMD2,YMD
   logical Update_Model,Update_Nudge,Sync_Error
   logical After_Beg   ,Before_End
   integer lchnk,ncol,i,j,k,n,m,indw
   integer nstep ! current timestep number
   character(len=2000) err_str

   !temporary working arrays 
   real(r8):: ftem(pcols,pver)
   real(r8):: ftem2(pcols) ! temporary workspace
   real(r8):: arr(pcols,begchunk:endchunk,pver)
   real(r8):: tmp_tend(pcols,1,begchunk:endchunk)

   nstep = get_nstep()

   ! Check if Nudging is initialized
   !---------------------------------
   if(.not.Nudge_Initialized) then
     call endrun('mltbc_timestep_init:: Nudging NOT Initialized')
   endif

   ! Get Current time
   !--------------------
   call get_curr_date(Year,Month,Day,Sec)
   YMD=(Year*10000) + (Month*100) + Day

   !-------------------------------------------------------
   ! Determine if the current time is AFTER the begining time
   ! and if it is BEFORE the ending time.
   !-------------------------------------------------------
   YMD1=(Nudge_Beg_Year*10000) + (Nudge_Beg_Month*100) + Nudge_Beg_Day
   call timemgr_time_ge(YMD1,Nudge_Beg_Sec,         &
                        YMD ,Sec          ,After_Beg)

   YMD1=(Nudge_End_Year*10000) + (Nudge_End_Month*100) + Nudge_End_Day
   call timemgr_time_ge(YMD ,Sec,                    &
                        YMD1,Nudge_End_Sec,Before_End)

   !--------------------------------------------------------------
   ! When past the NEXT time, Update Model Arrays and time indices
   !--------------------------------------------------------------
   YMD1=(Model_Next_Year*10000) + (Model_Next_Month*100) + Model_Next_Day
   call timemgr_time_ge(YMD1,Model_Next_Sec,            &
                        YMD ,Sec           ,Update_Model)

   if((Before_End).and.(Update_Model)) then
     ! Increment the Model times by the current interval
     !---------------------------------------------------
     Model_Curr_Year =Model_Next_Year
     Model_Curr_Month=Model_Next_Month
     Model_Curr_Day  =Model_Next_Day
     Model_Curr_Sec  =Model_Next_Sec
     YMD1=(Model_Curr_Year*10000) + (Model_Curr_Month*100) + Model_Curr_Day
     call timemgr_time_inc(YMD1,Model_Curr_Sec,              &
                           YMD2,Model_Next_Sec,Model_Step,0,0)

     ! Check for Sync Error where NEXT model time after the update
     ! is before the current time. If so, reset the next model
     ! time to a Model_Step after the current time.
     !--------------------------------------------------------------
     call timemgr_time_ge(YMD2,Model_Next_Sec,            &
                          YMD ,Sec           ,Sync_Error)
     if(Sync_Error) then
       Model_Curr_Year =Year
       Model_Curr_Month=Month
       Model_Curr_Day  =Day
       Model_Curr_Sec  =Sec
       call timemgr_time_inc(YMD ,Model_Curr_Sec,              &
                             YMD2,Model_Next_Sec,Model_Step,0,0)
       write(iulog,*) 'NUDGING: WARNING - Model_Time Sync ERROR... CORRECTED'
     endif
     Model_Next_Year =(YMD2/10000)
     YMD2            = YMD2-(Model_Next_Year*10000)
     Model_Next_Month=(YMD2/100)
     Model_Next_Day  = YMD2-(Model_Next_Month*100)

     call cnst_get_ind('Q',indw)
     do lchnk=begchunk,endchunk

       ncol=state(lchnk)%ncol

       if ( Nudge_Uprof .ne. 0 ) then
          Model_U(:ncol,:pver,lchnk)=state(lchnk)%u(:ncol,:pver)
       end if

       if ( Nudge_Vprof .ne. 0 ) then
          Model_V(:ncol,:pver,lchnk)=state(lchnk)%v(:ncol,:pver)
       end if

       if ( Nudge_Tprof .ne. 0 ) then
          Model_T(:ncol,:pver,lchnk)=state(lchnk)%t(:ncol,:pver)
       end if

       if ( Nudge_Qprof .ne. 0 ) then
          Model_Q(:ncol,:pver,lchnk)=state(lchnk)%q(:ncol,:pver,indw)
       end if

       if ( Nudge_PSprof .ne. 0 ) then
          Model_PS(:ncol,lchnk)=state(lchnk)%ps(:ncol)
       end if

     end do
   end if 

   !----------------------------------------------------------------
   ! When past the NEXT time, Update Nudging Arrays and time indices
   !----------------------------------------------------------------
   YMD1=(Nudge_Next_Year*10000) + (Nudge_Next_Month*100) + Nudge_Next_Day
   call timemgr_time_ge(YMD1,Nudge_Next_Sec,            &
                        YMD ,Sec           ,Update_Nudge)


   !-------------------------------------------------------
   ! Toggle Nudging flag when the time interval is between
   ! beginning and ending times, and the analyses file exists.
   !-------------------------------------------------------
   if((After_Beg).and.(Before_End)) then
     Nudge_ON=.true.
     if(Nudge_Land) then
       Nudge_SRF_On=.true.
     end if
     mltbc_nudge=.true.
   else
     Nudge_ON=.false.
     Nudge_SRF_On=.false.
     mltbc_nudge=.false.
   end if

   Nudge_Ustep(:,:,:) = 0._r8
   Nudge_Vstep(:,:,:) = 0._r8
   Nudge_Tstep(:,:,:) = 0._r8
   Nudge_Qstep(:,:,:) = 0._r8
   Nudge_PSstep(:,:)  = 0._r8
   Model_Var(:,:)     = fillvalue
   tmp_tend(:,1,:)    = 0._r8

   if ((nstep > 0).and.(Before_End).and.((Update_Nudge).or.(Update_Model))) then
     !call deepONet to predict nudging tendency 
     if (Nudge_Uprof .ne. 0) then
       do lchnk=begchunk,endchunk
         ncol = state(lchnk)%ncol
         arr(:ncol,lchnk,:pver) = state(lchnk)%u(:ncol,:pver)
       end do 
       call mltbc_gather_data(arr,pver,Nudge_ncol,Model_Var)
       call mltbc_advance(mltbc_model_path,'U',pver,Nudge_ncol,Model_Var,Nudge_Ustep) 
     end if
     if (Nudge_Vprof .ne. 0) then
       do lchnk=begchunk,endchunk
         ncol = state(lchnk)%ncol
         arr(:ncol,lchnk,:pver) = state(lchnk)%v(:ncol,:pver)
       end do
       call mltbc_gather_data(arr,pver,Nudge_ncol,Model_Var)
       call mltbc_advance(mltbc_model_path,'V',pver,Nudge_ncol,Model_Var,Nudge_Vstep) 
     end if
     if (Nudge_Tprof .ne. 0) then
       do lchnk=begchunk,endchunk
         ncol = state(lchnk)%ncol
         arr(:ncol,lchnk,:pver) = state(lchnk)%t(:ncol,:pver)
       end do
       call mltbc_gather_data(arr,pver,Nudge_ncol,Model_Var)
       call mltbc_advance(mltbc_model_path,'T',pver,Nudge_ncol,Model_Var,Nudge_Tstep) 
     end if
     if (Nudge_Qprof .ne. 0) then
       do lchnk=begchunk,endchunk
         ncol = state(lchnk)%ncol
         arr(:ncol,lchnk,:pver) = state(lchnk)%q(:ncol,:pver,indw)
       end do
       call mltbc_gather_data(arr,pver,Nudge_ncol,Model_Var)
       call mltbc_advance(mltbc_model_path,'Q',Nudge_nlev,Nudge_ncol,Model_Var,Nudge_Qstep) 
     end if
     if (Nudge_PSprof .ne. 0) then
       do lchnk=begchunk,endchunk
         ncol = state(lchnk)%ncol
         arr(:ncol,lchnk,1) = state(lchnk)%ps(:ncol)
       end do
       call mltbc_gather_data(arr,1,Nudge_ncol,Model_Var(:,1))
       call mltbc_advance(mltbc_model_path,'PS',1,Nudge_ncol,Model_Var(:,1),tmp_tend)
       Nudge_PSstep(:,:) = tmp_tend(:,1,:)
     end if
   end if 

   !Add a linear relexation of the nudging tendency on the upper layer 
   if (Nudge_Lin_Relax_On) then
     do lchnk=begchunk,endchunk
       ncol = state(lchnk)%ncol
       do i = 1, ncol 
         do k = pver, 1, -1     
           if ( state(lchnk)%pmid(i,k) < p_uv_relax ) then
             Nudge_Utau(i,k,lchnk) = Nudge_Utau(i,k,lchnk) * max(0.01_r8, state(lchnk)%pmid(i,k)/p_uv_relax)
             Nudge_Vtau(i,k,lchnk) = Nudge_Vtau(i,k,lchnk) * max(0.01_r8, state(lchnk)%pmid(i,k)/p_uv_relax)
           end if                    
           if ( state(lchnk)%pmid(i,k) < p_T_relax ) then
             Nudge_Ttau(i,k,lchnk) = Nudge_Ttau(i,k,lchnk) * max(0.01_r8, state(lchnk)%pmid(i,k)/p_T_relax)
           end if                    
           if ( state(lchnk)%pmid(i,k) < p_q_relax ) then
             Nudge_Qtau(i,k,lchnk) = Nudge_Qtau(i,k,lchnk) * max(0.01_r8, state(lchnk)%pmid(i,k)/p_Q_relax)
           end if                    
           if (state(lchnk)%pmid(i,k) < p_norelax ) then
             Nudge_Utau(i,k,lchnk) = 0._r8
             Nudge_Vtau(i,k,lchnk) = 0._r8
             Nudge_Ttau(i,k,lchnk) = 0._r8
             Nudge_Qtau(i,k,lchnk) = 0._r8
           end if
         end do 
       end do
     end do
   end if

   !apply the coefficient 
   do lchnk=begchunk,endchunk

     ncol = state(lchnk)%ncol

     Nudge_Ustep(:ncol,:pver,lchnk) = Nudge_Ustep(:ncol,:pver,lchnk) &
                                       * Nudge_Utau(:ncol,:pver,lchnk)

     Nudge_Vstep(:ncol,:pver,lchnk) = Nudge_Vstep(:ncol,:pver,lchnk) &
                                       * Nudge_Vtau(:ncol,:pver,lchnk)

     Nudge_Tstep(:ncol,:pver,lchnk) = Nudge_Tstep(:ncol,:pver,lchnk) &
                                       * Nudge_Ttau(:ncol,:pver,lchnk)

     Nudge_Qstep(:ncol,:pver,lchnk) = Nudge_Qstep(:ncol,:pver,lchnk) &
                                       * Nudge_Qtau(:ncol,:pver,lchnk)

     Nudge_PSstep(:ncol,lchnk)      = Nudge_PSstep(:ncol,lchnk) &
                                       * Nudge_PStau(:ncol,lchnk)

   end do

   !output Diagnostics 
   do lchnk=begchunk,endchunk

     call outfld('U_bf_ndg',  Model_U(:,:,lchnk), pcols,lchnk)
     call outfld('V_bf_ndg',  Model_V(:,:,lchnk), pcols,lchnk)
     call outfld('T_bf_ndg',  Model_T(:,:,lchnk), pcols,lchnk)
     call outfld('Q_bf_ndg',  Model_Q(:,:,lchnk), pcols,lchnk)
     call outfld('PS_bf_ndg', Model_PS(:,lchnk),  pcols,lchnk)

     call outfld('Nudge_U',   Nudge_Ustep(:,:,lchnk),pcols,lchnk)
     call outfld('Nudge_V',   Nudge_Vstep(:,:,lchnk),pcols,lchnk)
     call outfld('Nudge_T',   Nudge_Tstep(:,:,lchnk),pcols,lchnk)
     call outfld('Nudge_Q',   Nudge_Qstep(:,:,lchnk),pcols,lchnk)
     call outfld('Nudge_PS',  Nudge_PSstep(:,lchnk), pcols,lchnk)

     ncol=state(lchnk)%ncol
     ftem(:ncol,:pver) = Nudge_Ustep(:ncol,:pver,lchnk)*(state(lchnk)%pdel(:ncol,:pver)/gravit)
     ftem2(1:ncol)     = sum( ftem(1:ncol,:), 2 )
     call outfld('Nudge_U_vint',ftem2,pcols,lchnk)
     ftem(:ncol,:pver) = Nudge_Vstep(:ncol,:pver,lchnk)*(state(lchnk)%pdel(:ncol,:pver)/gravit)
     ftem2(1:ncol)     = sum( ftem(1:ncol,:), 2 )
     call outfld('Nudge_V_vint',ftem2,pcols,lchnk)
     ftem(:ncol,:pver) = Nudge_Tstep(:ncol,:pver,lchnk)*(state(lchnk)%pdel(:ncol,:pver)/gravit) 
     ftem2(1:ncol)     = sum( ftem(1:ncol,:), 2 )
     call outfld('Nudge_T_vint',ftem2,pcols,lchnk)
     ftem(:ncol,:pver) = Nudge_Qstep(:ncol,:pver,lchnk)*(state(lchnk)%pdel(:ncol,:pver)/gravit) 
     ftem2(1:ncol)     = sum( ftem(1:ncol,:), 2 )
     call outfld('Nudge_Q_vint',ftem2,pcols,lchnk)

   end do 

   ! End Routine
   !------------
   return
  end subroutine ! mltbc_timestep_init

  !===========================================================================
  ! JS - 11/05/2019: Based on Shixuan Zhang's suggestion, the calculation of 
  !                  nudging tendency can be done at the same location where 
  !                  the model state variables are written out.
  !                  This subroutine is called by tphysbc and works on a chunk 
  !===========================================================================
  subroutine nudging_calc_tend(state,pbuf,dtime)
   use physics_types,only: physics_state
   use constituents ,only: cnst_get_ind,pcnst
   use ppgrid       ,only: pver,pcols
   use cam_history  ,only: outfld
   use physconst    ,only: cpair, gravit, rga
   use physics_buffer, only: pbuf_get_index, pbuf_old_tim_idx, pbuf_get_field, &
                             pbuf_set_field, physics_buffer_desc
   use phys_grid     , only: get_rlat_all_p, get_rlon_all_p

   ! Arguments
   !-----------
   type(physics_state),intent(in):: state
   type(physics_buffer_desc), pointer :: pbuf(:)
   real(r8), intent(in):: dtime

   ! Local values
   !----------------
   integer :: pblh_idx        ! PBL pbuf
   integer :: lchnk,ncol,indw
   integer :: i,k,m
   real(r8):: ftem(pcols,pver), ftem2(pcols) ! temporary workspace
   real(r8):: rlats(pcols)
   real(r8):: rlons(pcols)

   real(r8):: u_tend(pcols,pver)
   real(r8):: v_tend(pcols,pver)
   real(r8):: t_tend(pcols,pver)
   real(r8):: q_tend(pcols,pver)
   real(r8):: ps_tend(pcols)

   real(r8):: zm_mod(pcols,pver)
   real(r8):: zm_obs(pcols,pver)

   real(r8), pointer :: PBLH(:)                ! Planetary boundary layer height

   if (mltbc_nudge) return

   lchnk = state%lchnk
   ncol  = state%ncol

   pblh_idx  = pbuf_get_index('pblh')
   call pbuf_get_field(pbuf, pblh_idx,  PBLH)

   ! Load values at Current into the Model arrays
   !-----------------------------------------------
   call cnst_get_ind('Q',indw)
   call get_rlon_all_p(lchnk,pcols,rlons)
   call get_rlat_all_p(lchnk,pcols,rlats)

   Model_U(:ncol,:pver,lchnk)=state%u(:ncol,:pver)
   Model_V(:ncol,:pver,lchnk)=state%v(:ncol,:pver)
   Model_T(:ncol,:pver,lchnk)=state%t(:ncol,:pver)
   Model_Q(:ncol,:pver,lchnk)=state%q(:ncol,:pver,indw)
   Model_PS(:ncol,lchnk)=state%ps(:ncol)
   Model_PHIS(:ncol,lchnk)=state%phis(:ncol)

   if ((l_Before_End).and.((l_Update_Nudge).or.(l_Update_Model))) then
     ! Update the nudging tendency
     ps_tend(:ncol)             = 0._r8
     u_tend(:ncol,:pver)        = 0._r8
     v_tend(:ncol,:pver)        = 0._r8
     t_tend(:ncol,:pver)        = 0._r8
     q_tend(:ncol,:pver)        = 0._r8
     call update_nudging_tend (ncol, dtime, Model_PS(:,lchnk), Target_PS(:,lchnk), Nudge_PStau(:,lchnk),   & !In 
                               Model_U(:,:,lchnk), Target_U(:,:,lchnk), Nudge_Utau(:,:,lchnk),             & !In 
                               Model_V(:,:,lchnk), Target_V(:,:,lchnk), Nudge_Vtau(:,:,lchnk),             & !In 
                               Model_T(:,:,lchnk), Target_T(:,:,lchnk), Nudge_Ttau(:,:,lchnk),             & !In  
                               Model_Q(:,:,lchnk), Target_Q(:,:,lchnk), Nudge_Qtau(:,:,lchnk),             & !In 
                               Target_U10(:,lchnk),Target_V10(:,lchnk), Target_T2(:,lchnk),                & !In 
                               Target_TD2(:,lchnk),Target_Q2(:,lchnk), Nudge_SRFtau(:,lchnk),              & !In  
                               Nudge_SRF_On, Nudge_SRF_State_On, Nudge_SRF_Q_On,                           & !In
                               Model_PHIS(:,lchnk), Target_PHIS(:,lchnk), PBLH,                            & !In
                               Nudge_PS_Adjust_On, Nudge_Q_Adjust_On, Nudge_Pdep_Weight_On,                & !In
                               Nudge_Lin_Relax_On, Nudge_NO_PBL_UV, Nudge_NO_PBL_T, Nudge_NO_PBL_Q,        & !In 
                               Nudge_PSprof, Nudge_PS_OPT, Nudge_Uprof, Nudge_Vprof, Nudge_UV_OPT,         & !In 
                               Nudge_Tprof,  Nudge_T_OPT, Nudge_Qprof, Nudge_Q_OPT,                        & !In 
                               zm_obs, zm_mod, ps_tend, u_tend, v_tend, t_tend, q_tend )                     !Out

     Nudge_PSstep(:ncol,lchnk)      = ps_tend(:ncol)
     Nudge_Ustep(:ncol,:pver,lchnk) = u_tend(:ncol,:pver)
     Nudge_Vstep(:ncol,:pver,lchnk) = v_tend(:ncol,:pver)
     Nudge_Tstep(:ncol,:pver,lchnk) = t_tend(:ncol,:pver)
     Nudge_Qstep(:ncol,:pver,lchnk) = q_tend(:ncol,:pver)
   else
     ! The following lines are used to reset the nudging tendency
     ! to zero in order to perform an intermittent simulation
     Nudge_Ustep(:ncol,:pver,lchnk) = 0._r8
     Nudge_Vstep(:ncol,:pver,lchnk) = 0._r8
     Nudge_Tstep(:ncol,:pver,lchnk) = 0._r8
     Nudge_Qstep(:ncol,:pver,lchnk) = 0._r8
     Nudge_PSstep(:ncol,lchnk)      = 0._r8
   end if         ! update nudging tendency

   call outfld('U_bf_ndg',  Model_U(:,:,lchnk), pcols, lchnk)
   call outfld('V_bf_ndg',  Model_V(:,:,lchnk), pcols, lchnk)
   call outfld('T_bf_ndg',  Model_T(:,:,lchnk), pcols, lchnk)
   call outfld('Q_bf_ndg',  Model_Q(:,:,lchnk), pcols, lchnk)
   call outfld('PS_bf_ndg', Model_PS(:,lchnk),  pcols, lchnk)

   do k = 1, pver
     ftem(:ncol,k) = zm_mod(:ncol,k) + state%phis(:ncol)*rga
   end do
   call outfld('Z3_bf_ndg', ftem,  pcols,lchnk)

   call outfld('U_ref',  Target_U(:,:,lchnk), pcols, lchnk)
   call outfld('V_ref',  Target_V(:,:,lchnk), pcols, lchnk)
   call outfld('T_ref',  Target_T(:,:,lchnk), pcols, lchnk)
   call outfld('Q_ref',  Target_Q(:,:,lchnk), pcols, lchnk)
   call outfld('PS_ref', Target_PS(:,lchnk),  pcols, lchnk)

   ftem2(:ncol) = Target_PHIS(:ncol,lchnk)
   call outfld('PHIS_ref', ftem2, pcols,lchnk)

   do k = 1, pver
     ftem(:ncol,k) = zm_obs(:ncol,k) + Target_PHIS(:ncol,lchnk)*rga
   end do
   call outfld('Z3_ref', ftem,  pcols,lchnk)

   call outfld('Nudge_PS', Nudge_PSstep(:,lchnk),  pcols, lchnk)
   call outfld('Nudge_U',  Nudge_Ustep(:,:,lchnk), pcols, lchnk)
   call outfld('Nudge_V',  Nudge_Vstep(:,:,lchnk), pcols, lchnk)

   ftem(:ncol,:pver) = Nudge_Ustep(:ncol,:pver,lchnk)*(state%pdel(:ncol,:pver)/gravit)
   ftem2(1:ncol)     = sum( ftem(1:ncol,:), 2 )
   call outfld('Nudge_U_vint',ftem2,pcols,lchnk)

   ftem(:ncol,:pver) = Nudge_Vstep(:ncol,:pver,lchnk)*(state%pdel(:ncol,:pver)/gravit)
   ftem2(1:ncol)     = sum( ftem(1:ncol,:), 2 )
   call outfld('Nudge_V_vint',ftem2,pcols,lchnk)

   if((Nudge_Tprof.eq.0).and.(Nudge_Uprof.ne.0)) then 
     ftem(:ncol,:pver) = ( Target_T(:ncol,:pver,lchnk)  &
                           -Model_T(:ncol,:pver,lchnk)) &
                         *Nudge_Utau(:ncol,:pver,lchnk)
     call outfld('Nudge_T',ftem,pcols,lchnk)

     ftem(:ncol,:pver) = ftem(:ncol,:pver)*(state%pdel(:ncol,:pver)/gravit)
     ftem2(1:ncol)     = sum( ftem(1:ncol,:), 2 )
     call outfld('Nudge_T_vint',ftem2,pcols,lchnk)
   else 
     call outfld('Nudge_T', Nudge_Tstep(:,:,lchnk), pcols, lchnk)

     ftem(:ncol,:pver) = Nudge_Tstep(:ncol,:pver,lchnk)*(state%pdel(:ncol,:pver)/gravit) 
     ftem2(1:ncol)     = sum( ftem(1:ncol,:), 2 )
     call outfld('Nudge_T_vint',ftem2,pcols,lchnk)
   end if 
   if((Nudge_Qprof.eq.0).and.(Nudge_Uprof.ne.0)) then
     ftem(:ncol,:pver) = ( Target_Q(:ncol,:pver,lchnk)  &
                           -Model_Q(:ncol,:pver,lchnk)) &
                         *Nudge_Utau(:ncol,:pver,lchnk)
     call outfld('Nudge_Q', ftem,pcols,lchnk)

     ftem(:ncol,:pver) = ftem(:ncol,:pver)*(state%pdel(:ncol,:pver)/gravit)
     ftem2(1:ncol)     = sum( ftem(1:ncol,:), 2 )
     call outfld('Nudge_Q_vint',ftem2,pcols,lchnk)
   else
     call outfld('Nudge_Q', Nudge_Qstep(:,:,lchnk), pcols, lchnk)

     ftem(:ncol,:pver) = Nudge_Qstep(:ncol,:pver,lchnk)*(state%pdel(:ncol,:pver)/gravit)
     ftem2(1:ncol)     = sum( ftem(1:ncol,:), 2 )
     call outfld('Nudge_Q_vint',ftem2,pcols,lchnk)
   end if

   return
  end subroutine  ! nudging_calc_tend

  !================================================================
  subroutine nudging_update_land_surface(state, pbuf, cam_in, cam_out, dtime)
   use physics_types,  only: physics_state
   use camsrfexch   ,  only: cam_out_t, cam_in_t
   use ppgrid       ,  only: pver,pcols,pverp
   use cam_history  ,  only: outfld
   use physconst    ,  only: cpair, gravit, rair, rga, stebol, zvir
   use physics_buffer, only: pbuf_get_index, pbuf_old_tim_idx, pbuf_get_field, &
                             pbuf_set_field, physics_buffer_desc
   use geopotential,   only: geopotential_t
   use comsrf,         only: prcsnw
   use time_manager,   only: get_nstep

   ! Arguments
   !-----------
   type(physics_buffer_desc), pointer :: pbuf(:)
   type(physics_state), intent(in)    :: state
   type(cam_out_t),     intent(inout) :: cam_out
   type(cam_in_t),      intent(in)    :: cam_in
   real(r8),            intent(in)    :: dtime

   !convective precipitation variables
   real(r8), pointer :: prec_dp(:)     ! total precipitation   from ZM convection
   real(r8), pointer :: snow_dp(:)     ! snow from ZM   convection
   real(r8), pointer :: prec_sh(:)     ! total precipitation   from Hack convection
   real(r8), pointer :: snow_sh(:)     ! snow from   Hack convection
   real(r8), pointer :: prec_sed(:)    ! total precipitation   from ZM convection
   real(r8), pointer :: snow_sed(:)    ! snow from ZM   convection
   real(r8), pointer :: prec_pcw(:)    ! total precipitation   from Hack convection
   real(r8), pointer :: snow_pcw(:)    ! snow from Hack   convection
   real(r8), pointer :: vmag_gust(:)

   real(r8) :: rad_fac    = 1.0_r8     ! scale factor for radiative flux nudging 

   !Local variables 
   integer  :: i, k, m, lchnk, ncol
   integer  :: nstep                ! current timestep number
   real(r8) :: frac1, frac2 
   real(r8) :: dflwds_mod(pcols)    ! Downward longwave flux at surface
   real(r8) :: dnetsw_mod(pcols)    ! Surface solar absorbed flux (Total sky; downward - upward)
   real(r8) :: sols_mod(pcols)      ! Direct downward shortwave flux, UV/vis
   real(r8) :: soll_mod(pcols)      ! Direct downward shortwave flux, near-IR
   real(r8) :: solsd_mod(pcols)     ! Diffuse downward shortwave flux, UV/vis
   real(r8) :: solld_mod(pcols)     ! Diffuse downward shortwave flux, near-IR
   real(r8) :: sols_obs(pcols)      ! Direct downward shortwave flux, UV/vis
   real(r8) :: soll_obs(pcols)      ! Direct downward shortwave flux, near-IR
   real(r8) :: solsd_obs(pcols)     ! Diffuse downward shortwave flux, UV/vis
   real(r8) :: solld_obs(pcols)     ! Diffuse downward shortwave flux, near-IR

   real(r8) :: precc_mod(pcols)     ! convective precip rate
   real(r8) :: precl_mod(pcols)     ! stratiform precip rate
   real(r8) :: snowc_mod(pcols)     ! convective snow rate
   real(r8) :: snowl_mod(pcols)     ! stratiform snow rate

   real(r8) :: radmask(pcols)
   real(r8) :: factor(pcols)       
   real(r8) :: psmod(pcols) 
   real(r8) :: zsmod(pcols)  
   real(r8) :: psobs(pcols)   
   real(r8) :: zsobs(pcols)       

   real(r8) :: zbmod(pcols)
   real(r8) :: pbmod(pcols)
   real(r8) :: ubmod(pcols)
   real(r8) :: vbmod(pcols)
   real(r8) :: tbmod(pcols)
   real(r8) :: thbmod(pcols)
   real(r8) :: qbmod(pcols)
   real(r8) :: vmag(pcols)
   real(r8) :: pslmod(pcols)

   real(r8) :: umod(pcols,pver)
   real(r8) :: vmod(pcols,pver)
   real(r8) :: tmod(pcols,pver)
   real(r8) :: qmod(pcols,pver)
   real(r8) :: zi_mod(pcols,pverp)
   real(r8) :: zm_mod(pcols,pver)
   real(r8) :: rairv(pcols,pver)
   real(r8) :: zvirv(pcols,pver)

   real(r8) :: precc_tend(pcols)
   real(r8) :: precl_tend(pcols)
   real(r8) :: snowc_tend(pcols)
   real(r8) :: snowl_tend(pcols)
   real(r8) :: soll_tend(pcols)
   real(r8) :: sols_tend(pcols)
   real(r8) :: solld_tend(pcols)
   real(r8) :: solsd_tend(pcols)
   real(r8) :: netsw_tend(pcols)
   real(r8) :: flwds_tend(pcols)

   lchnk = state%lchnk
   ncol  = state%ncol
   nstep = get_nstep()

   precc_tend = 0._r8 
   precl_tend = 0._r8
   snowc_tend = 0._r8
   snowl_tend = 0._r8
   netsw_tend = 0._r8
   flwds_tend = 0._r8
   soll_tend  = 0._r8
   sols_tend  = 0._r8
   solld_tend = 0._r8
   solsd_tend = 0._r8
   soll_obs   = 0._r8
   sols_obs   = 0._r8
   solld_obs  = 0._r8
   solsd_obs  = 0._r8
   factor     = 1.0_r8

   psmod(:ncol) = state%ps(:ncol)
   zsmod(:ncol) = state%phis(:ncol)
   psobs(:ncol) = Target_PS(:ncol,lchnk)
   zsobs(:ncol) = Target_PHIS(:ncol,lchnk)

   ! scaling factor based on surface pressure 
   do i = 1, ncol
     if ((Nudge_SRF_PSWgt_On) .and. (abs(psmod(i)-psobs(i)) > 2.E2_r8)) then 
         factor(i) = factor(i) * 2.E2_r8 / abs(psmod(i)-psobs(i))
     end if 
   end do 

   !surface state variables 
   if (Nudge_SRF_On .and. Nudge_SRF_State_On) then
     !updated the state due to nudging 
     do k=1,pver
       do i=1,ncol
         umod(i,k) = state%u(i,k)   + Nudge_Ustep(i,k,lchnk)*factor(i)*dtime 
         vmod(i,k) = state%v(i,k)   + Nudge_Vstep(i,k,lchnk)*factor(i)*dtime 
         tmod(i,k) = state%t(i,k)   + Nudge_Tstep(i,k,lchnk)*factor(i)*dtime 
         qmod(i,k) = state%q(i,k,1) + Nudge_Qstep(i,k,lchnk)*factor(i)*dtime 
       end do
     end do 

     !calculate geopotential height from updated temperature 
     do i = 1, ncol
       do k = 1, pver
         rairv(i,k) = rair
         zvirv(i,k) = zvir
       end do
     end do
     call geopotential_t(state%lnpint, state%lnpmid, state%pint, state%pmid,&
                         state%pdel  , state%rpdel ,  tmod     , qmod      ,&
                         rairv, gravit, zvirv, zi_mod, zm_mod, ncol)

     !update sea level pressure 
     call cpslec (ncol, state%pmid, state%phis, state%ps, tmod, pslmod, gravit, rair)

     !extract the bottom model level values required by land model 
     do i=1,ncol
       ubmod(i)  = umod(i,pver)
       vbmod(i)  = vmod(i,pver)
       tbmod(i)  = tmod(i,pver)
       qbmod(i)  = qmod(i,pver)
       zbmod(i)  = zm_mod(i,pver)
       pbmod(i)  = state%pmid(i,pver)
       thbmod(i) = tbmod(i)*state%exner(i,pver)
     end do 

     !adds gustiness (adapted from camsrfexch.F90 )
     vmag_gust_idx = pbuf_get_index('vmag_gust')
     call pbuf_get_field(pbuf, vmag_gust_idx, vmag_gust)
     do i=1,ncol
        vmag(i) = max(1.e-5_r8,sqrt( ubmod(i)**2._r8 + vbmod(i)**2._r8))
     end do

   end if

   !get the model precipitation 
   if (Nudge_SRF_On .and. Nudge_SRF_Prec_On) then
     prec_dp_idx  = pbuf_get_index('PREC_DP')
     snow_dp_idx  = pbuf_get_index('SNOW_DP')
     prec_sh_idx  = pbuf_get_index('PREC_SH')
     snow_sh_idx  = pbuf_get_index('SNOW_SH')
     prec_sed_idx = pbuf_get_index('PREC_SED')
     snow_sed_idx = pbuf_get_index('SNOW_SED')
     prec_pcw_idx = pbuf_get_index('PREC_PCW')
     snow_pcw_idx = pbuf_get_index('SNOW_PCW')

     call pbuf_get_field(pbuf, prec_dp_idx, prec_dp)
     call pbuf_get_field(pbuf, snow_dp_idx, snow_dp)
     call pbuf_get_field(pbuf, prec_sh_idx, prec_sh)
     call pbuf_get_field(pbuf, snow_sh_idx, snow_sh)
     call pbuf_get_field(pbuf, prec_sed_idx, prec_sed)
     call pbuf_get_field(pbuf, snow_sed_idx, snow_sed)
     call pbuf_get_field(pbuf, prec_pcw_idx, prec_pcw)
     call pbuf_get_field(pbuf, snow_pcw_idx, snow_pcw)

     !Precipitation rates (multi-process)
     precc_mod(:ncol) = prec_dp(:ncol)  + prec_sh(:ncol)
     precl_mod(:ncol) = prec_sed(:ncol) + prec_pcw(:ncol)
     snowc_mod(:ncol) = snow_dp(:ncol)  + snow_sh(:ncol)
     snowl_mod(:ncol) = snow_sed(:ncol) + snow_pcw(:ncol)

     !Derive the nudging tendency
     do i = 1, ncol
       !Note: *10-3 is to convert kg/m2/s (default unit in ERA5) to m/s  (unit in  model) 
       !Note: this conversion depends on the unit in the reanalysis data
       if (Target_PRECC(i,lchnk).ge.0._r8 ) then
         precc_tend(i) = (Target_PRECC(i,lchnk)*1.e-3_r8 - precc_mod(i))  * Nudge_SRFtau(i,lchnk)
       end if 
       if (Target_PRECL(i,lchnk).ge.0._r8 ) then
         precl_tend(i) = (Target_PRECL(i,lchnk)*1.e-3_r8 - precl_mod(i))  * Nudge_SRFtau(i,lchnk)
       end if 
       if (Target_PRECSC(i,lchnk).ge.0._r8 ) then
         snowc_tend(i) = (Target_PRECSC(i,lchnk)*1.e-3_r8 - snowc_mod(i)) * Nudge_SRFtau(i,lchnk)
       end if 
       if (Target_PRECSL(i,lchnk).ge.0._r8 ) then
         snowl_tend(i) = (Target_PRECSL(i,lchnk)*1.e-3_r8 - snowl_mod(i)) * Nudge_SRFtau(i,lchnk)
       end if 
     end do 

   end if 

   if (Nudge_SRF_On .and. Nudge_SRF_RadFlux_On) then

     !Specify a scaling factor for the nudging 
     do i = 1, ncol
       radmask(i) = rad_fac
       !Mask out the ocean region
       !if (cam_in%landfrac(i).ge.0.001) then
       !  radmask(i)=radmask(i)*0._r8
       !endif
     end do 

     sols_mod(1:ncol)   = cam_out%sols(1:ncol)   ! Direct downward shortwave flux, UV/vis
     soll_mod(1:ncol)   = cam_out%soll(1:ncol)   ! Direct downward shortwave flux, near-IR
     solsd_mod(1:ncol)  = cam_out%solsd(1:ncol)  ! Diffuse downward shortwave flux, UV/vis
     solld_mod(1:ncol)  = cam_out%solld(1:ncol)  ! Diffuse downward shortwave flux, near-IR 
     dflwds_mod(1:ncol) = cam_out%flwds(1:ncol)  ! Downward longwave flux at surface
     dnetsw_mod(1:ncol) = cam_out%netsw(1:ncol)  ! Surface solar absorbed flux (Total sky; downward - upward)

     !Derived the SOLS, SOLL, SOLSD and SOLLD. Assume that the ratio of 
     !SOLL/(SOLL+SOLS), SOLS/(SOLL+SOLS), SOLLD/(SOLLD+SOLSD), SOLSD/(SOLLD+SOLSD) 
     !are the same in model and observations 
     do i = 1,ncol
       frac1 = cam_out%sols(i)/(cam_out%sols(i)+cam_out%soll(i))
       frac2 = (cam_out%sols(i)+cam_out%solsd(i)) & 
               /(cam_out%solsd(i)+cam_out%solld(i)+ cam_out%sols(i)+cam_out%soll(i))
       sols_obs(i)  = Target_FSDSD(i,lchnk)*frac1            ! direct UV/vis 
       soll_obs(i)  = Target_FSDSD(i,lchnk)*(1.0_r8 - frac1) ! direct near-IR 
       solsd_obs(i) = Target_FSDS(i,lchnk)*frac2 &           ! diffuse UV/vis
                      - Target_FSDSD(i,lchnk)*frac1 
       solld_obs(i) = Target_FSDS(i,lchnk)*(1.0_r8 - frac2) &! diffuse near-IR
                      - Target_FSDSD(i,lchnk)*(1.0_r8 - frac1) 
       if (sols_obs(i).lt.0._r8)  sols_obs(i)  = 0._r8
       if (soll_obs(i).lt.0._r8)  soll_obs(i)  = 0._r8
       if (solsd_obs(i).lt.0._r8) solsd_obs(i) = 0._r8
       if (solld_obs(i).lt.0._r8) solld_obs(i) = 0._r8
     end do 

     !Derive the nudging tendency
     do i = 1, ncol

       if((sols_mod(i).ge.0._r8).and.(sols_obs(i).ge.0._r8)) then 
         sols_tend(i)  = (sols_obs(i)  - sols_mod(i)) * radmask(i) * Nudge_SRFtau(i,lchnk)
       end if 

       if((soll_mod(i).ge.0._r8).and.(soll_obs(i).ge.0._r8)) then
         soll_tend(i)  = (soll_obs(i)  - soll_mod(i)) * radmask(i) * Nudge_SRFtau(i,lchnk)
       end if 

       if((solsd_mod(i).ge.0._r8).and.(solsd_obs(i).ge.0._r8)) then
         solsd_tend(i) = (solsd_obs(i) - solsd_mod(i)) * radmask(i) * Nudge_SRFtau(i,lchnk)
       end if 

       if((solld_mod(i).ge.0._r8).and.(solld_obs(i).ge.0._r8)) then
         solld_tend(i) = (solld_obs(i) - solld_mod(i)) * radmask(i) * Nudge_SRFtau(i,lchnk)
       end if 

       if((dnetsw_mod(i).ge.0._r8).and.(Target_NETSW(i,lchnk).ge.0._r8)) then
         netsw_tend(i) = (Target_NETSW(i,lchnk) - dnetsw_mod(i)) * radmask(i) * Nudge_SRFtau(i,lchnk)
       end if 

       flwds_tend(i) = (Target_FLWDS(i,lchnk) - dflwds_mod(i)) * radmask(i) * Nudge_SRFtau(i,lchnk)

     end do

   end if

   !Start to update the variables for land surface model
   !surface state variables 
   if ((nstep > 0) .and. (Nudge_SRF_On) .and. (Nudge_SRF_State_On)) then
     do i=1,ncol
        cam_out%ubot(i)  = ubmod(i)*((vmag_gust(i)+vmag(i))/vmag(i))
        cam_out%vbot(i)  = vbmod(i)*((vmag_gust(i)+vmag(i))/vmag(i))
        cam_out%tbot(i)  = tbmod(i)  
        cam_out%thbot(i) = thbmod(i)
        cam_out%zbot(i)  = zbmod(i)
        cam_out%pbot(i)  = pbmod(i)
        cam_out%rho(i)   = cam_out%pbot(i)/(rair*cam_out%tbot(i))
        cam_out%qbot(i,1)= qbmod(i)
     end do
   end if 

   !Set the UV/vis and near-IR direct and dirruse downward shortwave flux at surface
   if ((nstep > 0) .and. (Nudge_SRF_On) .and. (Nudge_SRF_RadFlux_On)) then
     !Update the quantity in cam_out
     do i = 1,ncol
       cam_out%sols(i)  = sols_mod(i)    + sols_tend(i)  * factor(i) * dtime 
       cam_out%soll(i)  = soll_mod(i)    + soll_tend(i)  * factor(i) * dtime 
       cam_out%solsd(i) = solsd_mod(i)   + solsd_tend(i) * factor(i) * dtime 
       cam_out%solld(i) = solld_mod(i)   + solld_tend(i) * factor(i) * dtime 
       cam_out%netsw(i) = dnetsw_mod(i)  + netsw_tend(i) * factor(i) * dtime 
       cam_out%flwds(i) = dflwds_mod(i)  + flwds_tend(i) * factor(i) * dtime 
     end do
   end if

   ! Precipation and snow rates from shallow convection, deep convection and stratiform processes.
   ! Compute total convective and stratiform precipitation and snow rates
   if ((nstep > 0) .and. Nudge_SRF_On .and. Nudge_SRF_Prec_On) then 
     do i=1,ncol
       cam_out%precc (i) = precc_mod(i) + precc_tend(i) * factor(i) * dtime
       cam_out%precl (i) = precl_mod(i) + precl_tend(i) * factor(i) * dtime
       cam_out%precsc(i) = snowc_mod(i) + snowc_tend(i) * factor(i) * dtime
       cam_out%precsl(i) = snowl_mod(i) + snowl_tend(i) * factor(i) * dtime
       ! Safeguard check to ensure precipitation rate > 0
       if (cam_out%precc(i) .lt.0._r8) cam_out%precc(i)=0._r8
       if (cam_out%precl(i) .lt.0._r8) cam_out%precl(i)=0._r8
       if (cam_out%precsc(i).lt.0._r8) cam_out%precsc(i)=0._r8
       if (cam_out%precsl(i).lt.0._r8) cam_out%precsl(i)=0._r8
       if (cam_out%precsc(i).gt.cam_out%precc(i)) cam_out%precsc(i)=cam_out%precc(i)
       if (cam_out%precsl(i).gt.cam_out%precl(i)) cam_out%precsl(i)=cam_out%precl(i)
     end do

     ! total snowfall rate: needed by slab ocean model
     prcsnw(:ncol,lchnk) = cam_out%precsc(:ncol) + cam_out%precsl(:ncol)
   end if

   !Output the diagnostics 
   call outfld('PRECC_bf_ndg' ,precc_mod(:),pcols,lchnk)
   call outfld('PRECL_bf_ndg' ,precl_mod(:),pcols,lchnk)
   call outfld('PRECSC_bf_ndg',snowc_mod(:),pcols,lchnk)
   call outfld('PRECSL_bf_ndg',snowl_mod(:),pcols,lchnk)
   call outfld('SOLL_bf_ndg'  ,soll_mod(:),pcols,lchnk)
   call outfld('SOLS_bf_ndg'  ,sols_mod(:),pcols,lchnk)
   call outfld('SOLLD_bf_ndg' ,solld_mod(:),pcols,lchnk)
   call outfld('SOLSD_bf_ndg' ,solsd_mod(:),pcols,lchnk)
   call outfld('NETSW_bf_ndg' ,dnetsw_mod(:),pcols,lchnk)
   call outfld('FLWDS_bf_ndg' ,dflwds_mod(:),pcols,lchnk)

   call outfld('PRECC_af_ndg' ,cam_out%precc,pcols,lchnk)
   call outfld('PRECL_af_ndg' ,cam_out%precl,pcols,lchnk)
   call outfld('PRECSC_af_ndg',cam_out%precsc,pcols,lchnk)
   call outfld('PRECSL_af_ndg',cam_out%precsl,pcols,lchnk)
   call outfld('SOLL_af_ndg'  ,cam_out%soll,pcols,lchnk)
   call outfld('SOLS_af_ndg'  ,cam_out%sols,pcols,lchnk)
   call outfld('SOLLD_af_ndg' ,cam_out%solld,pcols,lchnk)
   call outfld('SOLSD_af_ndg' ,cam_out%solsd,pcols,lchnk)
   call outfld('NETSW_af_ndg' ,cam_out%netsw,pcols,lchnk)
   call outfld('FLWDS_af_ndg' ,cam_out%flwds,pcols,lchnk)

   call outfld('PRECC_ref' ,Target_PRECC(:,lchnk),pcols,lchnk)
   call outfld('PRECL_ref' ,Target_PRECL(:,lchnk),pcols,lchnk)
   call outfld('PRECSC_ref',Target_PRECSC(:,lchnk),pcols,lchnk)
   call outfld('PRECSL_ref',Target_PRECSL(:,lchnk),pcols,lchnk)
   call outfld('NETSW_ref' ,Target_NETSW(:,lchnk),pcols,lchnk)
   call outfld('FLWDS_ref' ,Target_FLWDS(:,lchnk),pcols,lchnk)
   call outfld('SOLL_ref'  ,soll_obs(:),pcols,lchnk)
   call outfld('SOLS_ref'  ,sols_obs(:),pcols,lchnk)
   call outfld('SOLLD_ref' ,solld_obs(:),pcols,lchnk)
   call outfld('SOLSD_ref' ,solsd_obs(:),pcols,lchnk)

   call outfld('Nudge_PRECC' ,precc_tend,pcols,lchnk)
   call outfld('Nudge_PRECL' ,precl_tend,pcols,lchnk)
   call outfld('Nudge_PRECSC',snowc_tend,pcols,lchnk)
   call outfld('Nudge_PRECSL',snowl_tend,pcols,lchnk)
   call outfld('Nudge_SOLL'  ,soll_tend,pcols,lchnk)
   call outfld('Nudge_SOLS'  ,sols_tend,pcols,lchnk)
   call outfld('Nudge_SOLLD' ,solld_tend,pcols,lchnk)
   call outfld('Nudge_SOLSD' ,solsd_tend,pcols,lchnk)
   call outfld('Nudge_NETSW' ,netsw_tend,pcols,lchnk)
   call outfld('Nudge_FLWDS' ,flwds_tend,pcols,lchnk)

   return
  end subroutine  !nudging_update_land_surface 

  !================================================================
  subroutine nudging_update_srf_flux(state, cam_in, cam_out, dtime)
   use physics_types,  only: physics_state
   use camsrfexch   ,  only: cam_out_t, cam_in_t
   use ppgrid       ,  only: pver,pcols,pverp
   use cam_history  ,  only: outfld
   use physconst    ,  only: cpair, gravit, rair, rga, stebol, zvir
   use physics_buffer, only: pbuf_get_index, pbuf_old_tim_idx, pbuf_get_field, &
                             pbuf_set_field, physics_buffer_desc
   use comsrf,         only: prcsnw


   ! Arguments
   !-----------
   type(physics_buffer_desc), pointer :: pbuf(:)
   type(physics_state), intent(in)    :: state
   type(cam_out_t),     intent(in)    :: cam_out
   type(cam_in_t),      intent(inout) :: cam_in
   real(r8),            intent(in)    :: dtime

   !Local variables 
   integer :: i, k, m, lchnk, ncol

   real(r8) :: lhf_mod(pcols)      ! latent heat fluxes 
   real(r8) :: shf_mod(pcols)      ! sensible heat fluxes
   real(r8) :: qflx_mod(pcols)     ! water vapor fluxes 

   real(r8) :: factor(pcols)
   real(r8) :: psmod(pcols)
   real(r8) :: zsmod(pcols)
   real(r8) :: psobs(pcols)
   real(r8) :: zsobs(pcols)

   real(r8) :: lhf_tend(pcols)
   real(r8) :: shf_tend(pcols)
   real(r8) :: qflx_tend(pcols)

   lchnk = state%lchnk
   ncol  = state%ncol

   lhf_tend(:ncol)  = 0._r8
   shf_tend(:ncol)  = 0._r8
   qflx_tend(:ncol) = 0._r8

   psmod(:ncol)  = state%ps(:ncol)
   zsmod(:ncol)  = state%phis(:ncol)
   psobs(:ncol)  = Target_PS(:ncol,lchnk)
   zsobs(:ncol)  = Target_PHIS(:ncol,lchnk)

   qflx_mod(:ncol) = cam_in%cflx(:ncol,1) 
   lhf_mod(:ncol)  = cam_in%lhf(:ncol) 
   shf_mod(:ncol)  = cam_in%shf(:ncol)  

   factor(:ncol) = 1.0_r8
   if (Nudge_SRF_PSWgt_On) then
     do i = 1, ncol
       if ( abs(psmod(i)-psobs(i)) > 2.E2_r8 ) then
         factor(i) = 2.E2_r8 / abs(psmod(i)-psobs(i))
       end if
     enddo
   end if

   !Derive the nudging tendency 
   do i = 1,ncol
     lhf_tend(i)  = (-Target_LHFLX(i,lchnk) - lhf_mod(i)) * Nudge_SRFtau(i,lchnk)
     shf_tend(i)  = (-Target_SHFLX(i,lchnk) - shf_mod(i)) * Nudge_SRFtau(i,lchnk)    
     qflx_tend(i) = (-Target_EVAP(i,lchnk) - qflx_mod(i)) * Nudge_SRFtau(i,lchnk)
   end do 
   
   !Update the quantity in cam_out
   do i = 1,ncol
     cam_in%cflx(i,1) = qflx_mod(i) + qflx_tend(i) * factor(i) * dtime
     cam_in%shf(i)    = shf_mod(i)  + shf_tend(i)  * factor(i) * dtime 
     cam_in%lhf(i)    = lhf_mod(i)  + lhf_tend(i)  * factor(i) * dtime
   end do

   call outfld('LHFLX_bf_ndg' ,lhf_mod,pcols,lchnk)
   call outfld('SHFLX_bf_ndg' ,shf_mod,pcols,lchnk)
   call outfld('QFLX_bf_ndg'  ,qflx_mod,pcols,lchnk)

   call outfld('LHFLX_af_ndg' ,cam_in%lhf,pcols,lchnk)
   call outfld('SHFLX_af_ndg' ,cam_in%shf,pcols,lchnk)
   call outfld('QFLX_af_ndg'  ,cam_in%cflx(:,1),pcols,lchnk)

   call outfld('LHFLX_ref' ,-Target_SHFLX(:,lchnk),pcols,lchnk)
   call outfld('SHFLX_ref' ,-Target_LHFLX(:,lchnk),pcols,lchnk)
   call outfld('QFLX_ref'  ,-Target_EVAP(:,lchnk),pcols,lchnk)

   call outfld('Nudge_LHFLX' ,lhf_mod(:),pcols,lchnk)
   call outfld('Nudge_SHFLX' ,shf_mod(:),pcols,lchnk)
   call outfld('Nudge_QFLX'  ,qflx_mod(:),pcols,lchnk)

   return
  end subroutine  !nudging_update_srf_flux

  !================================================================
  subroutine nudging_timestep_tend(phys_state,phys_tend, dtime)
   !
   ! NUDGING_TIMESTEP_TEND:
   !                If Nudging is ON, return the Nudging contributions
   !                to forcing using the current contents of the Nudge
   !                arrays. Send output to the cam history module as well.
   !===============================================================
   use physics_types,only: physics_state,physics_ptend,physics_ptend_init
   use constituents ,only: cnst_get_ind,pcnst
   use ppgrid       ,only: pver,pverp,pcols,begchunk,endchunk
   use cam_history  ,only: outfld
   use physconst    ,only: rga, cpair, gravit, rair, zvir, cappa
   use hycoef       ,only: hycoef_init, hyam, hybm, hyai, hybi, ps0

   ! Arguments
   !-------------
   type(physics_state), intent(inout) :: phys_state
   type(physics_ptend), intent(out):: phys_tend
  
   real(r8), intent(in):: dtime 

   ! Local values
   !--------------------
   real(r8):: nudge_q
   integer ixcldliq, ixcldice
   integer indw,ncol,lchnk
   logical lq(pcnst)
   integer Year, Month, Day, Sec
   integer i,k,m

   call cnst_get_ind('Q',indw)
   call cnst_get_ind('CLDLIQ', ixcldliq)
   call cnst_get_ind('CLDICE', ixcldice)

   lq(:)   =.false.
   lq(indw)=.true.
   call physics_ptend_init(phys_tend,phys_state%psetcols,'nudging',lu=.true.,lv=.true.,ls=.true.,lq=lq)

   if(Nudge_ON) then
     lchnk=phys_state%lchnk
     ncol =phys_state%ncol

     call get_curr_date(Year,Month,Day,Sec)
     ! the code below is necessary for a correct IMT nudging with a continue/restart run
     if ( Nudge_Method .eq. 'IMT' .and. (mod(Sec,Nudge_Step) .ne. 0) ) then
         Nudge_Ustep(:ncol,:pver,lchnk) =0._r8
         Nudge_Vstep(:ncol,:pver,lchnk) =0._r8
         Nudge_Tstep(:ncol,:pver,lchnk) =0._r8
         Nudge_Qstep(:ncol,:pver,lchnk) =0._r8
         Nudge_PSstep(:ncol,lchnk)      =0._r8
     end if

     if (Nudge_Uprof .ne. 0) then
        phys_tend%u(:ncol,:pver) = Nudge_Ustep(:ncol,:pver,lchnk)
     end if

     if (Nudge_Vprof .ne. 0) then
        phys_tend%v(:ncol,:pver) = Nudge_Vstep(:ncol,:pver,lchnk)
     end if

     if (Nudge_Tprof .ne. 0) then
        phys_tend%s(:ncol,:pver) = Nudge_Tstep(:ncol,:pver,lchnk)*cpair
     end if

     if (Nudge_Qprof .ne. 0) then
        phys_tend%q(:ncol,:pver,indw) = Nudge_Qstep(:ncol,:pver,lchnk)
     end if

     if (Nudge_PS_On) then
       !update ps 
       do i = 1,ncol
          phys_state%ps(i) = phys_state%ps(i) + Nudge_PSstep(i,lchnk)*dtime
       end do
       !Pdel etc.
       do k = 1, pver
         do i = 1, ncol
           phys_state%pdel(i,k)  = phys_state%pdel(i,k)  &
                                  + (hybi(k+1) -hybi(k))*Nudge_PSstep(i,lchnk)*dtime
           phys_state%rpdel(i,k) = 1._r8/phys_state%pdel(i,k)
         end do
       end do
     end if

     !adjust constitutes 
     if (Nudge_Q_Adjust_On) then
       do i = 1, ncol
         do k = 1, pver
           if (abs(Nudge_Qstep(i,k,lchnk)) > 0._r8) then
             nudge_q = phys_state%q(i,k,1)
             do m = 2, pcnst
               phys_state%q(i,k,m) = phys_state%q(i,k,m)*phys_state%pdel(i,k)
             end do
             phys_state%pdel(i,k)  = phys_state%pdel(i,k)*(1.0_r8 - nudge_q)
             nudge_q               = nudge_q + Nudge_Qstep(i,k,lchnk) * dtime
             phys_state%pdel(i,k)  = phys_state%pdel(i,k)/(1.0_r8 - nudge_q)
             phys_state%rpdel(i,k) = 1.0_r8/phys_state%pdel(i,k)
             do m = 2, pcnst
               phys_state%q(i,k,m) = phys_state%q(i,k,m)/phys_state%pdel(i,k)
             end do
           end if
         end do
       end do
     end if

   end if 

   ! End Routine
   !------------
   return
  end subroutine ! nudging_timestep_tend
  !================================================================


  !================================================================
  subroutine nudging_update_analyses_se(anal_file)
   !
   ! NUDGING_UPDATE_ANALYSES_SE:
   !                 Open the given analyses data file, read in
   !                 U,V,T,Q, and PS values and then distribute
   !                 the values to all of the chunks.
   !===============================================================
!   use wrap_nf
   use ppgrid ,only: pver
   use netcdf
   use filenames ,only: interpret_filename_spec

   ! Arguments
   !-------------
   character(len=*),intent(in):: anal_file

   ! Local values
   !-------------
   integer lev
   integer ncol,plev,istat
   integer ncid, varid
!   real(r8) PSanal(Nudge_ncol)
   real(r8) Lat_anal(Nudge_ncol)
   real(r8) Lon_anal(Nudge_ncol)
!   real(r8) :: Xanal(Nudge_ncol,Nudge_nlev)

   integer :: cnt3(3)               ! array of counts for each dimension
   integer :: strt3(3)              ! array of starting indices
   integer :: cnt2(2)               ! array of counts for each dimension
   integer :: strt2(2)              ! array of starting indices
   integer :: n, n_cnt, ncid1, ind
   integer :: timesiz               ! size of time dimension on dataset
   integer :: Year, Month, Day, Sec

   ! Check the existence of the analyses file; broadcast the file status to
   ! all the other MPI nodes. If the file is not there, then just return.
   !------------------------------------------------------------------------
   if(masterproc) then
     inquire(FILE=trim(anal_file),EXIST=Nudge_File_Present)
   endif
#ifdef SPMD
   call mpibcast(Nudge_File_Present, 1, mpilog, 0, mpicom)
#endif
   if (trim(Nudge_Method).eq. "MLTBC") then
     if(masterproc) then
        write(iulog,*) 'Warning: using Machine Learning model to predict nudging tendency'
        write(iulog,*) 'Warning: No need to read reference data, return'
     end if
     return
   else
     if(.not.Nudge_File_Present) return
   end if

   ! masterporc does all of the work here
   !-----------------------------------------
   if(masterproc) then

     ! Open the given file
     !-----------------------
     istat=nf90_open(trim(anal_file),NF90_NOWRITE,ncid)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*)'NF90_OPEN: failed for file ',trim(anal_file)
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_SE_OPEN')
     endif

     ! Read in Dimensions
     !--------------------
!     call wrap_inq_dimid (ncid,'ncol',varid)
!     call wrap_inq_dimlen(ncid,varid,ncol)
     istat=nf90_inq_dimid(ncid,'ncol',varid)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_SE')
     endif
     istat=nf90_inquire_dimension(ncid,varid,len=ncol)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_SE')
     endif

!     call wrap_inq_dimid (ncid,'lev',varid)
!     call wrap_inq_dimlen(ncid,varid,plev)
     istat=nf90_inq_dimid(ncid,'lev',varid)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_SE')
     endif
     istat=nf90_inquire_dimension(ncid,varid,len=plev)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_SE')
     endif

!     call wrap_inq_varid(ncid,'lon',varid)
!     call wrap_get_var_realx(ncid,varid,Lon_anal)
     istat=nf90_inq_varid(ncid,'lon',varid)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_SE')
     endif
     istat=nf90_get_var(ncid,varid,Lon_anal)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_SE')
     endif

!     call wrap_inq_varid(ncid,'lat',varid)
!     call wrap_get_var_realx(ncid,varid,Lat_anal)
     istat=nf90_inq_varid(ncid,'lat',varid)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_SE')
     endif
     istat=nf90_get_var(ncid,varid,Lat_anal)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_SE')
     endif

     if ((Nudge_ncol.ne.ncol) .or. (plev.ne.pver)) then
        write(iulog,*) 'ERROR: nudging_update_analyses_se: ncol=',ncol,' Nudge_ncol=',Nudge_ncol
        write(iulog,*) 'ERROR: nudging_update_analyses_se: plev=',plev,' pver=',pver
        call endrun('nudging_update_analyses_se: analyses dimension mismatch')
     end if

   end if ! (masterproc) then

   select case (Nudge_Method)
      case ('Step')
           if (masterproc) then
              call get_curr_date(Year,Month,Day,Sec)
              n_cnt = Sec/Nudge_Step + 1
              n_cnt = n_cnt + 1                      ! nudge to future model time step
              if (n_cnt .gt. Nudge_File_Ntime) then
                  n_cnt = 1
              end if
              strt3(1) = 1
              strt3(2) = 1
              strt3(3) = n_cnt
              cnt3(1)  = ncol
              cnt3(2)  = pver
              cnt3(3)  = 1
              strt2(1) = 1
              strt2(2) = n_cnt
              cnt2(1)  = ncol
              cnt2(2)  = 1
           end if
           call read_and_scatter_se(ncid, 'U', strt3, cnt3, Target_U)
           call read_and_scatter_se(ncid, 'V', strt3, cnt3, Target_V)
           call read_and_scatter_se(ncid, 'T', strt3, cnt3, Target_T)
           call read_and_scatter_se(ncid, 'Q', strt3, cnt3, Target_Q)
           call read_and_scatter_se_2d(ncid, 'PS', strt2, cnt2, Target_PS)
           call read_and_scatter_se_2d(ncid, 'PHIS', strt2, cnt2, Target_PHIS)

      case ('IMT')
           call get_curr_date(Year,Month,Day,Sec)
           if (mod(Sec,Nudge_Step) .ne. 0) return  ! ensure that intermittent simulations work for restart run
           if (masterproc) then
              call t_startf ('read_nudging_data')
              call open_netcdf (ncid1, -Nudge_Step, Nudge_File_Template)
              call t_stopf ('read_nudging_data')
              n_cnt = Sec/Nudge_Step + 1
              if (n_cnt .gt. Nudge_File_Ntime) then       ! account for one time slice per file
                 n_cnt = 1
              end if
              strt3(1) = 1
              strt3(2) = 1
              strt3(3) = n_cnt
              cnt3(1)  = ncol
              cnt3(2)  = pver
              cnt3(3)  = 1
              strt2(1) = 1
              strt2(2) = n_cnt
              cnt2(1)  = ncol
              cnt2(2)  = 1
           end if

           ! Use the CURR time slice for nudging
           !------------------------------------
           call read_and_scatter_se(ncid1, 'U', strt3, cnt3, Target_U)
           call read_and_scatter_se(ncid1, 'V', strt3, cnt3, Target_V)
           call read_and_scatter_se(ncid1, 'T', strt3, cnt3, Target_T)
           call read_and_scatter_se(ncid1, 'Q', strt3, cnt3, Target_Q)
           call read_and_scatter_se_2d(ncid1, 'PS', strt2, cnt2, Target_PS)
           call read_and_scatter_se_2d(ncid1, 'PHIS', strt2, cnt2, Target_PHIS)

      case ('Linear')
           ! Single time slice per file
           ! Need to open a new netcdf file to get the CURR time slice
           if (Nudge_File_Ntime .eq. 1) then
              if (first_file) then
                  first_file = .false.
                  if (masterproc) then
                     call t_startf ('read_nudging_data')
                     call open_netcdf (ncid1, -Nudge_Step, Nudge_File_Template)
                     call t_stopf ('read_nudging_data')
                     strt3(1) = 1
                     strt3(2) = 1
                     strt3(3) = 1 
                     cnt3(1)  = ncol
                     cnt3(2)  = pver
                     cnt3(3)  = 1
                     strt2(1) = 1
                     strt2(2) = 1
                     cnt2(1)  = ncol
                     cnt2(2)  = 1
                  end if

                  !-----------------------------------------                  
                  ! The start point uses the CURR time slice
                  !-----------------------------------------
                  call read_and_scatter_se(ncid1, 'U', strt3, cnt3, INTP_U(:,:,:,1))
                  call read_and_scatter_se(ncid1, 'V', strt3, cnt3, INTP_V(:,:,:,1))
                  call read_and_scatter_se(ncid1, 'T', strt3, cnt3, INTP_T(:,:,:,1))
                  call read_and_scatter_se(ncid1, 'Q', strt3, cnt3, INTP_Q(:,:,:,1))
                  call read_and_scatter_se_2d(ncid1, 'PS',   strt2, cnt2, INTP_PS(:,:,1))
                  call read_and_scatter_se_2d(ncid1, 'PHIS', strt2, cnt2, INTP_PHIS(:,:,1))

                  !-----------------------------------------               
                  ! The end point uses the NEXT time slice
                  !---------------------------------------
                  call read_and_scatter_se(ncid, 'U', strt3, cnt3, INTP_U(:,:,:,2))
                  call read_and_scatter_se(ncid, 'V', strt3, cnt3, INTP_V(:,:,:,2))
                  call read_and_scatter_se(ncid, 'T', strt3, cnt3, INTP_T(:,:,:,2))
                  call read_and_scatter_se(ncid, 'Q', strt3, cnt3, INTP_Q(:,:,:,2))
                  call read_and_scatter_se_2d(ncid, 'PS',   strt2, cnt2, INTP_PS(:,:,2))
                  call read_and_scatter_se_2d(ncid, 'PHIS', strt2, cnt2, INTP_PHIS(:,:,2))

                  call t_startf ('read_nudging_data')
                  if (masterproc) then
                      istat=nf90_close(ncid1)
                      if (istat.ne.NF90_NOERR) then
                         write(iulog,*) nf90_strerror(istat)
                         call endrun ('UPDATE_ANALYSES_SE_CLOSE_LINEAR_NETCDF')
                      end if
                  end if
                  call t_stopf ('read_nudging_data')

              else
                  ! The previous end point becomes the start point
                  ! Only need to read in the new end point
                  !-----------------------------------------------
                  INTP_U(:,:,:,1)  = INTP_U(:,:,:,2)
                  INTP_V(:,:,:,1)  = INTP_V(:,:,:,2)
                  INTP_Q(:,:,:,1)  = INTP_Q(:,:,:,2)
                  INTP_T(:,:,:,1)  = INTP_T(:,:,:,2)
                  INTP_PS(:,:,1)   = INTP_PS(:,:,2)
                  INTP_PHIS(:,:,1)   = INTP_PHIS(:,:,2)
                  if (masterproc) then
                     strt3(1) = 1
                     strt3(2) = 1
                     strt3(3) = 1
                     cnt3(1)  = ncol
                     cnt3(2)  = pver
                     cnt3(3)  = 1
                     strt2(1) = 1
                     strt2(2) = 1
                     cnt2(1)  = ncol
                     cnt2(2)  = 1
                  end if
                  call read_and_scatter_se(ncid, 'U', strt3, cnt3, INTP_U(:,:,:,2))
                  call read_and_scatter_se(ncid, 'V', strt3, cnt3, INTP_V(:,:,:,2))
                  call read_and_scatter_se(ncid, 'T', strt3, cnt3, INTP_T(:,:,:,2))
                  call read_and_scatter_se(ncid, 'Q', strt3, cnt3, INTP_Q(:,:,:,2))
                  call read_and_scatter_se_2d(ncid, 'PS',   strt2, cnt2, INTP_PS(:,:,2))
                  call read_and_scatter_se_2d(ncid, 'PHIS', strt2, cnt2, INTP_PHIS(:,:,2))

              end if    ! first file for single time slice per file

           else
              ! Multiple time slices per file
              ! The start point uses the CURR time slice
              ! The end point uses the NEXT time slice
              !-----------------------------------------
              if (first_file) then
                  first_file = .false.
                  call get_curr_date(Year,Month,Day,Sec)
                  n_cnt = Sec/Nudge_Step + 1
                  if (n_cnt .eq. Nudge_File_Ntime) then
                     if (masterproc) then
                        call t_startf ('read_nudging_data')
                        call open_netcdf (ncid1, -Nudge_Step, Nudge_File_Template)
                        call t_stopf ('read_nudging_data')
                        strt3(1) = 1
                        strt3(2) = 1
                        strt3(3) = n_cnt
                        cnt3(1)  = ncol
                        cnt3(2)  = pver
                        cnt3(3)  = 1
                        strt2(1) = 1
                        strt2(2) = n_cnt
                        cnt2(1)  = ncol
                        cnt2(2)  = 1
                     end if

                     !-----------------------------------------
                     ! The start point uses the CURR time slice
                     !-----------------------------------------
                     call read_and_scatter_se(ncid1, 'U', strt3, cnt3, INTP_U(:,:,:,1))
                     call read_and_scatter_se(ncid1, 'V', strt3, cnt3, INTP_V(:,:,:,1))
                     call read_and_scatter_se(ncid1, 'T', strt3, cnt3, INTP_T(:,:,:,1))
                     call read_and_scatter_se(ncid1, 'Q', strt3, cnt3, INTP_Q(:,:,:,1))
                     call read_and_scatter_se_2d(ncid1, 'PS',   strt2, cnt2, INTP_PS(:,:,1))
                     call read_and_scatter_se_2d(ncid1, 'PHIS', strt2, cnt2, INTP_PHIS(:,:,1))

                     if (masterproc) then
                        strt3(3) = 1
                        strt2(2) = 1
                     end if

                     !-----------------------------------------
                     ! The end point uses the NEXT time slice
                     !---------------------------------------
                     call read_and_scatter_se(ncid, 'U', strt3, cnt3, INTP_U(:,:,:,2))
                     call read_and_scatter_se(ncid, 'V', strt3, cnt3, INTP_V(:,:,:,2))
                     call read_and_scatter_se(ncid, 'T', strt3, cnt3, INTP_T(:,:,:,2))
                     call read_and_scatter_se(ncid, 'Q', strt3, cnt3, INTP_Q(:,:,:,2))
                     call read_and_scatter_se_2d(ncid, 'PS',   strt2, cnt2, INTP_PS(:,:,2))
                     call read_and_scatter_se_2d(ncid, 'PHIS', strt2, cnt2, INTP_PHIS(:,:,2))

                     call t_startf ('read_nudging_data')
                     if (masterproc) then
                         istat=nf90_close(ncid1)
                         if (istat.ne.NF90_NOERR) then
                            write(iulog,*) nf90_strerror(istat)
                            call endrun ('UPDATE_ANALYSES_SE_CLOSE_LINEAR_NETCDF')
                         end if
                     end if
                     call t_stopf ('read_nudging_data')
                  else
                     ! two time slices are in the same nudging data file
                     do n = n_cnt, n_cnt+1
                        if (masterproc) then
                           strt3(1) = 1
                           strt3(2) = 1
                           strt3(3) = n
                           cnt3(1)  = ncol
                           cnt3(2)  = pver
                           cnt3(3)  = 1
                           strt2(1) = 1
                           strt2(2) = n
                           cnt2(1)  = ncol
                           cnt2(2)  = 1
                        end if
                        ind = n - n_cnt + 1
                        call read_and_scatter_se(ncid, 'U', strt3, cnt3, INTP_U(:,:,:,ind))
                        call read_and_scatter_se(ncid, 'V', strt3, cnt3, INTP_V(:,:,:,ind))
                        call read_and_scatter_se(ncid, 'T', strt3, cnt3, INTP_T(:,:,:,ind))
                        call read_and_scatter_se(ncid, 'Q', strt3, cnt3, INTP_Q(:,:,:,ind))
                        call read_and_scatter_se_2d(ncid, 'PS',   strt2, cnt2, INTP_PS(:,:,ind))
                        call read_and_scatter_se_2d(ncid, 'PHIS', strt2, cnt2, INTP_PHIS(:,:,ind))
                     end do
                  end if

              else

                  call get_curr_date(Year,Month,Day,Sec)
                  n_cnt = Sec/Nudge_Step + 1
                  n_cnt = n_cnt + 1                ! open the nudging data at future model time step
                  if (n_cnt .gt. Nudge_File_Ntime) then
                      n_cnt = 1
                  end if
                  ! The previous end point becomes the start point
                  ! Only need to read in the new end point
                  !-----------------------------------------------
                  INTP_U(:,:,:,1)  = INTP_U(:,:,:,2)
                  INTP_V(:,:,:,1)  = INTP_V(:,:,:,2)
                  INTP_Q(:,:,:,1)  = INTP_Q(:,:,:,2)
                  INTP_T(:,:,:,1)  = INTP_T(:,:,:,2)
                  INTP_PS(:,:,1)   = INTP_PS(:,:,2)
                  INTP_PHIS(:,:,1) = INTP_PHIS(:,:,2)
                  if (masterproc) then
                     strt3(1) = 1
                     strt3(2) = 1
                     strt3(3) = n_cnt
                     cnt3(1)  = ncol
                     cnt3(2)  = pver
                     cnt3(3)  = 1
                     strt2(1) = 1
                     strt2(2) = n_cnt
                     cnt2(1)  = ncol
                     cnt2(2)  = 1
                  end if
                  call read_and_scatter_se(ncid, 'U', strt3, cnt3, INTP_U(:,:,:,2))
                  call read_and_scatter_se(ncid, 'V', strt3, cnt3, INTP_V(:,:,:,2))
                  call read_and_scatter_se(ncid, 'T', strt3, cnt3, INTP_T(:,:,:,2))
                  call read_and_scatter_se(ncid, 'Q', strt3, cnt3, INTP_Q(:,:,:,2))
                  call read_and_scatter_se_2d(ncid, 'PS',   strt2, cnt2, INTP_PS(:,:,2))
                  call read_and_scatter_se_2d(ncid, 'PHIS', strt2, cnt2, INTP_PHIS(:,:,2))

              end if ! first_file for multiple time slices

           end if     ! single vs. multiple time slices per file

      case ('MLTBC')
          write(iulog,*) 'ERROR: when machine learning method is called, code should not attempt to read reference data'
          call endrun('nudging_update_analyses_se: error in machine learning nudging')
      case default
          write(iulog,*) 'ERROR: Unknown Input Nudge Method'
          call endrun('nudging_update_analyses_se: bad input nudge method')
   end select

   if(masterproc) then
!!    call wrap_inq_varid    (ncid,'PS',varid)
!!    call wrap_get_var_realx(ncid,varid,PSanal)
!    istat=nf90_inq_varid(ncid,'PS',varid)
!    if(istat.ne.NF90_NOERR) then
!      write(iulog,*) nf90_strerror(istat)
!      call endrun ('UPDATE_ANALYSES_SE')
!    endif
!    istat=nf90_get_var(ncid,varid,PSanal)
!    if(istat.ne.NF90_NOERR) then
!      write(iulog,*) nf90_strerror(istat)
!      call endrun ('UPDATE_ANALYSES_SE')
!    endif

     ! Close the analyses file
     !------------------------
!     call wrap_close(ncid)
     istat=nf90_close(ncid)
     if (istat.ne.NF90_NOERR) then
        write(iulog,*) nf90_strerror(istat)
        call endrun ('UPDATE_ANALYSES_SE_CLOSE')
     end if
   end if ! (masterproc) then
 ! call scatter_field_to_chunk(1,         1,1,Nudge_ncol,PSanal,Target_PS)

   ! End Routine
   !------------
   return
  end subroutine ! nudging_update_analyses_se
  !================================================================

  !================================================================
  subroutine nudging_update_srf_analyses_se(anal_srf_file)
   !
   ! NUDGING_UPDATE_SRF_ANALYSES_SE:
   !                 Open the given analyses data file, read in
   !                 surface state, fluxes and radiation information 
   !                 and then distribute the values to all of the chunks.
   !===============================================================
   use ppgrid ,only: pver
   use netcdf
   use filenames ,only: interpret_filename_spec

   ! Arguments
   !-------------
   character(len=*),intent(in):: anal_srf_file

   ! Local values
   !-------------
   integer lev
   integer ncol,plev,istat
   integer ncid, varid

   real(r8) Lat_anal(Nudge_ncol)
   real(r8) Lon_anal(Nudge_ncol)

   integer :: cnt2(2)               ! array of counts for each dimension
   integer :: strt2(2)              ! array of starting indices
   integer :: n, n_cnt, ncid1, ind
   integer :: timesiz               ! size of time dimension on dataset
   integer :: Year, Month, Day, Sec

   ! Check the existence of the analyses file; broadcast the file status to
   ! all the other MPI nodes. If the file is not there, then just return.
   !------------------------------------------------------------------------
   if(masterproc) then
     inquire(FILE=trim(anal_srf_file),EXIST=Nudge_SRF_File_Present)
   endif
#ifdef SPMD
   call mpibcast(Nudge_SRF_File_Present, 1, mpilog, 0, mpicom)
#endif
   if(.not.Nudge_SRF_File_Present) return

   ! masterporc does all of the work here
   !-----------------------------------------
   if(masterproc) then

     ! Open the given file
     !-----------------------
     istat=nf90_open(trim(anal_srf_file),NF90_NOWRITE,ncid)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*)'NF90_OPEN: failed for file ',trim(anal_srf_file)
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_LAND_SURFACE_ANALYSES_SE_OPEN')
     endif
     
     ! Read in Dimensions
     !--------------------
     istat=nf90_inq_dimid(ncid,'ncol',varid)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_LAND_SURFACE_ANALYSES_SE')
     endif
     istat=nf90_inquire_dimension(ncid,varid,len=ncol)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_LAND_SURFACE_ANALYSES_SE')
     endif

     istat=nf90_inq_varid(ncid,'lon',varid)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_LAND_SURFACE_ANALYSES_SE')
     endif
     istat=nf90_get_var(ncid,varid,Lon_anal)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_LAND_SURFACE_ANALYSES_SE')
     endif

     if (Nudge_ncol.ne.ncol) then
        write(iulog,*) 'ERROR: nudging_update_srf_analyses_se: ncol=',ncol,' Nudge_ncol=',Nudge_ncol
        call endrun('nudging_update_srf_analyses_se: analyses dimension mismatch')
     end if

   end if ! (masterproc) then

   select case (Nudge_Method)
      case ('Step')
           if (masterproc) then
              call get_curr_date(Year,Month,Day,Sec)
              n_cnt = Sec/Nudge_Step + 1
              n_cnt = n_cnt + 1                      ! nudge to future model time step
              if (n_cnt .gt. Nudge_SRF_File_Ntime) then
                  n_cnt = 1
              end if
              strt2(1) = 1
              strt2(2) = n_cnt
              cnt2(1)  = ncol
              cnt2(2)  = 1
           end if
           call read_and_scatter_se_2d(ncid, 'U10',    strt2, cnt2, Target_U10)
           call read_and_scatter_se_2d(ncid, 'V10',    strt2, cnt2, Target_V10)
           call read_and_scatter_se_2d(ncid, 'T2',     strt2, cnt2, Target_T2)
           call read_and_scatter_se_2d(ncid, 'TD2',    strt2, cnt2, Target_TD2)
           call read_and_scatter_se_2d(ncid, 'TS',     strt2, cnt2, Target_TS)
           call read_and_scatter_se_2d(ncid, 'PRECC',  strt2, cnt2, Target_PRECC)
           call read_and_scatter_se_2d(ncid, 'PRECL',  strt2, cnt2, Target_PRECL)
           call read_and_scatter_se_2d(ncid, 'PRECSC', strt2, cnt2, Target_PRECSC)
           call read_and_scatter_se_2d(ncid, 'PRECSL', strt2, cnt2, Target_PRECSL)
           call read_and_scatter_se_2d(ncid, 'ER',     strt2, cnt2, Target_EVAP)
           call read_and_scatter_se_2d(ncid, 'SHFLX',  strt2, cnt2, Target_SHFLX)
           call read_and_scatter_se_2d(ncid, 'LHFLX',  strt2, cnt2, Target_LHFLX)
           call read_and_scatter_se_2d(ncid, 'FSNS',   strt2, cnt2, Target_NETSW)
           call read_and_scatter_se_2d(ncid, 'FLDS',   strt2, cnt2, Target_FLWDS)
           call read_and_scatter_se_2d(ncid, 'FSDS',   strt2, cnt2, Target_FSDS)
           call read_and_scatter_se_2d(ncid, 'FSDSD',  strt2, cnt2, Target_FSDSD)
           call read_and_scatter_se_2d(ncid, 'FSDSUV', strt2, cnt2, Target_FSDSUV)

      case ('IMT')
           call get_curr_date(Year,Month,Day,Sec)
           if (mod(Sec,Nudge_Step) .ne. 0) return  ! ensure that intermittent simulations work for restart run
           if (masterproc) then
              call t_startf ('read_nudging_data')
              call open_netcdf (ncid1, -Nudge_Step, Nudge_SRF_File_Template)
              call t_stopf ('read_nudging_data')
              n_cnt = Sec/Nudge_Step + 1
              if (n_cnt .gt. Nudge_SRF_File_Ntime) then       ! account for one time slice per file
                 n_cnt = 1
              end if
              strt2(1) = 1
              strt2(2) = n_cnt
              cnt2(1)  = ncol
              cnt2(2)  = 1
           end if

           ! Use the CURR time slice for nudging
           !------------------------------------
           call read_and_scatter_se_2d(ncid1, 'U10',    strt2, cnt2, Target_U10)
           call read_and_scatter_se_2d(ncid1, 'V10',    strt2, cnt2, Target_V10)
           call read_and_scatter_se_2d(ncid1, 'T2',     strt2, cnt2, Target_T2)
           call read_and_scatter_se_2d(ncid1, 'TD2',    strt2, cnt2, Target_TD2)
           call read_and_scatter_se_2d(ncid1, 'TS',     strt2, cnt2, Target_TS)
           call read_and_scatter_se_2d(ncid1, 'PRECC',  strt2, cnt2, Target_PRECC)
           call read_and_scatter_se_2d(ncid1, 'PRECL',  strt2, cnt2, Target_PRECL)
           call read_and_scatter_se_2d(ncid1, 'PRECSC', strt2, cnt2, Target_PRECSC)
           call read_and_scatter_se_2d(ncid1, 'PRECSL', strt2, cnt2, Target_PRECSL)
           call read_and_scatter_se_2d(ncid1, 'ER',     strt2, cnt2, Target_EVAP)
           call read_and_scatter_se_2d(ncid1, 'SHFLX',  strt2, cnt2, Target_SHFLX)
           call read_and_scatter_se_2d(ncid1, 'LHFLX',  strt2, cnt2, Target_LHFLX)
           call read_and_scatter_se_2d(ncid1, 'FSNS',   strt2, cnt2, Target_NETSW)
           call read_and_scatter_se_2d(ncid1, 'FLDS',   strt2, cnt2, Target_FLWDS)
           call read_and_scatter_se_2d(ncid1, 'FSDS',   strt2, cnt2, Target_FSDS)
           call read_and_scatter_se_2d(ncid1, 'FSDSD',  strt2, cnt2, Target_FSDSD)
           call read_and_scatter_se_2d(ncid1, 'FSDSUV', strt2, cnt2, Target_FSDSUV)

      case ('Linear')
           ! Single time slice per file
           ! Need to open a new netcdf file to get the CURR time slice
           if (Nudge_SRF_File_Ntime .eq. 1) then
              if (first_srf_file) then
                  first_srf_file = .false.
                  if (masterproc) then
                     call t_startf ('read_nudging_data')
                     call open_netcdf (ncid1, -Nudge_Step, Nudge_SRF_File_Template)
                     call t_stopf ('read_nudging_data')
                     strt2(1) = 1
                     strt2(2) = 1
                     cnt2(1)  = ncol
                     cnt2(2)  = 1
                  end if

                  !-----------------------------------------                  
                  ! The start point uses the CURR time slice
                  !-----------------------------------------
                  call read_and_scatter_se_2d(ncid1, 'U10',    strt2, cnt2, INTP_U10(:,:,1))
                  call read_and_scatter_se_2d(ncid1, 'V10',    strt2, cnt2, INTP_V10(:,:,1))
                  call read_and_scatter_se_2d(ncid1, 'T2',     strt2, cnt2, INTP_T2(:,:,1))
                  call read_and_scatter_se_2d(ncid1, 'TD2',    strt2, cnt2, INTP_TD2(:,:,1))
                  call read_and_scatter_se_2d(ncid1, 'TS',     strt2, cnt2, INTP_TS(:,:,1))
                  call read_and_scatter_se_2d(ncid1, 'PRECC',  strt2, cnt2, INTP_PRECC(:,:,1))
                  call read_and_scatter_se_2d(ncid1, 'PRECL',  strt2, cnt2, INTP_PRECL(:,:,1))
                  call read_and_scatter_se_2d(ncid1, 'PRECSC', strt2, cnt2, INTP_PRECSC(:,:,1))
                  call read_and_scatter_se_2d(ncid1, 'PRECSL', strt2, cnt2, INTP_PRECSL(:,:,1))
                  call read_and_scatter_se_2d(ncid1, 'ER',     strt2, cnt2, INTP_EVAP(:,:,1))
                  call read_and_scatter_se_2d(ncid1, 'LHFLX',  strt2, cnt2, INTP_LHFLX(:,:,1))
                  call read_and_scatter_se_2d(ncid1, 'SHFLX',  strt2, cnt2, INTP_SHFLX(:,:,1))
                  call read_and_scatter_se_2d(ncid1, 'FSNS',   strt2, cnt2, INTP_FSNS(:,:,1))
                  call read_and_scatter_se_2d(ncid1, 'FLDS',   strt2, cnt2, INTP_FLDS(:,:,1))
                  call read_and_scatter_se_2d(ncid1, 'FSDS',   strt2, cnt2, INTP_FSDS(:,:,1))
                  call read_and_scatter_se_2d(ncid1, 'FSDSD',  strt2, cnt2, INTP_FSDSD(:,:,1))
                  call read_and_scatter_se_2d(ncid1, 'FSDSUV', strt2, cnt2, INTP_FSDSUV(:,:,1))

                  !-----------------------------------------               
                  ! The end point uses the NEXT time slice
                  !---------------------------------------
                  call read_and_scatter_se_2d(ncid, 'U10',    strt2, cnt2, INTP_U10(:,:,2))
                  call read_and_scatter_se_2d(ncid, 'V10',    strt2, cnt2, INTP_V10(:,:,2))
                  call read_and_scatter_se_2d(ncid, 'T2',     strt2, cnt2, INTP_T2(:,:,2))
                  call read_and_scatter_se_2d(ncid, 'TD2',    strt2, cnt2, INTP_TD2(:,:,2))
                  call read_and_scatter_se_2d(ncid, 'TS',     strt2, cnt2, INTP_TS(:,:,2))
                  call read_and_scatter_se_2d(ncid, 'PRECC',  strt2, cnt2, INTP_PRECC(:,:,2))
                  call read_and_scatter_se_2d(ncid, 'PRECL',  strt2, cnt2, INTP_PRECL(:,:,2))
                  call read_and_scatter_se_2d(ncid, 'PRECSC', strt2, cnt2, INTP_PRECSC(:,:,2))
                  call read_and_scatter_se_2d(ncid, 'PRECSL', strt2, cnt2, INTP_PRECSL(:,:,2))
                  call read_and_scatter_se_2d(ncid, 'ER',     strt2, cnt2, INTP_EVAP(:,:,2))
                  call read_and_scatter_se_2d(ncid, 'LHFLX',  strt2, cnt2, INTP_LHFLX(:,:,2))
                  call read_and_scatter_se_2d(ncid, 'SHFLX',  strt2, cnt2, INTP_SHFLX(:,:,2))
                  call read_and_scatter_se_2d(ncid, 'FSNS',   strt2, cnt2, INTP_FSNS(:,:,2))
                  call read_and_scatter_se_2d(ncid, 'FLDS',   strt2, cnt2, INTP_FLDS(:,:,2))
                  call read_and_scatter_se_2d(ncid, 'FSDS',   strt2, cnt2, INTP_FSDS(:,:,2))
                  call read_and_scatter_se_2d(ncid, 'FSDSD',  strt2, cnt2, INTP_FSDSD(:,:,2))
                  call read_and_scatter_se_2d(ncid, 'FSDSUV', strt2, cnt2, INTP_FSDSUV(:,:,2))

                  call t_startf ('read_nudging_data')
                  if (masterproc) then
                      istat=nf90_close(ncid1)
                      if (istat.ne.NF90_NOERR) then
                         write(iulog,*) nf90_strerror(istat)
                         call endrun ('UPDATE_LAND_SURFACE_ANALYSES_SE_CLOSE_LINEAR_NETCDF')
                      end if
                  end if
                  call t_stopf ('read_nudging_data')


              else
                  ! The previous end point becomes the start point
                  ! Only need to read in the new end point
                  !-----------------------------------------------
                  INTP_U10(:,:,1)    = INTP_U10(:,:,2)
                  INTP_V10(:,:,1)    = INTP_V10(:,:,2)
                  INTP_T2(:,:,1)     = INTP_T2(:,:,2)
                  INTP_TD2(:,:,1)    = INTP_TD2(:,:,2)
                  INTP_TS(:,:,1)     = INTP_TS(:,:,2)
                  INTP_PRECC(:,:,1)  = INTP_PRECC(:,:,2)
                  INTP_PRECL(:,:,1)  = INTP_PRECL(:,:,2)
                  INTP_PRECSC(:,:,1) = INTP_PRECSC(:,:,2)
                  INTP_PRECSL(:,:,1) = INTP_PRECSL(:,:,2)
                  INTP_EVAP(:,:,1)   = INTP_EVAP(:,:,2)
                  INTP_LHFLX(:,:,1)  = INTP_LHFLX(:,:,2)
                  INTP_SHFLX(:,:,1)  = INTP_SHFLX(:,:,2)
                  INTP_FSNS(:,:,1)   = INTP_FSNS(:,:,2)
                  INTP_FLDS(:,:,1)   = INTP_FLDS(:,:,2)
                  INTP_FSDS(:,:,1)   = INTP_FSDS(:,:,2)
                  INTP_FSDSD(:,:,1)  = INTP_FSDSD(:,:,2)
                  INTP_FSDSUV(:,:,1) = INTP_FSDSUV(:,:,2)

                  if (masterproc) then
                     strt2(1) = 1
                     strt2(2) = 1
                     cnt2(1)  = ncol
                     cnt2(2)  = 1
                  end if
                  call read_and_scatter_se_2d(ncid, 'U10',    strt2, cnt2, INTP_U10(:,:,2))
                  call read_and_scatter_se_2d(ncid, 'V10',    strt2, cnt2, INTP_V10(:,:,2))
                  call read_and_scatter_se_2d(ncid, 'T2',     strt2, cnt2, INTP_T2(:,:,2))
                  call read_and_scatter_se_2d(ncid, 'TD2',    strt2, cnt2, INTP_TD2(:,:,2))
                  call read_and_scatter_se_2d(ncid, 'TS',     strt2, cnt2, INTP_TS(:,:,2))
                  call read_and_scatter_se_2d(ncid, 'PRECC',  strt2, cnt2, INTP_PRECC(:,:,2))
                  call read_and_scatter_se_2d(ncid, 'PRECL',  strt2, cnt2, INTP_PRECL(:,:,2))
                  call read_and_scatter_se_2d(ncid, 'PRECSC', strt2, cnt2, INTP_PRECSC(:,:,2))
                  call read_and_scatter_se_2d(ncid, 'PRECSL', strt2, cnt2, INTP_PRECSL(:,:,2))
                  call read_and_scatter_se_2d(ncid, 'ER',     strt2, cnt2, INTP_EVAP(:,:,2))
                  call read_and_scatter_se_2d(ncid, 'LHFLX',  strt2, cnt2, INTP_LHFLX(:,:,2))
                  call read_and_scatter_se_2d(ncid, 'SHFLX',  strt2, cnt2, INTP_SHFLX(:,:,2))
                  call read_and_scatter_se_2d(ncid, 'FSNS',   strt2, cnt2, INTP_FSNS(:,:,2))
                  call read_and_scatter_se_2d(ncid, 'FLDS',   strt2, cnt2, INTP_FLDS(:,:,2))
                  call read_and_scatter_se_2d(ncid, 'FSDS',   strt2, cnt2, INTP_FSDS(:,:,2))
                  call read_and_scatter_se_2d(ncid, 'FSDSD',  strt2, cnt2, INTP_FSDSD(:,:,2))
                  call read_and_scatter_se_2d(ncid, 'FSDSUV', strt2, cnt2, INTP_FSDSUV(:,:,2))

              end if    ! first file for single time slice per file

           else
              ! Multiple time slices per file
              ! The start point uses the CURR time slice
              ! The end point uses the NEXT time slice
              !-----------------------------------------
              if (first_srf_file) then
                  first_srf_file = .false.
                  call get_curr_date(Year,Month,Day,Sec)
                  n_cnt = Sec/Nudge_Step + 1
                  if (n_cnt .eq. Nudge_SRF_File_Ntime) then
                     if (masterproc) then
                        call t_startf ('read_nudging_data')
                        call open_netcdf (ncid1, -Nudge_Step, Nudge_SRF_File_Template)
                        call t_stopf ('read_nudging_data')
                        strt2(1) = 1
                        strt2(2) = n_cnt
                        cnt2(1)  = ncol
                        cnt2(2)  = 1
                     end if

                     !-----------------------------------------
                     ! The start point uses the CURR time slice
                     !-----------------------------------------
                     call read_and_scatter_se_2d(ncid1, 'U10',    strt2, cnt2, INTP_U10(:,:,1))
                     call read_and_scatter_se_2d(ncid1, 'V10',    strt2, cnt2, INTP_V10(:,:,1))
                     call read_and_scatter_se_2d(ncid1, 'T2',     strt2, cnt2, INTP_T2(:,:,1))
                     call read_and_scatter_se_2d(ncid1, 'TD2',    strt2, cnt2, INTP_TD2(:,:,1))
                     call read_and_scatter_se_2d(ncid1, 'TS',     strt2, cnt2, INTP_TS(:,:,1))
                     call read_and_scatter_se_2d(ncid1, 'PRECC',  strt2, cnt2, INTP_PRECC(:,:,1))
                     call read_and_scatter_se_2d(ncid1, 'PRECL',  strt2, cnt2, INTP_PRECL(:,:,1))
                     call read_and_scatter_se_2d(ncid1, 'PRECSC', strt2, cnt2, INTP_PRECSC(:,:,1))
                     call read_and_scatter_se_2d(ncid1, 'PRECSL', strt2, cnt2, INTP_PRECSL(:,:,1))
                     call read_and_scatter_se_2d(ncid1, 'ER',     strt2, cnt2, INTP_EVAP(:,:,1))
                     call read_and_scatter_se_2d(ncid1, 'LHFLX',  strt2, cnt2, INTP_LHFLX(:,:,1))
                     call read_and_scatter_se_2d(ncid1, 'SHFLX',  strt2, cnt2, INTP_SHFLX(:,:,1))
                     call read_and_scatter_se_2d(ncid1, 'FSNS',   strt2, cnt2, INTP_FSNS(:,:,1))
                     call read_and_scatter_se_2d(ncid1, 'FLDS',   strt2, cnt2, INTP_FLDS(:,:,1))
                     call read_and_scatter_se_2d(ncid1, 'FSDS',   strt2, cnt2, INTP_FSDS(:,:,1))
                     call read_and_scatter_se_2d(ncid1, 'FSDSD',  strt2, cnt2, INTP_FSDSD(:,:,1))
                     call read_and_scatter_se_2d(ncid1, 'FSDSUV', strt2, cnt2, INTP_FSDSUV(:,:,1))

                     if (masterproc) then
                        strt2(2) = 1
                     end if

                     !-----------------------------------------
                     ! The end point uses the NEXT time slice
                     !---------------------------------------
                     call read_and_scatter_se_2d(ncid, 'U10',    strt2, cnt2, INTP_U10(:,:,2))
                     call read_and_scatter_se_2d(ncid, 'V10',    strt2, cnt2, INTP_V10(:,:,2))
                     call read_and_scatter_se_2d(ncid, 'T2',     strt2, cnt2, INTP_T2(:,:,2))
                     call read_and_scatter_se_2d(ncid, 'TD2',    strt2, cnt2, INTP_TD2(:,:,2))
                     call read_and_scatter_se_2d(ncid, 'TS',     strt2, cnt2, INTP_TS(:,:,2))
                     call read_and_scatter_se_2d(ncid, 'PRECC',  strt2, cnt2, INTP_PRECC(:,:,2))
                     call read_and_scatter_se_2d(ncid, 'PRECL',  strt2, cnt2, INTP_PRECL(:,:,2))
                     call read_and_scatter_se_2d(ncid, 'PRECSC', strt2, cnt2, INTP_PRECSC(:,:,2))
                     call read_and_scatter_se_2d(ncid, 'PRECSL', strt2, cnt2, INTP_PRECSL(:,:,2))
                     call read_and_scatter_se_2d(ncid, 'ER',     strt2, cnt2, INTP_EVAP(:,:,2))
                     call read_and_scatter_se_2d(ncid, 'LHFLX',  strt2, cnt2, INTP_LHFLX(:,:,2))
                     call read_and_scatter_se_2d(ncid, 'SHFLX',  strt2, cnt2, INTP_SHFLX(:,:,2))
                     call read_and_scatter_se_2d(ncid, 'FSNS',   strt2, cnt2, INTP_FSNS(:,:,2))
                     call read_and_scatter_se_2d(ncid, 'FLDS',   strt2, cnt2, INTP_FLDS(:,:,2))
                     call read_and_scatter_se_2d(ncid, 'FSDS',   strt2, cnt2, INTP_FSDS(:,:,2))
                     call read_and_scatter_se_2d(ncid, 'FSDSD',  strt2, cnt2, INTP_FSDSD(:,:,2))
                     call read_and_scatter_se_2d(ncid, 'FSDSUV', strt2, cnt2, INTP_FSDSUV(:,:,2))

                     call t_startf ('read_nudging_data')
                     if (masterproc) then
                         istat=nf90_close(ncid1)
                         if (istat.ne.NF90_NOERR) then
                            write(iulog,*) nf90_strerror(istat)
                            call endrun ('UPDATE_LAND_SURFACE_ANALYSES_SE_CLOSE_LINEAR_NETCDF')
                         end if
                     end if
                     call t_stopf ('read_nudging_data')
                  else
                     ! two time slices are in the same nudging data file
                     do n = n_cnt, n_cnt+1
                        if (masterproc) then
                           strt2(1) = 1
                           strt2(2) = n
                           cnt2(1)  = ncol
                           cnt2(2)  = 1
                        end if
                        ind = n - n_cnt + 1
                        call read_and_scatter_se_2d(ncid, 'U10',    strt2, cnt2, INTP_U10(:,:,ind))
                        call read_and_scatter_se_2d(ncid, 'V10',    strt2, cnt2, INTP_V10(:,:,ind))
                        call read_and_scatter_se_2d(ncid, 'T2',     strt2, cnt2, INTP_T2(:,:,ind))
                        call read_and_scatter_se_2d(ncid, 'TD2',    strt2, cnt2, INTP_TD2(:,:,ind))
                        call read_and_scatter_se_2d(ncid, 'TS',     strt2, cnt2, INTP_TS(:,:,ind))
                        call read_and_scatter_se_2d(ncid, 'PRECC',  strt2, cnt2, INTP_PRECC(:,:,ind))
                        call read_and_scatter_se_2d(ncid, 'PRECL',  strt2, cnt2, INTP_PRECL(:,:,ind))
                        call read_and_scatter_se_2d(ncid, 'PRECSC', strt2, cnt2, INTP_PRECSC(:,:,ind))
                        call read_and_scatter_se_2d(ncid, 'PRECSL', strt2, cnt2, INTP_PRECSL(:,:,ind))
                        call read_and_scatter_se_2d(ncid, 'ER',     strt2, cnt2, INTP_EVAP(:,:,ind))
                        call read_and_scatter_se_2d(ncid, 'LHFLX',  strt2, cnt2, INTP_LHFLX(:,:,ind))
                        call read_and_scatter_se_2d(ncid, 'SHFLX',  strt2, cnt2, INTP_SHFLX(:,:,ind))
                        call read_and_scatter_se_2d(ncid, 'FSNS',   strt2, cnt2, INTP_FSNS(:,:,ind))
                        call read_and_scatter_se_2d(ncid, 'FLDS',   strt2, cnt2, INTP_FLDS(:,:,ind))
                        call read_and_scatter_se_2d(ncid, 'FSDS',   strt2, cnt2, INTP_FSDS(:,:,ind))
                        call read_and_scatter_se_2d(ncid, 'FSDSD',  strt2, cnt2, INTP_FSDSD(:,:,ind))
                        call read_and_scatter_se_2d(ncid, 'FSDSUV', strt2, cnt2, INTP_FSDSUV(:,:,ind))

                     end do
                  end if

              else

                  call get_curr_date(Year,Month,Day,Sec)
                  n_cnt = Sec/Nudge_Step + 1
                  n_cnt = n_cnt + 1                ! open the nudging data at future model time step
                  if (n_cnt .gt. Nudge_SRF_File_Ntime) then
                      n_cnt = 1
                  end if
                  ! The previous end point becomes the start point
                  ! Only need to read in the new end point
                  !-----------------------------------------------
                  INTP_U10(:,:,1)    = INTP_U10(:,:,2)
                  INTP_V10(:,:,1)    = INTP_V10(:,:,2)
                  INTP_T2(:,:,1)     = INTP_T2(:,:,2)
                  INTP_TD2(:,:,1)    = INTP_TD2(:,:,2)
                  INTP_TS(:,:,1)     = INTP_TS(:,:,2)
                  INTP_PRECC(:,:,1)  = INTP_PRECC(:,:,2)
                  INTP_PRECL(:,:,1)  = INTP_PRECL(:,:,2)
                  INTP_PRECSC(:,:,1) = INTP_PRECSC(:,:,2)
                  INTP_PRECSL(:,:,1) = INTP_PRECSL(:,:,2)
                  INTP_EVAP(:,:,1)   = INTP_EVAP(:,:,2)
                  INTP_LHFLX(:,:,1)  = INTP_LHFLX(:,:,2)
                  INTP_SHFLX(:,:,1)  = INTP_SHFLX(:,:,2)
                  INTP_FSNS(:,:,1)   = INTP_FSNS(:,:,2)
                  INTP_FLDS(:,:,1)   = INTP_FLDS(:,:,2)
                  INTP_FSDS(:,:,1)   = INTP_FSDS(:,:,2)
                  INTP_FSDSD(:,:,1)  = INTP_FSDSD(:,:,2)
                  INTP_FSDSUV(:,:,1) = INTP_FSDSUV(:,:,2)

                  if (masterproc) then
                     strt2(1) = 1
                     strt2(2) = n_cnt
                     cnt2(1)  = ncol
                     cnt2(2)  = 1
                  end if

                  call read_and_scatter_se_2d(ncid, 'U10',    strt2, cnt2, INTP_U10(:,:,2))
                  call read_and_scatter_se_2d(ncid, 'V10',    strt2, cnt2, INTP_V10(:,:,2))
                  call read_and_scatter_se_2d(ncid, 'T2',     strt2, cnt2, INTP_T2(:,:,2))
                  call read_and_scatter_se_2d(ncid, 'TD2',    strt2, cnt2, INTP_TD2(:,:,2))
                  call read_and_scatter_se_2d(ncid, 'TS',     strt2, cnt2, INTP_TS(:,:,2))
                  call read_and_scatter_se_2d(ncid, 'PRECC',  strt2, cnt2, INTP_PRECC(:,:,2))
                  call read_and_scatter_se_2d(ncid, 'PRECL',  strt2, cnt2, INTP_PRECL(:,:,2))
                  call read_and_scatter_se_2d(ncid, 'PRECSC', strt2, cnt2, INTP_PRECSC(:,:,2))
                  call read_and_scatter_se_2d(ncid, 'PRECSL', strt2, cnt2, INTP_PRECSL(:,:,2))
                  call read_and_scatter_se_2d(ncid, 'ER',     strt2, cnt2, INTP_EVAP(:,:,2))
                  call read_and_scatter_se_2d(ncid, 'LHFLX',  strt2, cnt2, INTP_LHFLX(:,:,2))
                  call read_and_scatter_se_2d(ncid, 'SHFLX',  strt2, cnt2, INTP_SHFLX(:,:,2))
                  call read_and_scatter_se_2d(ncid, 'FSNS',   strt2, cnt2, INTP_FSNS(:,:,2))
                  call read_and_scatter_se_2d(ncid, 'FLDS',   strt2, cnt2, INTP_FLDS(:,:,2))
                  call read_and_scatter_se_2d(ncid, 'FSDS',   strt2, cnt2, INTP_FSDS(:,:,2))
                  call read_and_scatter_se_2d(ncid, 'FSDSD',  strt2, cnt2, INTP_FSDSD(:,:,2))
                  call read_and_scatter_se_2d(ncid, 'FSDSUV', strt2, cnt2, INTP_FSDSUV(:,:,2))

              end if ! first_srf_file for multiple time slices

           end if     ! single vs. multiple time slices per file

      case default
          write(iulog,*) 'ERROR: Unknown Input Nudge Method'
          call endrun('nudging_update_srf_analyses_se: bad input nudge method')
   end select

   if(masterproc) then
     istat=nf90_close(ncid)
     if (istat.ne.NF90_NOERR) then
        write(iulog,*) nf90_strerror(istat)
        call endrun ('UPDATE_LAND_SURFACE_ANALYSES_SE_CLOSE')
     end if
   end if ! (masterproc) then

   ! End Routine
   !------------
   return
  end subroutine ! nudging_update_srf_analyses_se

  !================================================================
  subroutine nudging_update_analyses_eul(anal_file)
   !
   ! NUDGING_UPDATE_ANALYSES_EUL:
   !                 Open the given analyses data file, read in
   !                 U,V,T,Q, and PS values and then distribute
   !                 the values to all of the chunks.
   !===============================================================
!   use wrap_nf
   use ppgrid ,only: pver
   use netcdf

   ! Arguments
   !-------------
   character(len=*),intent(in):: anal_file

   ! Local values
   !-------------
   integer lev
   integer nlon,nlat,plev,istat
   integer ncid,varid
   integer ilat,ilon,ilev
   real(r8) Xanal(Nudge_nlon,Nudge_nlat,Nudge_nlev)
   real(r8) PSanal(Nudge_nlon,Nudge_nlat)
   real(r8) Lat_anal(Nudge_nlat)
   real(r8) Lon_anal(Nudge_nlon)
   real(r8) Xtrans(Nudge_nlon,Nudge_nlev,Nudge_nlat)

   ! Check the existence of the analyses file; broadcast the file status to
   ! all the other MPI nodes. If the file is not there, then just return.
   !------------------------------------------------------------------------
   if(masterproc) then
     inquire(FILE=trim(anal_file),EXIST=Nudge_File_Present)
   endif
#ifdef SPMD
   call mpibcast(Nudge_File_Present, 1, mpilog, 0, mpicom)
#endif
   if(.not.Nudge_File_Present) return

   ! masterporc does all of the work here
   !-----------------------------------------
   if(masterproc) then

     ! Open the given file
     !-----------------------
     istat=nf90_open(trim(anal_file),NF90_NOWRITE,ncid)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*)'NF90_OPEN: failed for file ',trim(anal_file)
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_EUL')
     endif

     ! Read in Dimensions
     !--------------------
!     call wrap_inq_dimid (ncid,'lon',varid)
!     call wrap_inq_dimlen(ncid,varid,nlon)
     istat=nf90_inq_dimid(ncid,'lon',varid)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_EUL')
     endif
     istat=nf90_inquire_dimension(ncid,varid,len=nlon)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_EUL')
     endif

!     call wrap_inq_dimid (ncid,'lat',varid)
!     call wrap_inq_dimlen(ncid,varid,nlat)
     istat=nf90_inq_dimid(ncid,'lat',varid)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_EUL')
     endif
     istat=nf90_inquire_dimension(ncid,varid,len=nlat)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_EUL')
     endif

!     call wrap_inq_dimid (ncid,'lev',varid)
!     call wrap_inq_dimlen(ncid,varid,plev)
     istat=nf90_inq_dimid(ncid,'lev',varid)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_EUL')
     endif
     istat=nf90_inquire_dimension(ncid,varid,len=plev)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_EUL')
     endif

!     call wrap_inq_varid(ncid,'lon',varid)
!     call wrap_get_var_realx(ncid,varid,Lon_anal)
     istat=nf90_inq_varid(ncid,'lon',varid)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_EUL')
     endif
     istat=nf90_get_var(ncid,varid,Lon_anal)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_EUL')
     endif

!     call wrap_inq_varid(ncid,'lat',varid)
!     call wrap_get_var_realx(ncid,varid,Lat_anal)
     istat=nf90_inq_varid(ncid,'lat',varid)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_EUL')
     endif
     istat=nf90_get_var(ncid,varid,Lat_anal)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_EUL')
     endif

     if((Nudge_nlon.ne.nlon).or.(Nudge_nlat.ne.nlat).or.(plev.ne.pver)) then
      write(iulog,*) 'ERROR: nudging_update_analyses_eul: nlon=',nlon,' Nudge_nlon=',Nudge_nlon
      write(iulog,*) 'ERROR: nudging_update_analyses_eul: nlat=',nlat,' Nudge_nlat=',Nudge_nlat
      write(iulog,*) 'ERROR: nudging_update_analyses_eul: plev=',plev,' pver=',pver
      call endrun('nudging_update_analyses_eul: analyses dimension mismatch')
     endif

     ! Read in, transpose lat/lev indices,
     ! and scatter data arrays
     !----------------------------------
!     call wrap_inq_varid    (ncid,'U',varid)
!     call wrap_get_var_realx(ncid,varid,Xanal)
     istat=nf90_inq_varid(ncid,'U',varid)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_EUL')
     endif
     istat=nf90_get_var(ncid,varid,Xanal)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_EUL')
     endif
     do ilat=1,nlat
     do ilev=1,plev
     do ilon=1,nlon
       Xtrans(ilon,ilev,ilat)=Xanal(ilon,ilat,ilev)
     end do
     end do
     end do
   endif ! (masterproc) then
   call scatter_field_to_chunk(1,Nudge_nlev,1,Nudge_nlon,Xtrans ,Target_U)

   if(masterproc) then
!     call wrap_inq_varid    (ncid,'V',varid)
!     call wrap_get_var_realx(ncid,varid,Xanal)
     istat=nf90_inq_varid(ncid,'V',varid)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_EUL')
     endif
     istat=nf90_get_var(ncid,varid,Xanal)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_EUL')
     endif
     do ilat=1,nlat
     do ilev=1,plev
     do ilon=1,nlon
       Xtrans(ilon,ilev,ilat)=Xanal(ilon,ilat,ilev)
     end do
     end do
     end do
   endif ! (masterproc) then
   call scatter_field_to_chunk(1,Nudge_nlev,1,Nudge_nlon,Xtrans ,Target_V)

   if(masterproc) then
!     call wrap_inq_varid    (ncid,'T',varid)
!     call wrap_get_var_realx(ncid,varid,Xanal)
     istat=nf90_inq_varid(ncid,'T',varid)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_EUL')
     endif
     istat=nf90_get_var(ncid,varid,Xanal)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_EUL')
     endif
     do ilat=1,nlat
     do ilev=1,plev
     do ilon=1,nlon
       Xtrans(ilon,ilev,ilat)=Xanal(ilon,ilat,ilev)
     end do
     end do
     end do
   endif ! (masterproc) then
   call scatter_field_to_chunk(1,Nudge_nlev,1,Nudge_nlon,Xtrans ,Target_T)

   if(masterproc) then
!     call wrap_inq_varid    (ncid,'Q',varid)
!     call wrap_get_var_realx(ncid,varid,Xanal)
     istat=nf90_inq_varid(ncid,'Q',varid)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_EUL')
     endif
     istat=nf90_get_var(ncid,varid,Xanal)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_EUL')
     endif
     do ilat=1,nlat
     do ilev=1,plev
     do ilon=1,nlon
       Xtrans(ilon,ilev,ilat)=Xanal(ilon,ilat,ilev)
     end do
     end do
     end do
   endif ! (masterproc) then
   call scatter_field_to_chunk(1,Nudge_nlev,1,Nudge_nlon,Xtrans ,Target_Q)

   if(masterproc) then
!!    call wrap_inq_varid    (ncid,'PS',varid)
!!    call wrap_get_var_realx(ncid,varid,PSanal)
!    istat=nf90_inq_varid(ncid,'PS',varid)
!    if(istat.ne.NF90_NOERR) then
!      write(iulog,*) nf90_strerror(istat)
!      call endrun ('UPDATE_ANALYSES_SE')
!    endif
!    istat=nf90_get_var(ncid,varid,PSanal)
!    if(istat.ne.NF90_NOERR) then
!      write(iulog,*) nf90_strerror(istat)
!      call endrun ('UPDATE_ANALYSES_SE')
!    endif

     ! Close the analyses file
     !-----------------------
!     call wrap_close(ncid)
     istat=nf90_close(ncid)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_EUL')
     endif
   endif ! (masterproc) then
!  call scatter_field_to_chunk(1,         1,1,Nudge_nlon,PSanal,Target_PS)

   ! End Routine
   !------------
   return
  end subroutine ! nudging_update_analyses_eul
  !================================================================


  !================================================================
  subroutine nudging_update_analyses_fv(anal_file)
   !
   ! NUDGING_UPDATE_ANALYSES_FV:
   !                 Open the given analyses data file, read in
   !                 U,V,T,Q, and PS values and then distribute
   !                 the values to all of the chunks.
   !===============================================================
!   use wrap_nf
   use ppgrid ,only: pver
   use netcdf
   use filenames ,only: interpret_filename_spec

   ! Arguments
   !-------------
   character(len=*),intent(in):: anal_file

   ! Local values
   !-------------
   integer lev
   integer nlon,nlat,plev,istat
   integer ncid,varid
   integer ilat,ilon,ilev
   real(r8) Xanal(Nudge_nlon,Nudge_nlat,Nudge_nlev)
   real(r8) PSanal(Nudge_nlon,Nudge_nlat)
   real(r8) Lat_anal(Nudge_nlat)
   real(r8) Lon_anal(Nudge_nlon)
   real(r8) Xtrans(Nudge_nlon,Nudge_nlev,Nudge_nlat)
!DIAG
   real(r8) Uanal(Nudge_nlon,Nudge_slat,Nudge_nlev)
!DIAG

   integer :: cnt4(4)               ! array of counts for each dimension
   integer :: strt4(4)              ! array of starting indices
   integer :: n, n_cnt, ncid1, ind
   integer :: Year, Month, Day, Sec

   ! Check the existence of the analyses file; broadcast the file status to
   ! all the other MPI nodes. If the file is not there, then just return.
   !------------------------------------------------------------------------
   if(masterproc) then
     inquire(FILE=trim(anal_file),EXIST=Nudge_File_Present)
   endif
#ifdef SPMD
   call mpibcast(Nudge_File_Present, 1, mpilog, 0, mpicom)
#endif
   if (trim(Nudge_Method).eq. "MLTBC") then
     if(masterproc) then 
        write(iulog,*) 'Warning: using Machine Learning model to predict nudging tendency'  
        write(iulog,*) 'Warning: No need to read reference data, return'
     end if 
     return       
   else
     if(.not.Nudge_File_Present) return
   end if

   ! masterporc does all of the work here
   !-----------------------------------------
   if(masterproc) then

     ! Open the given file
     !-----------------------
     istat=nf90_open(trim(anal_file),NF90_NOWRITE,ncid)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*)'NF90_OPEN: failed for file ',trim(anal_file)
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_FV')
     endif

     ! Read in Dimensions
     !--------------------
!     call wrap_inq_dimid (ncid,'lon',varid)
!     call wrap_inq_dimlen(ncid,varid,nlon)
     istat=nf90_inq_dimid(ncid,'lon',varid)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_FV')
     endif
     istat=nf90_inquire_dimension(ncid,varid,len=nlon)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_FV')
     endif

!     call wrap_inq_dimid (ncid,'lat',varid)
!     call wrap_inq_dimlen(ncid,varid,nlat)
     istat=nf90_inq_dimid(ncid,'lat',varid)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_FV')
     endif
     istat=nf90_inquire_dimension(ncid,varid,len=nlat)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_FV')
     endif

!     call wrap_inq_dimid (ncid,'lev',varid)
!     call wrap_inq_dimlen(ncid,varid,plev)
     istat=nf90_inq_dimid(ncid,'lev',varid)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_FV')
     endif
     istat=nf90_inquire_dimension(ncid,varid,len=plev)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_FV')
     endif

!     call wrap_inq_varid(ncid,'lon',varid)
!     call wrap_get_var_realx(ncid,varid,Lon_anal)
     istat=nf90_inq_varid(ncid,'lon',varid)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_FV')
     endif
     istat=nf90_get_var(ncid,varid,Lon_anal)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_FV')
     endif

!     call wrap_inq_varid(ncid,'lat',varid)
!     call wrap_get_var_realx(ncid,varid,Lat_anal)
     istat=nf90_inq_varid(ncid,'lat',varid)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_FV')
     endif
     istat=nf90_get_var(ncid,varid,Lat_anal)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_FV')
     endif

     if((Nudge_nlon.ne.nlon).or.(Nudge_nlat.ne.nlat).or.(plev.ne.pver)) then
      write(iulog,*) 'ERROR: nudging_update_analyses_fv: nlon=',nlon,' Nudge_nlon=',Nudge_nlon
      write(iulog,*) 'ERROR: nudging_update_analyses_fv: nlat=',nlat,' Nudge_nlat=',Nudge_nlat
      write(iulog,*) 'ERROR: nudging_update_analyses_fv: plev=',plev,' pver=',pver
      call endrun('nudging_update_analyses_fv: analyses dimension mismatch')
     endif

   end if ! (masterproc) then

   select case (Nudge_Method)
      case ('Step')
           if (masterproc) then
              call get_curr_date(Year,Month,Day,Sec)
              n_cnt = Sec/Nudge_Step + 1
              n_cnt = n_cnt + 1                      ! nudge to future model time step
              if (n_cnt .gt. Nudge_File_Ntime) then
                  n_cnt = 1
              end if
              strt4(1) = 1
              strt4(2) = 1
              strt4(3) = 1
              strt4(4) = n_cnt
              cnt4(1)  = nlon
              cnt4(2)  = nlat
              cnt4(3)  = pver
              cnt4(4)  = 1
           end if
           call read_and_scatter_fv(ncid, 'U', strt4, cnt4, Target_U)
           call read_and_scatter_fv(ncid, 'V', strt4, cnt4, Target_V)
           call read_and_scatter_fv(ncid, 'T', strt4, cnt4, Target_T)
           call read_and_scatter_fv(ncid, 'Q', strt4, cnt4, Target_Q)

      case ('IMT')
           call get_curr_date(Year,Month,Day,Sec)
           if (mod(Sec,Nudge_Step) .ne. 0) return  ! ensure that intermittent simulations work for restart run
           if (masterproc) then
              call t_startf ('read_nudging_data')
              call open_netcdf (ncid1, -Nudge_Step, Nudge_File_Template)
              call t_stopf ('read_nudging_data')
              n_cnt = Sec/Nudge_Step + 1
              if (n_cnt .gt. Nudge_File_Ntime) then       ! account for one time slice per file
                  n_cnt = 1
              end if
              strt4(1) = 1
              strt4(2) = 1
              strt4(3) = 1
              strt4(4) = n_cnt
              cnt4(1)  = nlon
              cnt4(2)  = nlat
              cnt4(3)  = pver
              cnt4(4)  = 1
           end if
           ! Use the CURR time slice for nudging
           !------------------------------------
           call read_and_scatter_fv(ncid1, 'U', strt4, cnt4, Target_U)
           call read_and_scatter_fv(ncid1, 'V', strt4, cnt4, Target_V)
           call read_and_scatter_fv(ncid1, 'T', strt4, cnt4, Target_T)
           call read_and_scatter_fv(ncid1, 'Q', strt4, cnt4, Target_Q)

      case ('Linear')
           ! Single time slice per file
           ! Need to open a new netcdf file to get the CURR time slice
           if (Nudge_File_Ntime .eq. 1) then
              if (first_file) then
                  first_file = .false.
                  if (masterproc) then
                     call t_startf ('read_nudging_data')
                     call open_netcdf (ncid1, -Nudge_Step, Nudge_File_Template)
                     call t_stopf ('read_nudging_data')
                     strt4(1) = 1
                     strt4(2) = 1
                     strt4(3) = 1
                     strt4(4) = 1 
                     cnt4(1)  = nlon
                     cnt4(2)  = nlat
                     cnt4(3)  = pver
                     cnt4(4)  = 1
                  end if

                  !---------------------------------------
                  ! The start point uses the CURR time slice
                  !-----------------------------------------
                  call read_and_scatter_fv(ncid1, 'U', strt4, cnt4, INTP_U(:,:,:,1))
                  call read_and_scatter_fv(ncid1, 'V', strt4, cnt4, INTP_V(:,:,:,1))
                  call read_and_scatter_fv(ncid1, 'T', strt4, cnt4, INTP_T(:,:,:,1))
                  call read_and_scatter_fv(ncid1, 'Q', strt4, cnt4, INTP_Q(:,:,:,1))

                  !---------------------------------------
                  ! The end point uses the NEXT time slice
                  !---------------------------------------
                  call read_and_scatter_fv(ncid, 'U', strt4, cnt4, INTP_U(:,:,:,2))
                  call read_and_scatter_fv(ncid, 'V', strt4, cnt4, INTP_V(:,:,:,2))
                  call read_and_scatter_fv(ncid, 'T', strt4, cnt4, INTP_T(:,:,:,2))
                  call read_and_scatter_fv(ncid, 'Q', strt4, cnt4, INTP_Q(:,:,:,2))

                  call t_startf ('read_nudging_data')
                  if (masterproc) then
                      istat=nf90_close(ncid1)
                      if (istat.ne.NF90_NOERR) then
                         write(iulog,*) nf90_strerror(istat)
                         call endrun ('UPDATE_ANALYSES_FV_CLOSE_LINEAR_NETCDF')
                      end if
                  end if
                  call t_stopf ('read_nudging_data')

              else
                  !-----------------------------------------------
                  ! The previous end point becomes the start point
                  ! Only need to read in the new end point
                  !-----------------------------------------------
                  INTP_U(:,:,:,1) = INTP_U(:,:,:,2)
                  INTP_V(:,:,:,1) = INTP_V(:,:,:,2)
                  INTP_Q(:,:,:,1) = INTP_Q(:,:,:,2)
                  INTP_T(:,:,:,1) = INTP_T(:,:,:,2)
                  if (masterproc) then
                     strt4(1) = 1
                     strt4(2) = 1
                     strt4(3) = 1
                     strt4(4) = 1
                     cnt4(1)  = nlon
                     cnt4(2)  = nlat
                     cnt4(3)  = pver
                     cnt4(4)  = 1
                  end if
                  call read_and_scatter_fv(ncid, 'U', strt4, cnt4, INTP_U(:,:,:,2))
                  call read_and_scatter_fv(ncid, 'V', strt4, cnt4, INTP_V(:,:,:,2))
                  call read_and_scatter_fv(ncid, 'T', strt4, cnt4, INTP_T(:,:,:,2))
                  call read_and_scatter_fv(ncid, 'Q', strt4, cnt4, INTP_Q(:,:,:,2))

              end if    ! first file for single time slice per file

           else
              ! Multiple time slices per file
              ! The start point uses the CURR time slice
              ! The end point uses the NEXT time slice
              !-----------------------------------------
              if (first_file) then
                  first_file = .false.
                  call get_curr_date(Year,Month,Day,Sec)
                  n_cnt = Sec/Nudge_Step + 1
                  if (n_cnt .eq. Nudge_File_Ntime) then
                     if (masterproc) then
                        call t_startf ('read_nudging_data')
                        call open_netcdf (ncid1, -Nudge_Step, Nudge_File_Template)
                        call t_stopf ('read_nudging_data')
                        strt4(1) = 1
                        strt4(2) = 1
                        strt4(3) = 1
                        strt4(4) = n_cnt
                        cnt4(1)  = nlon
                        cnt4(2)  = nlat
                        cnt4(3)  = pver
                        cnt4(4)  = 1
                     end if

                     !-----------------------------------------
                     ! The start point uses the CURR time slice
                     !-----------------------------------------
                     call read_and_scatter_fv(ncid1, 'U', strt4, cnt4, INTP_U(:,:,:,1))
                     call read_and_scatter_fv(ncid1, 'V', strt4, cnt4, INTP_V(:,:,:,1))
                     call read_and_scatter_fv(ncid1, 'T', strt4, cnt4, INTP_T(:,:,:,1))
                     call read_and_scatter_fv(ncid1, 'Q', strt4, cnt4, INTP_Q(:,:,:,1))

                     if (masterproc) then
                        strt4(4) = 1
                     end if

                     !-----------------------------------------
                     ! The end point uses the NEXT time slice
                     !---------------------------------------
                     call read_and_scatter_fv(ncid, 'U', strt4, cnt4, INTP_U(:,:,:,2))
                     call read_and_scatter_fv(ncid, 'V', strt4, cnt4, INTP_V(:,:,:,2))
                     call read_and_scatter_fv(ncid, 'T', strt4, cnt4, INTP_T(:,:,:,2))
                     call read_and_scatter_fv(ncid, 'Q', strt4, cnt4, INTP_Q(:,:,:,2))

                     call t_startf ('read_nudging_data')
                     if (masterproc) then
                         istat=nf90_close(ncid1)
                         if (istat.ne.NF90_NOERR) then
                            write(iulog,*) nf90_strerror(istat)
                            call endrun ('UPDATE_ANALYSES_SE_CLOSE_LINEAR_NETCDF')
                         end if
                     end if
                     call t_stopf ('read_nudging_data')
                  else
                     ! two time slices are in the same nudging data file
                     do n = n_cnt, n_cnt+1
                        if (masterproc) then
                           strt4(1) = 1
                           strt4(2) = 1
                           strt4(3) = 1
                           strt4(4) = n
                           cnt4(1)  = nlon
                           cnt4(2)  = nlat
                           cnt4(3)  = pver
                           cnt4(4)  = 1
                        end if
                        ind = n - n_cnt + 1
                        call read_and_scatter_fv(ncid, 'U', strt4, cnt4, INTP_U(:,:,:,ind))
                        call read_and_scatter_fv(ncid, 'V', strt4, cnt4, INTP_V(:,:,:,ind))
                        call read_and_scatter_fv(ncid, 'T', strt4, cnt4, INTP_T(:,:,:,ind))
                        call read_and_scatter_fv(ncid, 'Q', strt4, cnt4, INTP_Q(:,:,:,ind))
                     end do
                  end if

              else

                  call get_curr_date(Year,Month,Day,Sec)
                  n_cnt = Sec/Nudge_Step + 1
                  n_cnt = n_cnt + 1                ! open the nudging data at future model time step
                  if (n_cnt .gt. Nudge_File_Ntime) then
                      n_cnt = 1
                  end if

                  ! The previous end point becomes the start point
                  ! Only need to read in the new end point
                  !-----------------------------------------------
                  INTP_U(:,:,:,1) = INTP_U(:,:,:,2)
                  INTP_V(:,:,:,1) = INTP_V(:,:,:,2)
                  INTP_Q(:,:,:,1) = INTP_Q(:,:,:,2)
                  INTP_T(:,:,:,1) = INTP_T(:,:,:,2)
                  if (masterproc) then
                     strt4(1)        = 1
                     strt4(2)        = 1
                     strt4(3)        = 1
                     strt4(4)        = n_cnt
                     cnt4(1)         = nlon
                     cnt4(2)         = nlat
                     cnt4(3)         = pver
                     cnt4(4)         = 1
                  end if

                  call read_and_scatter_fv(ncid, 'U', strt4, cnt4, INTP_U(:,:,:,2))
                  call read_and_scatter_fv(ncid, 'V', strt4, cnt4, INTP_V(:,:,:,2))
                  call read_and_scatter_fv(ncid, 'T', strt4, cnt4, INTP_T(:,:,:,2))
                  call read_and_scatter_fv(ncid, 'Q', strt4, cnt4, INTP_Q(:,:,:,2))

              end if ! first_file for multiple time slices

           end if     ! single vs. multiple time slices per file

      case default
          write(iulog,*) 'ERROR: Unknown Input Nudge Method'
          call endrun('nudging_update_analyses_fv: bad input nudge method')
   end select

!! JS 04/11/2020 - Comment out the old versions for FV dycore, 
!!                 which uses US and VS instead U and V.
!!     ! Read in, transpose lat/lev indices,
!!     ! and scatter data arrays
!!     !----------------------------------
!!!DIAG:  Dont have U, so jam US into U so tests can proceed:
!!!DIAG     call wrap_inq_varid    (ncid,'U',varid)
!!!DIAG     call wrap_get_var_realx(ncid,varid,Xanal)
!!!DIAG     do ilat=1,nlat
!!!DIAG     do ilev=1,plev
!!!DIAG     do ilon=1,nlon
!!!DIAG       Xtrans(ilon,ilev,ilat)=Xanal(ilon,ilat,ilev)
!!!DIAG     end do
!!!DIAG     end do
!!!DIAG     end do
!!!     call wrap_inq_varid    (ncid,'US',varid)
!!!     call wrap_get_var_realx(ncid,varid,Uanal)
!!     istat=nf90_inq_varid(ncid,'US',varid)
!!     if(istat.ne.NF90_NOERR) then
!!       write(iulog,*) nf90_strerror(istat)
!!       call endrun ('UPDATE_ANALYSES_FV')
!!     endif
!!     istat=nf90_get_var(ncid,varid,Uanal)
!!     if(istat.ne.NF90_NOERR) then
!!       write(iulog,*) nf90_strerror(istat)
!!       call endrun ('UPDATE_ANALYSES_FV')
!!     endif
!!     do ilat=1,(nlat-1)
!!     do ilev=1,plev
!!     do ilon=1,nlon
!!       Xtrans(ilon,ilev,ilat)=Uanal(ilon,ilat,ilev)
!!     end do
!!     end do
!!     end do
!!     Xtrans(:,:,ilat)=Xtrans(:,:,ilat-1)
!!   endif ! (masterproc) then
!!   call scatter_field_to_chunk(1,Nudge_nlev,1,Nudge_nlon,Xtrans ,Target_U)
!!
!!   if(masterproc) then
!!!DIAG:  Dont have V, so jam VS into V so tests can proceed:
!!!DIAG     call wrap_inq_varid    (ncid,'V',varid)
!!!DIAG     call wrap_get_var_realx(ncid,varid,Xanal)
!!!DIAG     do ilat=1,nlat
!!!DIAG     do ilev=1,plev
!!!DIAG     do ilon=1,nlon
!!!DIAG       Xtrans(ilon,ilev,ilat)=Xanal(ilon,ilat,ilev)
!!!DIAG     end do
!!!DIAG     end do
!!!DIAG     end do
!!!     call wrap_inq_varid    (ncid,'VS',varid)
!!!     call wrap_get_var_realx(ncid,varid,Xanal)
!!     istat=nf90_inq_varid(ncid,'VS',varid)
!!     if(istat.ne.NF90_NOERR) then
!!       write(iulog,*) nf90_strerror(istat)
!!       call endrun ('UPDATE_ANALYSES_FV')
!!     endif
!!     istat=nf90_get_var(ncid,varid,Xanal)
!!     if(istat.ne.NF90_NOERR) then
!!       write(iulog,*) nf90_strerror(istat)
!!       call endrun ('UPDATE_ANALYSES_FV')
!!     endif
!!     do ilat=1,nlat
!!     do ilev=1,plev
!!     do ilon=1,nlon
!!       Xtrans(ilon,ilev,ilat)=Xanal(ilon,ilat,ilev)
!!     end do
!!     end do
!!     end do
!!   endif ! (masterproc) then
!!   call scatter_field_to_chunk(1,Nudge_nlev,1,Nudge_nlon,Xtrans ,Target_V)
!!
!!   if(masterproc) then
!!!     call wrap_inq_varid    (ncid,'T',varid)
!!!     call wrap_get_var_realx(ncid,varid,Xanal)
!!     istat=nf90_inq_varid(ncid,'T',varid)
!!     if(istat.ne.NF90_NOERR) then
!!       write(iulog,*) nf90_strerror(istat)
!!       call endrun ('UPDATE_ANALYSES_FV')
!!     endif
!!     istat=nf90_get_var(ncid,varid,Xanal)
!!     if(istat.ne.NF90_NOERR) then
!!       write(iulog,*) nf90_strerror(istat)
!!       call endrun ('UPDATE_ANALYSES_FV')
!!     endif
!!     do ilat=1,nlat
!!     do ilev=1,plev
!!     do ilon=1,nlon
!!       Xtrans(ilon,ilev,ilat)=Xanal(ilon,ilat,ilev)
!!     end do
!!     end do
!!     end do
!!   endif ! (masterproc) then
!!   call scatter_field_to_chunk(1,Nudge_nlev,1,Nudge_nlon,Xtrans ,Target_T)
!!
!!   if(masterproc) then
!!!     call wrap_inq_varid    (ncid,'Q',varid)
!!!     call wrap_get_var_realx(ncid,varid,Xanal)
!!     istat=nf90_inq_varid(ncid,'Q',varid)
!!     if(istat.ne.NF90_NOERR) then
!!       write(iulog,*) nf90_strerror(istat)
!!       call endrun ('UPDATE_ANALYSES_FV')
!!     endif
!!     istat=nf90_get_var(ncid,varid,Xanal)
!!     if(istat.ne.NF90_NOERR) then
!!       write(iulog,*) nf90_strerror(istat)
!!       call endrun ('UPDATE_ANALYSES_FV')
!!     endif
!!     do ilat=1,nlat
!!     do ilev=1,plev
!!     do ilon=1,nlon
!!       Xtrans(ilon,ilev,ilat)=Xanal(ilon,ilat,ilev)
!!     end do
!!     end do
!!     end do
!!   endif ! (masterproc) then
!!   call scatter_field_to_chunk(1,Nudge_nlev,1,Nudge_nlon,Xtrans ,Target_Q)

   if(masterproc) then
!!    call wrap_inq_varid    (ncid,'PS',varid)
!!    call wrap_get_var_realx(ncid,varid,PSanal)
!    istat=nf90_inq_varid(ncid,'PS',varid)
!    if(istat.ne.NF90_NOERR) then
!      write(iulog,*) nf90_strerror(istat)
!      call endrun ('UPDATE_ANALYSES_SE')
!    endif
!    istat=nf90_get_var(ncid,varid,PSanal)
!    if(istat.ne.NF90_NOERR) then
!      write(iulog,*) nf90_strerror(istat)
!      call endrun ('UPDATE_ANALYSES_SE')
!    endif

     ! Close the analyses file
     !-----------------------
!     call wrap_close(ncid)
     istat=nf90_close(ncid)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_FV')
     endif
   endif ! (masterproc) then
!  call scatter_field_to_chunk(1,         1,1,Nudge_nlon,PSanal,Target_PS)

   ! End Routine
   !------------
   return
  end subroutine ! nudging_update_analyses_fv
  !================================================================


  !================================================================
  subroutine nudging_set_profile(rlat,rlon,Nudge_prof,Wprof,nlev)
   !
   ! NUDGING_SET_PROFILE: for the given lat,lon, and Nudging_prof, set
   !                      the verical profile of window coeffcients.
   !                      Values range from 0. to 1. to affect spatial
   !                      variations on nudging strength.
   !===============================================================

   ! Arguments
   !--------------
   integer  nlev,Nudge_prof
   real(r8) rlat,rlon
   real(r8) Wprof(nlev)

   ! Local values
   !----------------
   integer  ilev
   real(r8) Hcoef,latx,lonx,Vmax,Vmin
   real(r8) lon_lo,lon_hi,lat_lo,lat_hi,lev_lo,lev_hi

   !---------------
   ! set coeffcient
   !---------------
   if(Nudge_prof.eq.0) then
     ! No Nudging
     !-------------
     Wprof(:)=0.0
   elseif(Nudge_prof.eq.1) then
     ! Uniform Nudging
     !-----------------
     Wprof(:)=1.0
   elseif(Nudge_prof.eq.2) then
     ! Localized Nudging with specified Heaviside window function
     !------------------------------------------------------------
     if(Nudge_Hwin_max.le.Nudge_Hwin_min) then
       ! For a constant Horizontal window function,
       ! just set Hcoef to the maximum of Hlo/Hhi.
       !--------------------------------------------
       Hcoef=max(Nudge_Hwin_lo,Nudge_Hwin_hi)
     else
       ! get lat/lon relative to window center
       !------------------------------------------
       latx=rlat-Nudge_Hwin_lat0
       lonx=rlon-Nudge_Hwin_lon0
       if(lonx.gt. 180.) lonx=lonx-360.
       if(lonx.le.-180.) lonx=lonx+360.

       ! Calcualte RAW window value
       !-------------------------------
       lon_lo=(Nudge_Hwin_lonWidthH+lonx)/Nudge_Hwin_lonDelta
       lon_hi=(Nudge_Hwin_lonWidthH-lonx)/Nudge_Hwin_lonDelta
       lat_lo=(Nudge_Hwin_latWidthH+latx)/Nudge_Hwin_latDelta
       lat_hi=(Nudge_Hwin_latWidthH-latx)/Nudge_Hwin_latDelta
       Hcoef=((1.+tanh(lon_lo))/2.)*((1.+tanh(lon_hi))/2.) &
            *((1.+tanh(lat_lo))/2.)*((1.+tanh(lat_hi))/2.)

       ! Scale the horizontal window coef for specfied range of values.
       !--------------------------------------------------------
       Hcoef=(Hcoef-Nudge_Hwin_min)/(Nudge_Hwin_max-Nudge_Hwin_min)
       Hcoef=(1.-Hcoef)*Nudge_Hwin_lo + Hcoef*Nudge_Hwin_hi
     endif

     ! Load the RAW vertical window
     !------------------------------
     do ilev=1,nlev
       lev_lo=(float(ilev)-Nudge_Vwin_Lindex)/Nudge_Vwin_Ldelta
       lev_hi=(Nudge_Vwin_Hindex-float(ilev))/Nudge_Vwin_Hdelta
       Wprof(ilev)=((1.+tanh(lev_lo))/2.)*((1.+tanh(lev_hi))/2.)
     end do

     ! Scale the Window function to span the values between Vlo and Vhi:
     !-----------------------------------------------------------------
     Vmax=maxval(Wprof)
     Vmin=minval(Wprof)
     if(Vmax.le.Vmin) then
       ! For a constant Vertical window function,
       ! load maximum of Vlo/Vhi into Wprof()
       !--------------------------------------------
       Vmax=max(Nudge_Vwin_lo,Nudge_Vwin_hi)
       Wprof(:)=Vmax
     else
       ! Scale the RAW vertical window for specfied range of values.
       !--------------------------------------------------------
       Wprof(:)=(Wprof(:)-Vmin)/(Vmax-Vmin)
       Wprof(:)=Nudge_Vwin_lo + Wprof(:)*(Nudge_Vwin_hi-Nudge_Vwin_lo)
     endif

     ! The desired result is the product of the vertical profile
     ! and the horizontal window coeffcient.
     !----------------------------------------------------
     Wprof(:)=Hcoef*Wprof(:)
   else
     call endrun('nudging_set_profile:: Unknown Nudge_prof value')
   endif

   ! End Routine
   !------------
   return
  end subroutine ! nudging_set_profile
  !================================================================

  !================================================================
  real(r8) function nudging_set_PSprofile(rlat,rlon,Nudge_PSprof)
   !
   ! NUDGING_SET_PSPROFILE: for the given lat and lon set the surface
   !                      pressure profile value for the specified index.
   !                      Values range from 0. to 1. to affect spatial
   !                      variations on nudging strength.
   !===============================================================

   ! Arguments
   !--------------
   real(r8) rlat,rlon
   integer  Nudge_PSprof

   ! Local values
   !----------------

   !---------------
   ! set coeffcient
   !---------------
   if(Nudge_PSprof.eq.0) then
     ! No Nudging
     !-------------
     nudging_set_PSprofile=0.0
   elseif(Nudge_PSprof.eq.1) then
     ! Uniform Nudging
     !-----------------
     nudging_set_PSprofile=1.0
   else
     call endrun('nudging_set_PSprofile:: Unknown Nudge_prof value')
   endif

   ! End Routine
   !------------
   return
  end function ! nudging_set_PSprofile
  !================================================================

  !================================================================
  real(r8) function nudging_set_SRFprofile(rlat,rlon,Nudge_SRFprof)
   !
   ! NUDGING_SET_SRFPROFILE: for the given lat and lon set the land 
   !                      surface profile value for the specified index.
   !                      Values range from 0. to 1. to affect spatial
   !                      variations on nudging strength.
   !===============================================================

   ! Arguments
   !--------------
   real(r8) rlat,rlon
   integer  Nudge_SRFprof

   ! Local values
   !----------------

   !---------------
   ! set coeffcient
   !---------------
   if(Nudge_SRFprof.eq.0) then
     ! No Nudging
     !-------------
     nudging_set_SRFprofile=0.0
   elseif(Nudge_SRFprof.eq.1) then
     ! Uniform Nudging
     !-----------------
     nudging_set_SRFprofile=1.0
   else
     call endrun('nudging_set_SRFprofile:: Unknown Nudge_prof value')
   endif

   ! End Routine
   !------------
   return
  end function ! nudging_seat_SRFprofile
  !================================================================

  !================================================================
  ! JS - 04/11/2020 : subroutines for new functions and flexibility
  !================================================================

  !-----------------------
  ! open a new netcdf file
  !-----------------------
  subroutine open_netcdf (ncid, incre_step, Nudge_File_Template)
  use cam_abortutils, only : endrun
  use perf_mod
  use netcdf
  use filenames ,only      : interpret_filename_spec

  integer, intent(out)          :: ncid
  integer, intent(in)           :: incre_step
  character(len=cl), intent(in) :: Nudge_File_Template

  ! local variable  
  integer :: YMD3, YMD4, Nudge_Next1_Sec, Nudge_Next1_Year, &
             Nudge_Next1_Month, Nudge_Next1_Day, istat
  character(len=cl)       :: nudge_file

  YMD3 = (Nudge_Next_Year*10000) + &
         (Nudge_NEXT_Month*100) + Nudge_Next_Day

  call timemgr_time_inc(YMD3,Nudge_Next_Sec,  &
                        YMD4,Nudge_Next1_Sec, &
                        incre_step,0,0)

  Nudge_Next1_Year = YMD4/10000
  YMD4 = YMD4-(Nudge_Next1_Year*10000)
  Nudge_Next1_Month = YMD4/100
  Nudge_Next1_Day   = YMD4-(Nudge_Next1_Month*100)

  nudge_file = interpret_filename_spec (  &
                    Nudge_File_Template,  &
               yr_spec=Nudge_Next1_Year,  &
              mon_spec=Nudge_Next1_Month, &
              day_spec=Nudge_Next1_Day  , &
              sec_spec=Nudge_Next1_Sec    )

  istat=nf90_open(trim(Nudge_Path)//trim(nudge_file),NF90_NOWRITE,ncid)
  if (istat .ne. NF90_NOERR) then
      write(iulog,*) 'NF90_OPEN: failed for file',trim(nudge_file)
      write(iulog,*) nf90_strerror(istat)
      call endrun ('UPDATE_ANALYSES_OPEN_NETCDF')
  endif
  write(iulog,*) 'NUDGING: Reading new analyses:',trim(Nudge_Path)//trim(nudge_file)

  return
  end subroutine !open_netcdf

  !--------------------------------------------
  ! Get and scatter nudging data for SE 2d var 
  !--------------------------------------------
  subroutine read_and_scatter_se_2d (ncid, vname, strt2, cnt2, out_x)
  use ppgrid, only                 : pver,pcols,begchunk,endchunk
  use cam_abortutils, only         : endrun
  use perf_mod
  use netcdf

  integer, intent(in)             :: ncid
  integer, intent(in)             :: strt2(2), cnt2(2)
  character (len = *), intent(in) :: vname
  real(r8), intent(out)           :: out_x(pcols,begchunk:endchunk)

  ! local variables
  real(r8)                        :: Xanal(Nudge_ncol)
  real(r8)                        :: tinfo
  integer                         :: istat, varid, varid1

  if (masterproc) then
    call t_startf ('read_nudging_data')
    istat = nf90_inq_varid(ncid,vname,varid)
    if (istat .ne. NF90_NOERR) then
        write(iulog,*) nf90_strerror(istat)
        call endrun ('ANALYSES_SE_INQ_VARID')
    end if
    istat = nf90_get_var(ncid,varid,Xanal,strt2,cnt2)
    if (istat .ne. NF90_NOERR) then
        write(iulog,*) nf90_strerror(istat)
        call endrun ('ANALYSES_SE_GET_VAR')
    end if
    call t_stopf ('read_nudging_data')

    ! check whether the time slice is read in correctly
    istat = nf90_inq_varid(ncid,'time',varid1)
    if (istat .ne. NF90_NOERR) then
         write(iulog,*) nf90_strerror(istat)
         call endrun ('ANALYSES_SE_INQ_VARID')
     end if
     istat = nf90_get_var(ncid,varid1,tinfo,start=(/strt2(2)/))
     if (istat .ne. NF90_NOERR) then
         write(iulog,*) nf90_strerror(istat)
         call endrun ('ANALYSES_SE_GET_VAR')
     end if
     write(iulog,*) 'NUDGING: Current time slice is: ', tinfo, ', strt2(2) = ', strt2(2), ', reading variable: ', vname
  end if
  call t_startf ('distribute_data')
  call scatter_field_to_chunk(1,1,1,Nudge_ncol,Xanal,out_x)
  call t_stopf ('distribute_data')

  return
  end subroutine !read_and_scatter_se_2d

  !------------------------------------
  ! Get and scatter nudging data for SE
  !------------------------------------
  subroutine read_and_scatter_se (ncid, vname, strt3, cnt3, out_x)
  use ppgrid, only                 : pver,pcols,begchunk,endchunk
  use cam_abortutils, only         : endrun
  use perf_mod
  use netcdf

  integer, intent(in)             :: ncid
  integer, intent(in)             :: strt3(3), cnt3(3)
  character (len = *), intent(in) :: vname
  real(r8), intent(out)           :: out_x(pcols,pver,begchunk:endchunk)

  ! local variables
  real(r8)                        :: Xanal(Nudge_ncol,Nudge_nlev)
  real(r8)                        :: tinfo
  integer                         :: istat, varid, varid1

  if (masterproc) then
    call t_startf ('read_nudging_data')
    istat = nf90_inq_varid(ncid,vname,varid)
    if (istat .ne. NF90_NOERR) then
        write(iulog,*) nf90_strerror(istat)
        call endrun ('ANALYSES_SE_INQ_VARID')
    end if
    istat = nf90_get_var(ncid,varid,Xanal,strt3,cnt3)
    if (istat .ne. NF90_NOERR) then
        write(iulog,*) nf90_strerror(istat)
        call endrun ('ANALYSES_SE_GET_VAR')
    end if
    call t_stopf ('read_nudging_data')

    ! check whether the time slice is read in correctly
    istat = nf90_inq_varid(ncid,'time',varid1)
    if (istat .ne. NF90_NOERR) then         
         write(iulog,*) nf90_strerror(istat)                
         call endrun ('ANALYSES_SE_INQ_VARID')              
     end if         
     istat = nf90_get_var(ncid,varid1,tinfo,start=(/strt3(3)/))             
     if (istat .ne. NF90_NOERR) then                
         write(iulog,*) nf90_strerror(istat)                
         call endrun ('ANALYSES_SE_GET_VAR')                
     end if         
     write(iulog,*) 'NUDGING: Current time slice is: ', tinfo, ', strt3(3) = ', strt3(3), ', reading variable: ', vname
  end if            
  call t_startf ('distribute_data')         
  call scatter_field_to_chunk(1,Nudge_nlev,1,Nudge_ncol,Xanal,out_x)                
  call t_stopf ('distribute_data')          
               
  return
  end subroutine !read_and_scatter_se              
                
  !------------------------------------        
  ! Get and scatter nudging data for FV         
  !------------------------------------         
  subroutine read_and_scatter_fv (ncid, vname, strt4, cnt4, out_x)           
  use ppgrid, only                 : pver,pcols,begchunk,endchunk              
  use cam_abortutils, only         : endrun            
  use perf_mod         
  use netcdf           
               
  integer, intent(in)             :: ncid 
  integer, intent(in)             :: strt4(4), cnt4(4)         
  character (len = *), intent(in) :: vname             
  real(r8), intent(out)           :: out_x(pcols,pver,begchunk:endchunk)               
               
  ! local variables            
  real(r8)                        :: Xanal(Nudge_nlon,Nudge_nlat,Nudge_nlev)           
  real(r8)                        :: Xtrans(Nudge_nlon,Nudge_nlev,Nudge_nlat)          
  real(r8)                        :: tinfo             
  integer                         :: istat, varid, varid1, &           
                                     ilat, ilon, ilev          
               
  if (masterproc) then              
     call t_startf ('read_nudging_data')            
     istat = nf90_inq_varid(ncid,vname,varid)               
     if (istat .ne. NF90_NOERR) then                
         write(iulog,*) nf90_strerror(istat)                
         call endrun ('ANALYSES_FV_INQ_VARID')              
     end if         
     istat = nf90_get_var(ncid,varid,Xanal,strt4,cnt4)              
     if (istat .ne. NF90_NOERR) then                
         write(iulog,*) nf90_strerror(istat)                
         call endrun ('ANALYSES_FV_GET_VAR')                
     end if         
     call t_stopf ('read_nudging_data')             
            
     ! check whether the time slice is read in correctly            
     istat = nf90_inq_varid(ncid,'time',varid1)             
     if (istat .ne. NF90_NOERR) then                
         write(iulog,*) nf90_strerror(istat)                
         call endrun ('ANALYSES_FV_INQ_VARID')              
     end if         
     istat = nf90_get_var(ncid,varid1,tinfo,start=(/strt4(4)/))             
     if (istat .ne. NF90_NOERR) then                
         write(iulog,*) nf90_strerror(istat)
         call endrun ('ANALYSES_FV_GET_VAR')
    end if
    write(iulog,*) 'NUDGING: Current time slice is: ', tinfo, ', strt4(4) = ', strt4(4), ', reading variable: ', vname
    do ilat = 1, Nudge_nlat
       do ilev = 1, pver
          do ilon = 1, Nudge_nlon
             Xtrans(ilon,ilev,ilat) = Xanal(ilon,ilat,ilev)
          end do
       end do
    end do
  end if
  call t_startf ('distribute_data')
  call scatter_field_to_chunk(1,Nudge_nlev,1,Nudge_nlon,Xtrans,out_x)
  call t_stopf ('distribute_data')

  return 
  end subroutine !read_and_scatter_fv

  !---------------------
  ! Linear interpolation
  !---------------------
  subroutine linear_interpolation (input, output)
  use ppgrid,      only : pver,pcols,begchunk,endchunk
  real(r8), intent(in)  :: input(pcols,pver,begchunk:endchunk,2)
  real(r8), intent(out) :: output(pcols,pver,begchunk:endchunk)

  ! Local variable
  integer               :: Year, Month, Day, Sec
  real(r8)              :: factor
  call get_curr_date(Year,Month,Day,Sec)

  if ( Nudge_CurrentStep ) then
     factor = (Sec - (Sec/Nudge_Step)*Nudge_Step) * &
              1._r8 / (Nudge_Step * 1._r8) 
  else
     factor = (Sec - (Sec/Nudge_Step)*Nudge_Step + Model_Step) * &
              1._r8 / (Nudge_Step * 1._r8)
  end if

  if ( factor .eq. 1._r8 ) then
     output = input(:,:,:,2)
  else if ( factor .eq. 0._r8 ) then
     output = input(:,:,:,1)
  else
     output = ( input(:,:,:,2) - input(:,:,:,1) ) * &
              factor + input(:,:,:,1)
  end if
 
  return
  end subroutine !linear_interpolation

  !----------------------------------
  ! Linear interpolation 2d variable
  !----------------------------------
  subroutine linear_interpolation_2d (input, output)
  use ppgrid,      only : pver,pcols,begchunk,endchunk
  real(r8), intent(in)  :: input(pcols,begchunk:endchunk,2)
  real(r8), intent(out) :: output(pcols,begchunk:endchunk)

  ! Local variable
  integer               :: Year, Month, Day, Sec
  real(r8)              :: factor
  call get_curr_date(Year,Month,Day,Sec)

  if ( Nudge_CurrentStep ) then
     factor = (Sec - (Sec/Nudge_Step)*Nudge_Step) * &
              1._r8 / (Nudge_Step * 1._r8)
  else
     factor = (Sec - (Sec/Nudge_Step)*Nudge_Step + Model_Step) * &
              1._r8 / (Nudge_Step * 1._r8)
  end if

  if ( factor .eq. 1._r8 ) then
     output = input(:,:,2)
  else if ( factor .eq. 0._r8 ) then
     output = input(:,:,1)
  else
     output = ( input(:,:,2) - input(:,:,1) ) * &
              factor + input(:,:,1)
  end if

  return
  end subroutine ! linear_interpolation_2d

  !------------------------------------
  ! Update nudging tendency for 3d field SE 
  !------------------------------------
  subroutine update_nudging_tend (ncol, dtime, psmod, psobs, psfac,     & !In
                                  umod, uobs, ufac, vmod, vobs, vfac,   & !In 
                                  tmod, tobs, tfac, qmod, qobs, qfac,   & !In 
                                  ubobs, vbobs, tbobs, tdbobs, qbobs,   & !In 
                                  sfac, ndg_srf_on, ndg_srf_state_on,   & !In
                                  ndg_srf_q, phis_mod, phis_obs, pblh,  & !In
                                  use_ps_adj, use_q_adj, use_pdep_nudge,& !In 
                                  use_upp_lrelx, no_pbl_uv, no_pbl_t,   & !In 
                                  no_pbl_q, ndg_ps_flg, ndg_ps_opt,     & !In
                                  ndg_u_flg, ndg_v_flg, ndg_uv_opt,     & !In
                                  ndg_t_flg, ndg_t_opt,                 & !In 
                                  ndg_q_flg, ndg_q_opt,                 & !In 
                                  zm_obs, zm_mod, psdt,                 & !Out
                                  udt, vdt, tdt, qdt)                     !Out 

  use hycoef,        only  : hycoef_init, hyam, hybm, hyai, hybi, ps0
  use ppgrid,        only  : pver,pverp,pcols
  use shr_vmath_mod, only  : shr_vmath_log
  use physconst,     only  : rga, cpair, gravit, rair, latvap, rh2o, zvir, cappa
  use constituents,  only  : cnst_get_ind, pcnst
  use cam_abortutils,only  : endrun
  use wv_saturation, only  : qsat, qsat_water, svp_ice
  use geopotential,  only  : geopotential_t

  integer, intent(in)  :: ncol    ! number of columns
 
  logical, intent(in)  :: use_pdep_nudge
  logical, intent(in)  :: use_upp_lrelx
  logical, intent(in)  :: use_ps_adj
  logical, intent(in)  :: use_q_adj

  integer, intent(in)  :: no_pbl_uv
  integer, intent(in)  :: no_pbl_t
  integer, intent(in)  :: no_pbl_q

  real(r8), intent(in) :: dtime
  real(r8), intent(in) :: pblh(pcols)

  !!variables for land-surface nudging 
  logical, intent(in)  :: ndg_srf_on
  logical, intent(in)  :: ndg_srf_state_on
  logical, intent(in)  :: ndg_srf_q

  real(r8), intent(in) :: sfac(pcols)
  real(r8), intent(in) :: ubobs(pcols) ! 10m wind 
  real(r8), intent(in) :: vbobs(pcols) ! 10m wind 
  real(r8), intent(in) :: tbobs(pcols) ! 2m temperature 
  real(r8), intent(in) :: tdbobs(pcols)! 2m dewpoint tempeature 
  real(r8), intent(in) :: qbobs(pcols) ! 2m humidity

  !!variables for 3-D atmosphere nudging 
  integer, intent(in)  :: ndg_ps_opt
  integer, intent(in)  :: ndg_uv_opt
  integer, intent(in)  :: ndg_t_opt
  integer, intent(in)  :: ndg_q_opt

  integer, intent(in)  :: ndg_ps_flg
  integer, intent(in)  :: ndg_u_flg
  integer, intent(in)  :: ndg_v_flg
  integer, intent(in)  :: ndg_t_flg
  integer, intent(in)  :: ndg_q_flg

  real(r8), intent(in) :: phis_mod(pcols)
  real(r8), intent(in) :: phis_obs(pcols)

  real(r8), intent(inout) :: psmod(pcols)
  real(r8), intent(inout) :: umod(pcols,pver)
  real(r8), intent(inout) :: vmod(pcols,pver)
  real(r8), intent(inout) :: tmod(pcols,pver)
  real(r8), intent(inout) :: qmod(pcols,pver)

  real(r8), intent(inout) :: psobs(pcols)
  real(r8), intent(inout) :: uobs(pcols,pver)
  real(r8), intent(inout) :: vobs(pcols,pver)
  real(r8), intent(inout) :: tobs(pcols,pver)
  real(r8), intent(inout) :: qobs(pcols,pver)

  real(r8), intent(inout) :: psfac(pcols)
  real(r8), intent(inout) :: ufac(pcols,pver)
  real(r8), intent(inout) :: vfac(pcols,pver)
  real(r8), intent(inout) :: tfac(pcols,pver)
  real(r8), intent(inout) :: qfac(pcols,pver)

  real(r8), intent(inout) :: psdt(pcols)
  real(r8), intent(inout) :: udt(pcols,pver)
  real(r8), intent(inout) :: vdt(pcols,pver)
  real(r8), intent(inout) :: tdt(pcols,pver)
  real(r8), intent(inout) :: qdt(pcols,pver)

  real(r8), intent(out)   :: zm_mod(pcols,pver)
  real(r8), intent(out)   :: zm_obs(pcols,pver)

  ! local variables 
  real(r8) :: pint_obs(pcols,pverp)
  real(r8) :: pmid_obs(pcols,pver)
  real(r8) :: pdel_obs(pcols,pver)
  real(r8) :: rpdel_obs(pcols,pver)
  real(r8) :: lnpint_obs(pcols,pverp)
  real(r8) :: lnpmid_obs(pcols,pver)
  real(r8) :: exner_obs(pcols,pver)
  real(r8) :: zi_obs(pcols,pverp)

  real(r8) :: pint_mod(pcols,pverp)
  real(r8) :: pmid_mod(pcols,pver)
  real(r8) :: pdel_mod(pcols,pver)
  real(r8) :: rpdel_mod(pcols,pver)
  real(r8) :: lnpint_mod(pcols,pverp)
  real(r8) :: lnpmid_mod(pcols,pver)
  real(r8) :: exner_mod(pcols,pver)
  real(r8) :: zi_mod(pcols,pverp)

  real(r8) :: dqsdT_mod(pcols,pver)
  real(r8) :: dqsdT_obs(pcols,pver)
  real(r8) :: qsmod(pcols,pver)
  real(r8) :: qsobs(pcols,pver)
  real(r8) :: esmod(pcols,pver)
  real(r8) :: esobs(pcols,pver)
  real(r8) :: rhmod(pcols,pver)
  real(r8) :: rhobs(pcols,pver)
  real(r8) :: tvmod(pcols,pver)
  real(r8) :: tvobs(pcols,pver)

  real(r8) :: rairv(pcols,pver)
  real(r8) :: zvirv(pcols,pver)
  real(r8) :: qref(pcols,pver)
  real(r8) :: tref(pcols,pver)
  real(r8) :: uref(pcols,pver)
  real(r8) :: vref(pcols,pver)
  real(r8) :: psref(pcols)

  integer  :: i, k, m
  logical  :: l_adj_super_saturation 

  ! initialize the flag for supersaturation adjustment 
  if(ndg_t_opt > 0 .or. ndg_q_opt > 0 ) then
    l_adj_super_saturation = .true.
  else
    l_adj_super_saturation = .false.
  end if

  ! initialize all
  do i = 1, ncol
    psdt(i)  = 0._r8
    psref(i) = psobs(i)
    do k = 1, pver
      udt(i,k)  = 0._r8
      vdt(i,k)  = 0._r8
      tdt(i,k)  = 0._r8
      qdt(i,k)  = 0._r8
      rairv(i,k) = rair
      zvirv(i,k) = zvir
    end do
  end do

  !compute pressure and ln(pres) at layer interfaces
  do k = 1, pver
     do i = 1, ncol 
       pint_obs(i,k) = hyai(k)*ps0 + hybi(k)* psobs(i)
       pmid_obs(i,k) = hyam(k)*ps0 + hybm(k)* psobs(i)
       pint_mod(i,k) = hyai(k)*ps0 + hybi(k)* psmod(i)
       pmid_mod(i,k) = hyam(k)*ps0 + hybm(k)* psmod(i)
     end do 
     !logrithm of pressure 
     call shr_vmath_log(pint_obs(1:ncol,k),lnpint_obs(1:ncol,k),ncol)
     call shr_vmath_log(pmid_obs(1:ncol,k),lnpmid_obs(1:ncol,k),ncol)
     call shr_vmath_log(pint_mod(1:ncol,k),lnpint_mod(1:ncol,k),ncol)
     call shr_vmath_log(pmid_mod(1:ncol,k),lnpmid_mod(1:ncol,k),ncol)
  end do 

  !top interface level pver+1 
  do i=1,ncol
     pint_obs(i,pverp)=hyai(pverp)*ps0+hybi(pverp)*psobs(i)
     pint_mod(i,pverp)=hyai(pverp)*ps0+hybi(pverp)*psmod(i)
  end do
  !logrithm of pressure 
  call shr_vmath_log(pint_obs(1:ncol,pverp),lnpint_obs(1:ncol,pverp),ncol)
  call shr_vmath_log(pint_mod(1:ncol,pverp),lnpint_mod(1:ncol,pverp),ncol)

  !derive the layer thickness and exner
  do k = 1, pver
     do i = 1, ncol
        pdel_obs(i,k)  = pint_obs(i,k+1) - pint_obs(i,k)
        pdel_mod(i,k)  = pint_mod(i,k+1) - pint_mod(i,k)
        rpdel_obs(i,k) = 1._r8/pdel_obs(i,k)
        rpdel_mod(i,k) = 1._r8/pdel_mod(i,k)
        exner_mod(i,k) = (pint_mod(i,pverp)/pmid_mod(i,k))**cappa
        exner_obs(i,k) = (pint_obs(i,pverp)/pmid_obs(i,k))**cappa 
     end do
  end do

  !derive the geoptential height 
  call geopotential_t(lnpint_obs, lnpmid_obs, pint_obs, pmid_obs, pdel_obs, rpdel_obs, &
                      tobs, qobs, rairv,gravit,zvirv,zi_obs,zm_obs,ncol)

  call geopotential_t(lnpint_mod, lnpmid_mod, pint_mod, pmid_mod, pdel_mod, rpdel_mod, &
                      tmod, qmod, rairv,gravit,zvirv,zi_mod,zm_mod,ncol)

  !Nudge surface pressure 
  if ( ndg_ps_flg > 0 ) then
    select case (ndg_ps_opt)
       case (0)
        !mimic the ps nudging in GFDL model 
         do k = 1, pver
           do i = 1, ncol
             uref(i,k) = umod(i,k)
             vref(i,k) = vmod(i,k) 
             tref(i,k) = tmod(i,k)*(1.0_r8 + zvir * qmod(i,k)) 
             qref(i,k) = qmod(i,k)
           end do
         end do
         call ps_nudging(ncol, dtime, use_ps_adj, zi_obs, zm_obs, &
                         pint_obs, pint_mod, pmid_obs, pmid_mod, &
                         lnpint_obs, lnpint_mod, lnpmid_obs, lnpmid_mod, & 
                         pdel_obs, pdel_mod, rpdel_obs, rpdel_mod, &
                         uref, vref, tref, qref, tobs, psref, &
                         phis_obs, phis_mod, psobs, psmod, psfac, psdt )
         do k = 1, pver
           do i = 1, ncol
             udt(i,k)  = udt(i,k) + (uref(i,k) - umod(i,k)) / dtime
             vdt(i,k)  = vdt(i,k) + (vref(i,k) - vmod(i,k)) / dtime
             tdt(i,k)  = tdt(i,k) + (tref(i,k)/(1.0_r8 + zvir * qref(i,k)) - tmod(i,k)) / dtime
             qdt(i,k)  = qdt(i,k) + (qref(i,k) - qmod(i,k)) / dtime
           enddo
         enddo
       case default
         call endrun('nudging_tend error: invalid option for PS nudging')
    end select
  
  end if

  !Apply scaling on nudging strength in upper and bottom model layers 
  call update_nudge_prof (ncol, zm_mod, zi_mod, pmid_mod, pblh, & !In 
                          use_pdep_nudge, use_upp_lrelx,        & !In 
                          no_pbl_uv, no_pbl_t, no_pbl_q,        & !In 
                          ufac, vfac, tfac, qfac)                 !Out 

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Start to calculate nudging tendency 
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !calculate the virtual temperature, saturation mixing ratio and relative humidity 
  !these information are needed for option 1 and 2 of T nudging 
  !and for option 1 of Q nudging 
  if ( ndg_t_opt > 0 .or. ndg_q_opt > 0 .or. ndg_srf_on ) then
   ! calculate the virtual temperature, saturation mixing ratio and relative
   ! humidity 
     do k = 1, pver
      do i = 1, ncol
        call qsat(tobs(i,k), pmid_obs(i,k), esobs(i,k), qsobs(i,k), dqsdt=dqsdT_obs(i,k))
        call qsat(tmod(i,k), pmid_mod(i,k), esmod(i,k), qsmod(i,k), dqsdt=dqsdT_mod(i,k))
        rhobs(i,k) = qobs(i,k) / qsobs(i,k)
        rhmod(i,k) = qmod(i,k) / qsmod(i,k)
        tvobs(i,k) = tobs(i,k) * (1.0_r8 + zvir * qobs(i,k))
        tvmod(i,k) = tmod(i,k) * (1.0_r8 + zvir * qmod(i,k))
      end do
    end do
    ! limit humidity in the observational data 
    do k = 1, pver
      do i = 1, ncol
        qref(i,k)  = max(1.E-8_r8, qobs(i,k))
        ! adjust supersaturation
        if(l_adj_super_saturation .and. (qobs(i,k) > qsobs(i,k))) then
           qref(i,k)  = qsobs(i,k)
           rhobs(i,k) = 1._r8
        end if
      end do
    end do
  end if

  ! zonal wind component
  if ( ndg_u_flg > 0 ) then 
    select case (ndg_uv_opt)
       case (0)
         ! direct nudging 
         udt(:,:) = udt(:,:) + (uobs(:,:) - umod(:,:))*ufac(:,:)  
       case default
       call endrun('nudging_tend error: invalid option for U nudging')
    end select
  end if 

  ! meridional wind component
  if ( ndg_v_flg > 0 ) then
    select case (ndg_uv_opt)
       case (0)
         ! direct nudging 
         vdt(:,:) = vdt(:,:) + (vobs(:,:) - vmod(:,:))*vfac(:,:)
       case default
       call endrun('nudging_tend error: invalid option for V nudging')
    end select
  end if

  ! temperature 
  if ( ndg_t_flg > 0 ) then
    select case (ndg_t_opt)
       case (0)
         ! direct nudging 
         tdt(:,:) = tdt(:,:) + (tobs(:,:) - tmod(:,:))*tfac(:,:) 
       case (1)
         ! nudge virtual temperature (Tv) 
         ! Tangent linear of Tv with respect to T and q 
         ! dT/dTv = d(Tv/(1 + zvir * q))/dTv = 1/(1 + zvir * q)
         ! dT     = dTv / (1 + zvir * q) 
         do k = 1, pver
          do i = 1, ncol
            tdt(i,k) = tdt(i,k) + (tvobs(i,k) - tvmod(i,k))*tfac(i,k) &
                                  /(1.0_r8 + zvir * qmod(i,k)) 
          end do            
         end do 
       case (2)
         ! Reference: Yang et. al. (2021, JMR), Fan and Tilley (2005, MWR)
         ! nudge relative humidity (RH) to adjust T 
         ! Tangent linear of RH with respect to T and q 
         ! dRH/dT = d(q/qs)/dT = - q/qs^2 * dqs/dT 
         ! dT     = dRH / (-RH/qs * dqs/dT )
         do k = 1, pver
          do i = 1, ncol
            if ( (qmod(i,k) > 0._r8) .and. (dqsdT_mod(i,k) > 0._r8) ) then
              tdt(i,k) = tdt(i,k) - (rhobs(i,k)-rhmod(i,k))*qfac(i,k) &
                                   *qsmod(i,k)*qsmod(i,k)/(max(qmod(i,k),1.E-8_r8)*dqsdT_mod(i,k))
            else 
              tdt(i,k) = tdt(i,k) + (tvobs(i,k) - tvmod(i,k))*tfac(i,k) &
                                   /(1.0_r8 + zvir * qmod(i,k))
            end if     
          end do 
         end do
       case default
         call endrun('nudging_tend error: invalid option for T nudging')
    end select
  end if 

  ! humidity
  if ( ndg_q_flg > 0 ) then
    select case (ndg_q_opt)
       case (0)
         ! direct nudging 
         qdt(:,:) = qdt(:,:) + (qref(:,:) - qmod(:,:))*qfac(:,:) 
       case (1)
         ! Reference: Yang et. al. (2021, JMR), Fan and Tilley (2005, MWR)
         ! nudge relative humidity (RH), only adjust q
         ! Tangent linear of RH with respect to q 
         ! dRH/dT   = d(q/qs)/dT = - q/qs^2 * dqs/dT = -RH/qs * dqs/dT 
         ! dRH/dq   = d(q/qs)/dq = 1/qs = RH / q 
         ! dq       = qs*dRH
         do k = 1, pver
          do i = 1, ncol
             qdt(i,k) = qdt(i,k) + (rhobs(i,k)-rhmod(i,k))*qsmod(i,k)*qfac(i,k)
          end do
         end do
       case default
         call endrun('nudging_tend error: invalid option for Q nudging')
    end select
  end if   

  ! use the surface pressure to scale down the tendency at surface 
  if (ndg_ps_flg > 0 ) then
    do i = 1,ncol
      udt(i,pver) = udt(i,pver)* psfac(i)
      vdt(i,pver) = vdt(i,pver)* psfac(i)
      tdt(i,pver) = tdt(i,pver)* psfac(i)
      qdt(i,pver) = qdt(i,pver)* psfac(i)
    end do
  end if

  ! land surface nudging 
  if (ndg_srf_on .and. ndg_srf_state_on) then 

    call update_land_nudging_tend(ncol, umod, vmod, tvmod, qmod,    & ! In 
                                  ubobs, vbobs, tbobs, tdbobs,      & ! In 
                                  qbobs, pmid_obs, sfac, ndg_srf_q, & ! In           
                                  udt, vdt, tdt, qdt )                ! Out 

  end if 

  return
  end subroutine !update_nudging_tend

  subroutine mltbc_advance(file_path,varname,nlev,ngcol,model_state,nudging_tend) 
  !===========================================================================
  ! SZ - 06/05/2023: This subroutine attempt to read and calculate the nudging
  !                  tendency using the Machine Learning (ML) model,
  !                  the subroutine works on a specific model state variable
  ! SZ - 06/21/2023: Merge the DeepOnet model developed by Brown University 
  !===========================================================================
   use ppgrid,           only : pver,pverp,pcols,begchunk,endchunk
   use phys_grid,        only : scatter_field_to_chunk 
   use shr_const_mod,    only : SHR_CONST_PI, SHR_CONST_REARTH
   use cam_history  ,    only : outfld

   implicit none
   logical, parameter           :: l_print_diag = .false.
   character(len=*), intent(in) :: file_path !Path to Machine Learning model files 
   character(len=*), intent(in) :: varname   !nudge variable
   integer,intent(in)           :: nlev,ngcol 
   real(r8),intent(in)          :: model_state(ngcol,nlev)
   real(r8),intent(inout)       :: nudging_tend(pcols,nlev,begchunk:endchunk)

   ! Arguments
   !-----------
   type(torch_module)           :: torch_mod
   type(torch_tensor_wrap)      :: input_tensors
   type(torch_tensor)           :: out_tensor

   ! Local variables 
   !----------------
   integer  :: ncols,lchnk
   integer  :: i,j,n,m,k,ii,jj,ierr
   real(r8) :: vmin,vmax
   real(r8) :: sum_x(nlev)
   real(r8) :: wrk_tend(ngcol,nlev)

   !DeepONet convolution 2d model (regional)
   integer,  parameter :: don_conv2d_nx = 6      ! ML model trunk in X 6
   integer,  parameter :: don_conv2d_ny = 6      ! ML model trunk in Y 6
   real(r8), parameter :: don_tend_dtime = 3.0_r8 ! unit: hour

   !DeepONet global model 
   integer,  parameter :: don_glb_nt = 1      ! Dummy time array size
   integer,  parameter :: don_glb_ncol = 2700   ! ML model trunk in X 6
   integer,  parameter :: don_glb_nlev = 64     ! ML model trunk in Y 6

   !Temporary work array 
   real(r8), allocatable :: wrk_state(:,:,:)
   real(r8), allocatable :: enc_out(:,:,:)
   real(r8), allocatable :: don_out(:,:,:)
   real(r8), allocatable :: dcd_out(:,:,:)
   real(r8), allocatable :: wrk_out(:,:)
   real(r8), allocatable :: don_stat(:,:)

   !initialize tendency array 
   wrk_tend(:,:) = 0.0_r8

   !Currently the machine learning model was only trained for U, V only 
   !Return if other variables are passed to this subroutine
   if ((trim(varname) /= 'U') .and. (trim(varname) /= 'V')) then
     write(iulog,*) "Machine Learning Nudging Warning: working variable ", trim(varname)
     write(iulog,*) "Machine Learning Nudging Warning: Convolution 2D model  & 
                     only designed for U, V nudging, return"
     return
   end if

   !start to call deepONet model 
   call t_startf('mltbc_don_prediction')

   if (mltbc_regional_on) then

     allocate (wrk_out(mltbc_patch_nlon*mltbc_patch_nlat,nlev))
     allocate (wrk_state(mltbc_patch_nlon,mltbc_patch_nlat,nlev))

     !collect data for deepONet
     call mltbc_gather_patch(varname,nlev,ngcol,mltbc_patch_nlon,mltbc_patch_nlat, & 
                                model_state,wrk_state) !inout  
     if (masterproc .and. l_print_diag ) then
       write(iulog,*) 'shape of model_state = ', shape(model_state)
       write(iulog,*) 'model_state(min/max) = ', minval(model_state),maxval(model_state)
       write(iulog,*) 'shape of wrk_state   = ', shape(wrk_state)
       write(iulog,*) 'wrk_state(min/max)   = ', minval(wrk_state),maxval(wrk_state)
     end if

     !There are two options for convolution 2D model
     select case (mltbc_predict_option)
       case (0)      
         !option 0: encoder-->deepONet-->decoder approach to 
         !          predict nudging tendency from before nuding state
         allocate (enc_out(don_conv2d_nx,don_conv2d_ny,nlev))
         allocate (don_out(don_conv2d_nx,don_conv2d_ny,nlev))
         allocate (dcd_out(mltbc_patch_nlon,mltbc_patch_nlat,nlev))

         !call encoder 
         call mltbc_don_encoder(mltbc_regional_on,file_path,varname, & 
                                mltbc_patch_nlon,mltbc_patch_nlat,nlev, & 
                                don_conv2d_nx,don_conv2d_ny,nlev, & 
                                wrk_state,enc_out)
 
         !call deeponet               
         call mltbc_don_tendadv(mltbc_regional_on,file_path,varname, &
                                don_conv2d_nx,don_conv2d_ny,nlev, & 
                                enc_out,don_out)

         !call decoder                
         call mltbc_don_decoder(mltbc_regional_on,file_path,varname, & 
                                mltbc_patch_nlon,mltbc_patch_nlat,nlev, &
                                don_conv2d_nx,don_conv2d_ny,nlev, &
                                don_out,dcd_out)

         !prepare for final regrid 
         do j = 1, mltbc_patch_nlat
           do i = 1, mltbc_patch_nlon
             m = (j-1)*mltbc_patch_nlon + i
             wrk_out(m,:) = dcd_out(i,j,:) 
           end do
         end do
         deallocate(enc_out)
         deallocate(don_out)
         deallocate(dcd_out)

       case (1) 
         allocate (don_stat(mltbc_patch_nlon*mltbc_patch_nlat,nlev))
         !option 1: deepONet (before nudging state --> after nudging state) 
         !          nudging tedency = (After - Before) / 3*3600 (3hour)
         call mltbc_don_statadv(file_path,varname,mltbc_patch_nlon,mltbc_patch_nlat, &
                                nlev,wrk_state,don_stat)
         !compute nudging tendency 
         do j = 1, mltbc_patch_nlat
           do i = 1, mltbc_patch_nlon
             m = (j-1)*mltbc_patch_nlon + i 
             !calculate nudging tendency 
             wrk_out(m,:) = (don_stat(m,:) - wrk_state(i,j,:)) / don_tend_dtime / sec_per_hour
           end do
         end do
         deallocate(don_stat)

       case (2) 
         !option 2: ML(state)-->ML(tendency)
         !          predict nudging tendency from before nuding state
         !call machine learning model 
         call mltbc_reg_tendadv(file_path,varname,mltbc_patch_nlon,mltbc_patch_nlat, &
                                nlev,wrk_state,wrk_out)
       case default
         call endrun('Machine Learning Nudging Error: invalid option for regional model option')
     end select

     !Debug Diagnostics 
     if (masterproc .and. l_print_diag) then
       write(iulog,*) 'mltbc_advance: run deepONet successfully'
       write(iulog,*) 'predict variable = ',varname
       write(iulog,*) 'shape of wrk_in  = ',shape(wrk_state)
       write(iulog,*) 'wrk_in(min/max)  = ',minval(wrk_state),maxval(wrk_state)
       write(iulog,*) 'shape of wrk_out = ',shape(wrk_out)
       write(iulog,*) 'wrk_out(min/max) = ',minval(wrk_out),maxval(wrk_out)
     end if

     !convert data back to model grid
     if ( .not. mltbc_patch_biln ) then
       !method 1: nearest to destination grid 
       do j = 1, mltbc_patch_nlat
         do i = 1, mltbc_patch_nlon
           m = (j-1)*mltbc_patch_nlon + i
           n = mltbc_se2latlon_ind(m,5)
           wrk_tend(n,:) = wrk_out(m,:) 
         end do 
       end do
     else 
       !method 2: linear interpolation to destination grid 
       do n = 1, Nudge_ncol
         j = 0 
         sum_x(:) = 0.0_r8
         do i = 1, 4 
           if (mltbc_latlon2se_wgt(n,i) == 0.0_r8) j = j + 1
           m = mltbc_latlon2se_ind(n,i) 
           sum_x(:) = sum_x(:) + mltbc_latlon2se_wgt(n,i) * wrk_out(m,:) 
         end do
         if (j /= 4) then
           wrk_tend(n,:) = sum_x(:)
         end if 
       end do   
     end if

     !release array space   
     deallocate(wrk_out)
     deallocate(wrk_state)

   else ! global model 

     !There are two options for gobal model 
     select case (mltbc_predict_option)
       !option : ML(state)-->ML(tendency)
       !         predict nudging tendency from before nuding state
       !call machine learning model 
       case (0)
         !option 0 w/ ViTO,Unet,M&M: ML(X)-> X_tend approach
         call mltbc_glb_tendadv(file_path,varname,ngcol,nlev,model_state,wrk_tend)
       case (1)
         !option 1 w/ DeepONet: encoder(X)-->Y-->deepONet(Y)-->Z-->decoder(Z)--> approach 
         allocate (enc_out(don_glb_ncol,don_glb_nlev,don_glb_nt))
         allocate (don_out(don_glb_ncol,don_glb_nlev,don_glb_nt))
         allocate (dcd_out(ngcol,nlev,don_glb_nt))
         allocate (wrk_state(ngcol,nlev,don_glb_nt))
         wrk_state(1:ngcol,1:nlev,1) = model_state(1:ngcol,1:nlev)

         !call encoder 
         call mltbc_don_encoder(mltbc_regional_on,file_path,varname, & 
                                ngcol,nlev,don_glb_nt, &
                                don_glb_ncol,don_glb_nlev,don_glb_nt, & 
                                wrk_state,enc_out)
         !call deeponet               
         call mltbc_don_tendadv(mltbc_regional_on,file_path,varname, & 
                                don_glb_ncol,don_glb_nlev,don_glb_nt, & 
                                enc_out,don_out)
         !call decoder                
         call mltbc_don_decoder(mltbc_regional_on,file_path,varname, & 
                                ngcol,nlev,don_glb_nt, &
                                don_glb_ncol,don_glb_nlev,don_glb_nt, &
                                don_out,dcd_out)

         wrk_tend(1:ngcol,1:nlev) = dcd_out(1:ngcol,1:nlev,1) 
         deallocate(wrk_state)
         deallocate(enc_out)
         deallocate(don_out)
         deallocate(dcd_out)
       case default
         call endrun('Machine Learning Nudging Error: invalid option for global model option')
     end select

   end if
   call t_stopf('mltbc_don_prediction')

   !scatter filed to chunk (from wrk_tend to nudging_tend) 
   call scatter_field_to_chunk(1,nlev,1,Nudge_ncol,wrk_tend,nudging_tend)

   return
  end subroutine !mltbc_advance 

  subroutine mltbc_glb_tendadv(file_path,varname,ncol,nz,vari,varo)
  !===========================================================================
  ! SZ - 06/05/2023: This subroutine attempt to call the forecast for
  !                  the Machine Learning (ML) model,
  !===========================================================================
   use cam_abortutils  , only : endrun
   use shr_const_mod,    only : SHR_CONST_PI, SHR_CONST_REARTH

   implicit none
   character(len=*), intent(in) :: file_path !Path to machine learning model files
   character(len=*), intent(in) :: varname   !nudge variable
   integer, intent(in)          :: ncol,nz
   real(r8),intent(in)          :: vari(ncol,nz)
   real(r8),intent(inout)       :: varo(ncol,nz)

   ! Arguments
   !-----------
   type(torch_module)           :: torch_mod
   type(torch_tensor_wrap)      :: input_tensors
   type(torch_tensor)           :: out_tensor

   ! Local values
   !----------------
   logical, parameter           :: l_print_diag = .false.
   integer                      :: i,j,n,m,k,ii,jj
   real(r4)                     :: doninp(ncol,nz,1)
   real(r4), pointer            :: donout(:,:,:)

   if (masterproc) then
     !check if machine learning pt file exist 
     file_predictor = trim(varname)//'_DeepONet.pt'
     inquire(file=trim(file_path)//trim(file_predictor),exist=l_mltbc_predictor)
     if ( .not. l_mltbc_predictor) then
       write(iulog,*) "ERROR: "//trim(file_path)//trim(file_predictor)//" not found!"
       call endrun('Machine Learning Nudging Error: model file not exist')
     end if
   end if
#ifdef SPMD
   call mpibcast(file_predictor,len(file_predictor),mpichar,0,mpicom)
#endif

   call t_startf ('mltbc_input_reorg')
   !prepare input data
   !note: if use do loops, always remember to put lev as the innerest loop 
   doninp(:,:,1) = real(vari(:,:),kind=r4) ! single precision: float 32, 64 
   call t_stopf ('mltbc_input_reorg')

   call t_startf ('mltbc_create_array')
   call input_tensors%create
   call t_stopf ('mltbc_create_array')

   call t_startf ('mltbc_add_array')
   call input_tensors%add_array(doninp)
   call t_stopf ('mltbc_add_array')

   call t_startf ('mltbc_load_mlpt')
   call torch_mod%load(trim(file_path)//trim(file_predictor))
   call t_stopf ('mltbc_load_mlpt')

   call t_startf ('mltbc_forward_model')
   call torch_mod%forward(input_tensors,out_tensor,1)
   call t_stopf ('mltbc_forward_model')

   call t_startf ('mltbc_output_array')
   call out_tensor%to_array(donout)
   call t_stopf ('mltbc_output_array')

   if (masterproc.and.l_print_diag) then
     write(iulog,*) 'shape of doninp = ',shape(doninp)
     write(iulog,*) 'doninp(min/max) = ',minval(doninp),maxval(doninp)
     write(iulog,*) 'shape of donout = ',shape(donout)
     write(iulog,*) 'donout(min/max) = ',minval(donout),maxval(donout)
   end if

   call t_startf ('mltbc_output_reorg')
   !return output data
   varo(:,:) = real(donout(:,:,1),kind=r8)
   call t_stopf ('mltbc_output_reorg')

   return
  end subroutine !mltbc_glb_tendadv

  subroutine mltbc_reg_tendadv(file_path,varname,nx,ny,nz,vari,varo)
  !===========================================================================
  ! SZ - 06/05/2023: This subroutine attempt to call the forecast for 
  !                  the Machine Learning (ML) model,
  !===========================================================================
   use cam_abortutils  , only : endrun
   use shr_const_mod,    only : SHR_CONST_PI, SHR_CONST_REARTH

   implicit none
   character(len=*), intent(in) :: file_path !Path to machine learning model files 
   character(len=*), intent(in) :: varname   !nudge variable
   integer, intent(in)          :: nx,ny,nz
   real(r8),intent(in)          :: vari(nx,ny,nz)
   real(r8),intent(inout)       :: varo(nx*ny,nz)

   ! Arguments
   !-----------
   type(torch_module)           :: torch_mod
   type(torch_tensor_wrap)      :: input_tensors
   type(torch_tensor)           :: out_tensor

   ! Local values
   !----------------
   logical, parameter           :: l_print_diag = .false.
   integer                      :: i,j,n,m,k,ii,jj
   real(r4)                     :: doninp(nx,ny,nz,1)
   real(r4), pointer            :: donout(:,:,:,:)

   if (masterproc) then
     !check if machine learning pt file exist 
     file_predictor = trim(varname)//'_DeepONet.pt'
     inquire(file=trim(file_path)//trim(file_predictor),exist=l_mltbc_predictor)
     if ( .not. l_mltbc_predictor) then
       write(iulog,*) "ERROR: "//trim(file_path)//trim(file_predictor)//" not found!"
       call endrun('Machine Learning Nudging Error: model file not exist')
     end if
   end if
#ifdef SPMD
   call mpibcast(file_predictor,len(file_predictor),mpichar,0,mpicom)
#endif

   call t_startf ('mltbc_input_reorg')
   !prepare input data 
   !do k = 1,nz
   !  do j = 1,ny
   !    do i = 1,nx
   !      doninp(i,j,k,1) = real(vari(i,j,k),kind=r4)
   !    end do
   !  end do
   !end do
   !avoide do loops to obtain computational gain 
   doninp(:,:,:,1) = real(vari(:,:,:),kind=r4)

   call t_stopf ('mltbc_input_reorg')

   call t_startf ('mltbc_create_array')
   call input_tensors%create
   call t_stopf ('mltbc_create_array')

   call t_startf ('mltbc_add_array')
   call input_tensors%add_array(doninp)
   call t_stopf ('mltbc_add_array')

   call t_startf ('mltbc_load_pt')
   call torch_mod%load(trim(file_path)//trim(file_predictor))
   call t_stopf ('mltbc_load_pt')

   call t_startf ('mltbc_forward_model')
   call torch_mod%forward(input_tensors,out_tensor,1)
   call t_stopf ('mltbc_forward_model')

   call t_startf ('mltbc_ouput_array')
   call out_tensor%to_array(donout)
   call t_stopf ('mltbc_output_array')

   if (masterproc.and.l_print_diag) then
     write(iulog,*) 'shape of doninp = ',shape(doninp)
     write(iulog,*) 'doninp(min/max) = ',minval(doninp),maxval(doninp)
     write(iulog,*) 'shape of donout = ',shape(donout)
     write(iulog,*) 'donout(min/max) = ',minval(donout),maxval(donout)
   end if

   call t_startf ('mltbc_output_reorg')
   !return output data
   do j = 1, ny
     do i = 1, nx
       m = (j-1)*nx + i
       varo(m,:) = real(donout(i,j,:,1),kind=r8) 
     end do
   end do
   call t_stopf ('mltbc_output_reorg')

   return
  end subroutine !mltbc_reg_tendadv

  subroutine mltbc_don_tendadv(l_conv2d,file_path,varname,nx,ny,nz,vari,varo)
  !===========================================================================
  ! SZ - 06/05/2023: This subroutine attempt to call the forecast for 
  !                  the Machine Learning (ML) model,
  !===========================================================================
   use cam_abortutils  , only : endrun
   use shr_const_mod,    only : SHR_CONST_PI, SHR_CONST_REARTH

   implicit none
   logical, intent(in)          :: l_conv2d
   character(len=*), intent(in) :: file_path !Path to machine learning model files 
   character(len=*), intent(in) :: varname   !nudge variable
   integer, intent(in)          :: nx,ny,nz
   real(r8),intent(in)          :: vari(nx,ny,nz)
   real(r8),intent(inout)       :: varo(nx,ny,nz)

   ! Arguments
   !-----------
   type(torch_module)           :: torch_mod
   type(torch_tensor_wrap)      :: input_tensors
   type(torch_tensor)           :: out_tensor

   ! Local values
   !----------------
   logical, parameter           :: l_print_diag = .false.
   integer                      :: i,j,n,m,k,ii,jj

   ! convolution 2d (regional)  
   integer, parameter           :: dumy_nt = 248 ! Dummy time array size
   real(r8)                     :: donmin, donmax
   real(r4)                     :: vmax,vmin
   real(r4)                     :: x_trunk(1,dumy_nt)    
   real(r4)                     :: donin4d(nx*ny,dumy_nt,1,nz)
   real(r4), pointer            :: donout4d(:,:,:,:)

   ! global model  
   real(r4)                     :: donin3d(nx,ny,nz)
   real(r4), pointer            :: donout3d(:,:,:)

   if (masterproc) then
     !check if machine learning pt file exist 
     file_predictor = trim(varname)//'_DeepONet.pt'
     inquire(file=trim(file_path)//trim(file_predictor),exist=l_mltbc_predictor)
     if ( .not. l_mltbc_predictor) then
       write(iulog,*) "ERROR: "//trim(file_path)//trim(file_predictor)//" not found!"
       call endrun('Machine Learning Nudging Error: model file not exist')
     end if
   end if
#ifdef SPMD
   call mpibcast(file_predictor,len(file_predictor),mpichar,0,mpicom)
#endif

   if ( l_conv2d ) then 

     !parameters to denormalize DeepONet prediction
     if (trim(varname) == 'U') then
       donmax = 9.087432_r8
       donmin = -9.915943_r8
     else  ! "V"
       donmax = 9.116263_r8
       donmin = -8.079512_r8
     end if

     !prepare input data with normalization using min/max 
     vmin = real(minval(vari),kind=r4)
     vmax = real(maxval(vari),kind=r4)
     do j = 1,ny
       do i = 1,nx
         m = (j-1)*nx + i
         donin4d(m,1,1,:) = 2.0_r4*(real(vari(i,j,:),kind=r4)-vmin)/(vmax-vmin)-1.0_r4
       end do
     end do

     !convolution 2D model options 
     !0: before nudging state->nudging tendency 
     !dummy input time
     do n = 1,dumy_nt
       x_trunk(1,n) = 2.0_r4*(real(n,kind=r4)-1.0_r4)/(real(dumy_nt,kind=r4)-1.0_r4)-1.0_r4
       if ( n > 1 ) then 
         donin4d(:,n,1,:) = donin4d(:,1,1,:)
       end if 
     end do

     call input_tensors%create
     call input_tensors%add_array(donin4d)
     call input_tensors%add_array(x_trunk)
     call torch_mod%load(trim(file_path)//trim(file_predictor))
     call torch_mod%forward(input_tensors,out_tensor,1)
     call out_tensor%to_array(donout4d)

     if (masterproc.and.l_print_diag) then
       write(iulog,*) 'shape of donin  = ',shape(donin4d)
       write(iulog,*) 'donin(min/max)  = ',minval(donin4d),maxval(donin4d)
       write(iulog,*) 'shape of donout = ',shape(donout4d)
       write(iulog,*) 'donout(min/max) = ',minval(donout4d),maxval(donout4d)
     end if

     !denormalize the deeponet prediction and output data
     do j = 1, ny
       do i = 1, nx
         m = (j-1)*nx + i
         varo(i,j,:) = 0.5_r8*(real(donout4d(i,j,1,:),kind=r8)+1.0_r8)*(donmax-donmin)+donmin
       end do
     end do
   else 
     !prepare input data 
     donin3d(1:nx,1:ny,1:nz) = real(vari(1:nx,1:ny,1:nz),kind=r4)
     ! call deepOnet 
     call input_tensors%create
     call input_tensors%add_array(donin3d)
     call torch_mod%load(trim(file_path)//trim(file_predictor))
     call torch_mod%forward(input_tensors,out_tensor,1)
     call out_tensor%to_array(donout3d)

     if (masterproc.and.l_print_diag) then
       write(iulog,*) 'shape of donin  = ',shape(donin3d)
       write(iulog,*) 'donin(min/max)  = ',minval(donin3d),maxval(donin3d)
       write(iulog,*) 'shape of donout = ',shape(donout3d)
       write(iulog,*) 'donout(min/max) = ',minval(donout3d),maxval(donout3d)
     end if
     !output data 
     varo(1:nx,1:ny,1:nz) = real(donout3d(1:nx,1:ny,1:nz),kind=r8)
   end if 

   return
  end subroutine !mltbc_don_tendadv

  subroutine mltbc_don_statadv(file_path,varname,nx,ny,nz,vari,varo)
  !===========================================================================
  ! SZ - 06/05/2023: This subroutine attempt to call the forecast for 
  !                  the Machine Learning (ML) model,
  !===========================================================================
   use cam_abortutils  , only : endrun
   use shr_const_mod,    only : SHR_CONST_PI, SHR_CONST_REARTH

   implicit none
   character(len=*), intent(in) :: file_path !Path to machine learning model files 
   character(len=*), intent(in) :: varname   !nudge variable
   integer, intent(in)          :: nx,ny,nz
   real(r8),intent(in)          :: vari(nx,ny,nz)
   real(r8),intent(inout)       :: varo(nx*ny,nz)

   ! Arguments
   !-----------
   type(torch_module)           :: torch_mod
   type(torch_tensor_wrap)      :: input_tensors
   type(torch_tensor)           :: out_tensor

   ! Local values
   !----------------
   logical, parameter           :: l_print_diag = .false.
   integer                      :: i,j,n,m,k,ii,jj
   real(r4)                     :: x_trunk(2,nx*ny)
   real(r4)                     :: doninp(ny,nx,1,nz)
   real(r4), pointer            :: donout(:,:)
   real(r4)                     :: vmini,vmaxi
   real(r8)                     :: vmino,vmaxo

   if (masterproc) then
     !check if machine learning pt file exist 
     file_predictor = trim(varname)//'_DeepONet.pt'
     inquire(file=trim(file_path)//trim(file_predictor),exist=l_mltbc_predictor)
     if ( .not. l_mltbc_predictor) then
       write(iulog,*) "ERROR: "//trim(file_path)//trim(file_predictor)//" not found!"
       call endrun('Machine Learning Nudging Error: model file not exist')
     end if
   end if
#ifdef SPMD
   call mpibcast(file_predictor,len(file_predictor),mpichar,0,mpicom)
#endif

   !convolution 2D model options
   !before nudging state->after nudging state
   !input trunk data 
   do j = 1, ny
     do i = 1, nx
       n = (j-1)*nx + i
       x_trunk(1,n) = 2.0_r4*(real(i,kind=r4)-1.0_r4)/(real(nx,kind=r4)-1.0_r4)-1.0_r4
       x_trunk(2,n) = 2.0_r4*(real(j,kind=r4)-1.0_r4)/(real(ny,kind=r4)-1.0_r4)-1.0_r4
     end do
   end do

   !prepare input data 
   vmini = real(minval(vari(:,:,:)),kind=r4)
   vmaxi = real(maxval(vari(:,:,:)),kind=r4)
   do j = 1,ny
     do i = 1,nx
       !normalize input data
       doninp(i,j,1,:) = 2.0_r4*(real(vari(i,j,:),kind=r4)-vmini)/(vmaxi-vmini)-1.0_r4
     end do
   end do
   call input_tensors%create
   call input_tensors%add_array(doninp)
   call input_tensors%add_array(x_trunk)
   call torch_mod%load(trim(file_path)//trim(file_predictor))
   call torch_mod%forward(input_tensors,out_tensor,1)
   call out_tensor%to_array(donout)

   if (masterproc.and.l_print_diag) then
     write(iulog,*) 'shape of doninp = ',shape(doninp)
     write(iulog,*) 'doninp(min/max) = ',minval(doninp),maxval(doninp)
     write(iulog,*) 'shape of donout = ',shape(donout)
     write(iulog,*) 'donout(min/max) = ',minval(donout),maxval(donout)
   end if

  !output data (denormalize)
   vmino = minval(vari(:,:,:))
   vmaxo = maxval(vari(:,:,:))
   do j = 1, ny
     do i = 1, nx
       m = (j-1)* nx + i
       varo(m,:) = 0.5_r8*(real(donout(m,:),kind=r8)+1.0_r8)*(vmaxo-vmino)+vmino
     end do
   end do

   return
  end subroutine !mltbc_don_statadv

  subroutine mltbc_don_encoder(l_conv2d,file_path,varname,nx,ny,nz,nx1,ny1,nz1,vari,varo)
  !===========================================================================
  ! SZ - 06/05/2023: This subroutine attempt to call auto encoder for 
  !                  the Machine Learning (ML) model,
  !===========================================================================
   use cam_abortutils  , only : endrun
   use shr_const_mod,    only : SHR_CONST_PI, SHR_CONST_REARTH

   implicit none
   logical, intent(in)          :: l_conv2d
   character(len=*), intent(in) :: file_path !Path to machine learning model files 
   character(len=*), intent(in) :: varname   !nudge variable
   integer, intent(in)          :: nx,ny,nz,nx1,ny1,nz1
   real(r8),intent(in)          :: vari(nx,ny,nz) 
   real(r8),intent(inout)       :: varo(nx1,ny1,nz1)

   ! Arguments
   !-----------
   type(torch_module)           :: torch_mod
   type(torch_tensor_wrap)      :: input_tensors
   type(torch_tensor)           :: out_tensor

   ! Local values
   !----------------
   logical, parameter :: l_print_diag = .false.
   integer            :: i,j,n,m,k
   real(r4)           :: vmax,vmin 
   real(r4)           :: encin3d(nx,ny,nz)   ! input array to encoder 
   real(r4)           :: encin4d(nx,ny,1,nz) ! input array to encoder 
   real(r4), pointer  :: encout2d(:,:)       ! output array from encoder 
   real(r4), pointer  :: encout3d(:,:,:)     ! output array from encoder  

   if (masterproc) then
     !check if machine learning pt file exist 
     file_encoder = trim(varname)//'_Encoder.pt'
     inquire(file=trim(file_path)//trim(file_encoder),exist=l_mltbc_encoder)
     if ( .not. l_mltbc_encoder) then
       write(iulog,*) "ERROR: "//trim(file_path)//trim(file_encoder)//" not found!"
       call endrun('Machine Learning Nudging Error: model file not exist')
     end if
   end if
#ifdef SPMD
   call mpibcast(file_encoder,len(file_encoder),mpichar,0,mpicom)
#endif

   if ( l_conv2d ) then 
     !prepare data and run encoder for convolution 2d model 
     !note: transponse of E3SM(lon,lat,lev)->DeepONet(lat,lon,1,lev) 
     vmin = real(minval(vari),kind=r4)
     vmax = real(maxval(vari),kind=r4)
     do j = 1, ny 
       do i = 1, nx
         !normalize the data with max/min
         encin4d(i,j,1,:) = 2.0_r4*(real(vari(i,j,:),kind=r4)-vmin)/(vmax-vmin) - 1.0_r4
       end do
     end do

     !call deepOnet Encoder  
     call input_tensors%create
     call input_tensors%add_array(encin4d)
     call torch_mod%load(trim(file_path)//trim(file_encoder))
     call torch_mod%forward(input_tensors,out_tensor)
     call out_tensor%to_array(encout2d)  

     if (masterproc .and. l_print_diag) then
       write(iulog,*) 'shape of encin  = ',shape(encin4d)
       write(iulog,*) 'encin(min/max)  = ',minval(encin4d),maxval(encin4d)
       write(iulog,*) 'shape of encout = ',shape(encout2d)
       write(iulog,*) 'encout(min/max) = ',minval(encout2d),maxval(encout2d)
     end if

     !process to output array 
     do j = 1, ny1 
       do i = 1, nx1
         m = (j-1)*nx1 + i
         varo(i,j,1:nz1) = real(encout2d(m,1:nz1), kind = r8)
       end do 
     end do 

   else 

     !prepare data and run encoder for global model 
     encin3d(1:nx,1:ny,1:nz) = real(vari(1:nx,1:ny,1:nz),kind=r4) 

     !call deepOnet Encoder  
     call input_tensors%create
     call input_tensors%add_array(encin3d)
     call torch_mod%load(trim(file_path)//trim(file_encoder))
     call torch_mod%forward(input_tensors,out_tensor)
     call out_tensor%to_array(encout3d)  
     if (masterproc .and. l_print_diag) then
       write(iulog,*) 'shape of encin  = ',shape(encin3d)
       write(iulog,*) 'encin(min/max)  = ',minval(encin3d),maxval(encin3d)
       write(iulog,*) 'shape of encout = ',shape(encout3d)
       write(iulog,*) 'encout(min/max) = ',minval(encout3d),maxval(encout3d)
     end if

     !output array 
     varo(1:nx1,1:ny1,1:nz1) = real(encout3d(1:nx1,1:ny1,1:nz1), kind = r8)

   end if 

   return 
  end subroutine !mltbc_don_encoder 

  subroutine mltbc_don_decoder(l_conv2d,file_path,varname,nx,ny,nz,nx1,ny1,nz1,vari,varo)
  !===========================================================================
  ! SZ - 06/05/2023: This subroutine attempt to call auto decoder for 
  !                  the Machine Learning (ML) model,
  !===========================================================================
   use cam_abortutils  , only : endrun
   use shr_const_mod,    only : SHR_CONST_PI, SHR_CONST_REARTH

   implicit none
   logical, intent(in)          :: l_conv2d
   character(len=*), intent(in) :: file_path !Path to machine learning model files 
   character(len=*), intent(in) :: varname   !nudge variable
   integer, intent(in)          :: nx,ny,nz
   integer, intent(in)          :: nx1,ny1,nz1
   real(r8),intent(in)          :: vari(nx1,ny1,nz1)
   real(r8),intent(inout)       :: varo(nx,ny,nz)

   ! Arguments
   !-----------
   type(torch_module)           :: torch_mod
   type(torch_tensor_wrap)      :: input_tensors
   type(torch_tensor)           :: out_tensor

   ! Local values
   !----------------
   logical, parameter           :: l_print_diag = .false.
   integer                      :: i,j,n,m,k,ii,jj
   real(r8)                     :: dcdmin,dcdmax
   real(r4)                     :: dcdin2d(nx1*ny1,nz1)
   real(r4), pointer            :: dcdout2d(:,:)
   real(r4)                     :: dcdin3d(nx1,ny1,nz1)
   real(r4), pointer            :: dcdout3d(:,:,:)

   if (masterproc) then
     !check if machine learning pt file exist 
     file_decoder = trim(varname)//'_Decoder.pt'
     inquire(file=trim(file_path)//trim(file_decoder),exist=l_mltbc_decoder)
     if (.not. l_mltbc_decoder) then
       write(iulog,*) "ERROR: "//trim(file_path)//trim(file_decoder)//" not found!"
       call endrun('Machine Learning Nudging Error: model file not exist')
     end if
   end if
#ifdef SPMD
   call mpibcast(file_decoder,len(file_decoder),mpichar,0,mpicom)
#endif

   if (l_conv2d) then 
     !parameters to normalize/denormalize DeepONet prediction 
     if (trim(varname) == 'U') then
       dcdmax = 0.0013223718_r8
       dcdmin = -0.0010001023_r8
     else  ! "V"
       dcdmax = 0.0010921889_r8
       dcdmin = -0.0014089954_r8
     end if
     !prepare input data and run decoder 
     !need to convert to single precision 
     do j = 1, ny1
       do i = 1, nx1
         m = (j-1)*nx1 + i
         dcdin2d(m,1:nz1) = real(vari(i,j,1:nz1),kind=r4)
       end do 
     end do 
     !call decoder 
     call input_tensors%create
     call input_tensors%add_array(dcdin2d)
     call torch_mod%load(trim(file_path)//trim(file_decoder))
     call torch_mod%forward(input_tensors, out_tensor)
     call out_tensor%to_array(dcdout2d)
     if (masterproc .and. l_print_diag) then
       write(iulog,*) 'shape of dcdin  = ',shape(dcdin2d)
       write(iulog,*) 'dcdin(min/max)  = ',minval(dcdin2d),maxval(dcdin2d)
       write(iulog,*) 'shape of dcdout = ',shape(dcdout2d)
       write(iulog,*) 'dcdout(min/max) = ',minval(dcdout2d),maxval(dcdout2d)
     end if
     !denormalize ouput from decoder to convert it to tendency 
     do j = 1, ny 
       do i = 1, nx
         m = (j-1)*nx + i
         varo(i,j,1:nz) = 0.5_r8*(real(dcdout2d(m,1:nz),kind=r8)+1.0_r8)*(dcdmax-dcdmin)+dcdmin
       end do
     end do
   else 
     !input data 
     dcdin3d(1:nx1,1:ny1,1:nz1) = real(vari(1:nx1,1:ny1,1:nz1),kind=r4)
     !call decoder 
     call input_tensors%create
     call input_tensors%add_array(dcdin3d)
     call torch_mod%load(trim(file_path)//trim(file_decoder))
     call torch_mod%forward(input_tensors, out_tensor)
     call out_tensor%to_array(dcdout3d)
     if (masterproc .and. l_print_diag) then
       write(iulog,*) 'shape of dcdin  = ',shape(dcdin3d)
       write(iulog,*) 'dcdin(min/max)  = ',minval(dcdin3d),maxval(dcdin3d)
       write(iulog,*) 'shape of dcdout = ',shape(dcdout3d)
       write(iulog,*) 'dcdout(min/max) = ',minval(dcdout3d),maxval(dcdout3d)
     end if
     !output data 
     varo(1:nx,1:ny,1:nz) = real(dcdout3d(1:nx,1:ny,1:nz),kind=r8)
   end if 

   return 
  end subroutine !mltbc_don_decoder

  subroutine mltbc_gather_data(arr,nflds,ngcols,arro) !out
  !===========================================================================
  ! SZ - 06/05/2023: This subroutine attempt to collect the data at all chunks 
  !                  and convert to a sigle array
  !===========================================================================
   use ppgrid,        only: pcols, begchunk, endchunk
   use phys_grid,     only: gather_chunk_to_field
   use dyn_grid,      only: get_horiz_grid_dim_d, get_horiz_grid_d, get_dyn_grid_parm_real1d

   implicit none
   logical, parameter      :: l_print_diag = .false. 
   integer, intent(in)     :: ngcols,nflds  ! number of fields
   real(r8), intent(in)    :: arr(pcols,begchunk:endchunk,nflds) ! Input array, chunked
   real(r8), intent(inout) :: arro(ngcols,nflds) !Output array, global rectagular 

   integer :: i, j, n, ifld 
   integer :: ngtot  ! global column count (all)
   integer :: hdim1, hdim2  ! dimensions of rectangular horizontal 
                            !  grid data structure, If 1D data 
                            !  structure, then hdim2_d == 1.

   real(r8), allocatable :: arr_field(:,:,:)  ! rectangular version of arr

   call get_horiz_grid_dim_d(hdim1, hdim2)
   allocate(arr_field(hdim1,hdim2,nflds))

   arr_field(:,:,:) = 0.0_r8
   call gather_chunk_to_field (1, 1, nflds, hdim1, arr, arr_field)

   ngtot = hdim1*hdim2
   if ( ngtot /= ngcols ) then 
     write(iulog,*) 'Machine Learning Nudging Error: in/out size mismatch ngcols(in),ngcols(out) = ',ngtot,ngcols 
     call endrun ('Machine Learning Nudging Error: mltbc_gather_data failed')
   end if 

   !combine lat-lon array to ncol array
   do j = 1, hdim2
     do i = 1, hdim1
       n = (j-1)*hdim1 + i
       arro(n,:) = arr_field(i,j,:) 
     end do 
   end do 
 
   if (masterproc.and.l_print_diag) then
     write(iulog,*) 'shape of arr(in)  = ',shape(arr)
     write(iulog,*) 'min/max arr(in)   = ',minval(arr),maxval(arr)
     write(iulog,*) 'shape of arr(out) = ',shape(arro)
     write(iulog,*) 'min/max arr(out)  = ',minval(arro),maxval(arro)
   end if

   deallocate(arr_field)

   return
  end subroutine !mltbc_gather_data

  subroutine mltbc_gather_patch(vname,nlev,ngcol,nlon,nlat,model_state,out_state) !out
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!Following subroutine is used to gather global data and process it to      !!
  !!Model space input to Machine Learning model                      !! 
  !!Author: Shixuan Zhang (shixuan.zhang@pnnl.gov)                            !! 
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
   use ppgrid,           only : pver,pverp,pcols,begchunk,endchunk
   use phys_grid,        only : get_ncols_p, scatter_field_to_chunk 
   use cam_abortutils  , only : endrun
   use shr_const_mod,    only : SHR_CONST_PI, SHR_CONST_REARTH
   use cam_history  ,    only : outfld

   implicit none
   logical, parameter           :: l_print_diag = .false. 
   character(len=*), intent(in) :: vname   !nudge variable
   integer, intent(in)          :: nlev,ngcol,nlon,nlat 
   real(r8), intent(in)         :: model_state(ngcol,nlev)
   real(r8), intent(inout)      :: out_state(nlon,nlat,nlev)

   ! Local values
   !----------------
   integer  :: lchnk,ncol
   integer  :: i,j,k,ii,jj,m,n 
   real(r8) :: sum_x(nlev)
   real(r8), allocatable :: test_rgd(:,:)
   real(r8), allocatable :: test_dif(:,:)
   real(r8), allocatable :: ftem(:,:,:)
   real(r8), allocatable :: ftem2(:,:)

   ! initialize output array as fillvalue 
   out_state(:,:,:) = fillvalue

   !Interpolate to regional data 
   if (mltbc_patch_biln) then
     !bilinear interpolation      
     do j = 1, nlat
       do i = 1, nlon
         sum_x(:) = 0.0_r8
         ii = 0
         m = (j-1)*nlon + i
         do jj = 1, 4
           n = mltbc_se2latlon_ind(m,jj)
           if (mltbc_se2latlon_wgt(m,jj) == 0.0_r8) ii = ii + 1
           sum_x(:) = sum_x(:) + model_state(n,:) * mltbc_se2latlon_wgt(m,jj)
         end do
         if (ii /= 4) then
           out_state(i,j,:) = sum_x
         end if
       end do
     end do
   else 
     !neareast to destination 
     do j = 1, nlat
       do i = 1, nlon
         m = (j-1)*nlon + i
         n = mltbc_se2latlon_ind(m,5)
         out_state(i,j,:) = model_state(n,:) 
       end do
     end do 
   end if 

   if (masterproc.and.l_print_diag) then
     write(iulog,*) 'shape of state(in)  = ',shape(model_state)
     write(iulog,*) 'min/max state(in)   = ',minval(model_state),maxval(model_state)
     write(iulog,*) 'shape of state(out) = ',shape(out_state)
     write(iulog,*) 'min/max state(out)  = ',minval(out_state),maxval(out_state)
   end if

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
   !!Begin diagnose interpolation!! 
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   if (mltbc_interp_test) then 
     allocate (test_dif(ngcol,nlev))
     allocate (test_rgd(nlon*nlat,nlev))
     test_dif(:,:) = 0.0_r8
     !interpolated back to model grid and check  
     if (mltbc_patch_biln) then
       !method 1: linear interpolation to destination grid
       do j = 1, nlat
         do i = 1, nlon
           m = (j-1)*nlon+i
           test_rgd(m,:) = out_state(i,j,:)
         end do
       end do 
       do n = 1, ngcol
         sum_x(:) = 0.0_r8
         ii = 0
         do jj = 1, 4
           m = mltbc_latlon2se_ind(n,jj)
           if(mltbc_latlon2se_wgt(n,jj) == 0.0_r8) ii = ii + 1
           sum_x(:) = sum_x(:) + mltbc_latlon2se_wgt(n,jj)*test_rgd(m,:)
         end do
         if( ii /= 4 ) then
           test_dif(n,:) = model_state(n,:) - sum_x(:)
         end if
       end do   
       end do
     else 
       !method 2: nearest to destination grid
       do j = 1, nlat
         do i = 1, nlon
           m = (j-1)*nlon + i
           n = mltbc_se2latlon_ind(m,5)
           test_dif(n,:) = (model_state(n,:) - out_state(i,j,:))
         end do
       end do
     end if
     !output Diagnostics 
     !scatter filed to chunk (from wrk_tend to nudging_tend) 
     if (nlev > 1) then 
       allocate (ftem(pcols,pver,begchunk:endchunk))
       call scatter_field_to_chunk(1,nlev,1,ngcol,test_dif(:,:),ftem)
       do lchnk=begchunk,endchunk
         ncol=get_ncols_p(lchnk)
         if (trim(vname) == 'U') then 
           call outfld('DON_RGD_UERR',  ftem(:,:,lchnk), pcols,lchnk)
         end if 
         if (trim(vname) == 'V') then
           call outfld('DON_RGD_VERR',  ftem(:,:,lchnk), pcols,lchnk)
         end if
         if (trim(vname) == 'T') then
           call outfld('DON_RGD_TERR',  ftem(:,:,lchnk), pcols,lchnk)
         end if
         if (trim(vname) == 'Q') then
           call outfld('DON_RGD_QERR',  ftem(:,:,lchnk), pcols,lchnk)
         end if
       end do 
       deallocate(ftem)
     else
       allocate (ftem2(pcols,begchunk:endchunk))      
       call scatter_field_to_chunk(1,nlev,1,ngcol,test_dif(:,1),ftem2) 
       do lchnk = begchunk,endchunk
         ncol=get_ncols_p(lchnk)
         if (trim(vname) == 'PS') then
           call outfld('DON_RGD_PSERR',  ftem2(:,lchnk), pcols,lchnk)
         end if
       end do 
       deallocate(ftem2)
     end if 
     deallocate(test_dif)
     deallocate(test_rgd)

   end if
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
   !!End diagnose interpolation  !! 
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   !Sanity check
   if (any(out_state(:,:,:) == fillvalue)) then 
      write(iulog,*) 'Machine Learning Nudging Error: missing value appears'
      call endrun('Machine Learning Nudging Error: working data contains missing values')
   end if 
  end subroutine !mltbc_gather_patch

  subroutine update_land_nudging_tend(ncol, umod, vmod, tvmod, qmod,    & ! In 
                                      ubobs, vbobs, tbobs, tdbobs,      & ! In  
                                      qbobs, pmid_obs, sfac, ndg_srf_q, & ! In 
                                      udt, vdt, tdt, qdt )                ! Out 
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!Following subroutine is used to apply nudging at the surface layer only   !! 
  !! Author: Shixuan Zhang (shixuan.zhang@pnnl.gov)                           !! 
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
  use ppgrid,        only  : pver,pverp,pcols
  use physconst,     only  : rga, cpair, gravit, rair, latvap, rh2o, zvir, cappa
  use cam_abortutils,only  : endrun
  use wv_saturation, only  : qsat, qsat_water, svp_ice


  integer,  intent(in) :: ncol    ! number of columns
  logical,  intent(in) :: ndg_srf_q

  real(r8), intent(in) :: sfac(pcols)
  real(r8), intent(in) :: umod(pcols,pver)
  real(r8), intent(in) :: vmod(pcols,pver)
  real(r8), intent(in) :: tvmod(pcols,pver)
  real(r8), intent(in) :: qmod(pcols,pver)

  real(r8), intent(in) :: pmid_obs(pcols,pver)
  real(r8), intent(in) :: ubobs(pcols) ! 10m wind 
  real(r8), intent(in) :: vbobs(pcols) ! 10m wind 
  real(r8), intent(in) :: tbobs(pcols) ! 2m temperature 
  real(r8), intent(in) :: tdbobs(pcols)! 2m dewpoint tempeature 
  real(r8), intent(in) :: qbobs(pcols) ! 2m humidity

  real(r8), intent(inout) :: udt(pcols,pver)
  real(r8), intent(inout) :: vdt(pcols,pver)
  real(r8), intent(inout) :: tdt(pcols,pver)
  real(r8), intent(inout) :: qdt(pcols,pver)

  !Local variable 

  !use dewpoint temperature to derive specific humidity at 2-m 
  !because ERA5 reanalysis only provides dewpoint temperature at 2-m  
  logical, parameter :: l_q2_from_td = .true.

  !adjust supersaturation 
  logical, parameter :: l_adj_super_saturation = .false.

  integer  :: i, k, m

  real(r8) :: wuprof(pcols,pver)
  real(r8) :: wvprof(pcols,pver)
  real(r8) :: wtprof(pcols,pver)
  real(r8) :: wqprof(pcols,pver)

  real(r8) :: qsbobs(pcols)
  real(r8) :: esbobs(pcols)
  real(r8) :: tvbobs(pcols)
  real(r8) :: qbref(pcols)

  !specify the vertical weighting 
  do i = 1, ncol
    do k = 1, pver
     if( k < pver ) then
       wuprof(i,k) = 0._r8
       wvprof(i,k) = 0._r8
       wtprof(i,k) = 0._r8
       wqprof(i,k) = 0._r8
     else
       wuprof(i,k) = 1.0_r8
       wvprof(i,k) = 1.0_r8
       wtprof(i,k) = 1.0_r8
       wqprof(i,k) = 1.0_r8
     end if
    end do
  end do

  if ( l_q2_from_td ) then
    !derive the specific humidity from dewpoint temperature 
    !see Eq(11) in LAWRENCE(2005,MWR)
    do i = 1, ncol
      call qsat(tbobs(i), pmid_obs(i,pver), esbobs(i), qsbobs(i))
      qbref(i) = qsbobs(i) * exp((1.0_r8 - tbobs(i)/tdbobs(i))*latvap/rh2o/tbobs(i))
    end do
  else
    !limit humidity in the observational data 
    do i = 1, ncol
      qbref(i)  = max(1.E-8_r8, qbobs(i))
      ! adjust supersaturation
      call qsat(tbobs(i), pmid_obs(i,pver), esbobs(i), qsbobs(i))
      if(l_adj_super_saturation .and. (qbobs(i) > qsbobs(i))) then
         qbref(i) = qsbobs(i)
      end if
    end do
  end if

  !convert the obs temperture to virtual temperature 
  do i = 1, ncol
    tvbobs(i) = tbobs(i) * (1.0_r8 + zvir * qbref(i))
  end do

  !calcualte the nudging tendency 
  do k = 1, pver
    do i = 1, ncol
      udt(i,k) = udt(i,k) + (ubobs(i) - umod(i,pver))*wuprof(i,k)*sfac(i)
      vdt(i,k) = vdt(i,k) + (vbobs(i) - vmod(i,pver))*wvprof(i,k)*sfac(i)
      tdt(i,k) = tdt(i,k) + (tvbobs(i) - tvmod(i,pver))*wtprof(i,k)*sfac(i) &
                            /(1.0_r8 + zvir * qmod(i,pver))
      if (ndg_srf_q) then
        qdt(i,k) = qdt(i,k) + (qbref(i) - qmod(i,pver))*wqprof(i,k)*sfac(i)
      end if
    end do
  end do

  return 
  end subroutine !update_land_nudging_tend

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!Following subroutine is used to modify the weight profiles of nudging tendency with
  !!(1) a linear relaxation method to linearly relax the weight from 1 to 0 at upper   
  !!    model levels so that weaker nudging is applied. This is becase the
  !!    reanalysis and E3SM show large discrepencies in these layers                    
  !!(2) a method to exclude the nudging tendencies within boundary layers as the
  !!    nudging are found to significantly interfere model climate when the
  !!    moisture is nudged within the boundary layers.  
  !! Author: Shixuan Zhang (shixuan.zhang@pnnl.gov)   
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
  subroutine update_nudge_prof (ncol, zm_mod, zi_mod, pmid_mod, pblh, & !In 
                                use_pdep_nudge, use_upp_lrelx,        & !In 
                                no_pbl_uv, no_pbl_t, no_pbl_q,        & !In 
                                ufac, vfac, tfac, qfac)                 !Out 

  use ppgrid,        only  : pver,pverp,pcols
  use physconst,     only  : rga, cpair, gravit, rair, latvap, rh2o, zvir, cappa
  use cam_abortutils,only  : endrun

  logical, intent(in)  :: use_upp_lrelx
  logical, intent(in)  :: use_pdep_nudge

  integer, intent(in)  :: ncol    
  integer, intent(in)  :: no_pbl_uv
  integer, intent(in)  :: no_pbl_t
  integer, intent(in)  :: no_pbl_q

  real(r8), intent(in) :: pmid_mod(pcols,pver) 
  real(r8), intent(in) :: zi_mod(pcols,pverp) 
  real(r8), intent(in) :: zm_mod(pcols,pverp) 
  real(r8), intent(in) :: pblh(pcols)

  real(r8), intent(inout) :: ufac(pcols,pver)
  real(r8), intent(inout) :: vfac(pcols,pver)
  real(r8), intent(inout) :: tfac(pcols,pver)
  real(r8), intent(inout) :: qfac(pcols,pver)

  !local variables 
  integer  :: i, k, m
  real(r8) :: kpblt(pcols)
  real(r8) :: wuprof(pcols,pver)
  real(r8) :: wvprof(pcols,pver)
  real(r8) :: wtprof(pcols,pver)
  real(r8) :: wqprof(pcols,pver)
 
  ! initialize kpbl and weigthing profiles 
  do i = 1, ncol
    kpblt(i) = pver
    do k = 1,pver
      wuprof(i,k) = 0.0_r8 
      wvprof(i,k) = 0.0_r8
      wtprof(i,k) = 0.0_r8
      wqprof(i,k) = 0.0_r8
    end do 
  end do

  !Apply scaling such that nudging strength falls off with pressure.
  if (use_pdep_nudge) then
    do k = 1, pver
      do i = 1, ncol
        ufac(i,k) = ufac(i,k) * (pmid_mod(i,k) / pmid_mod(i,pver))
        vfac(i,k) = vfac(i,k) * (pmid_mod(i,k) / pmid_mod(i,pver))
        tfac(i,k) = tfac(i,k) * (pmid_mod(i,k) / pmid_mod(i,pver))
        qfac(i,k) = qfac(i,k) * (pmid_mod(i,k) / pmid_mod(i,pver))
      end do
    end do
  end if
   
  ! Add a linear relexation of the nudging tendency on the upper layer 
  if (use_upp_lrelx) then
    do i = 1, ncol
      do k = pver, 1, -1
        if ( pmid_mod(i,k) < p_uv_relax ) then
          ufac(i,k) = ufac(i,k) * max(0.01_r8, pmid_mod(i,k)/p_uv_relax)
          vfac(i,k) = vfac(i,k) * max(0.01_r8, pmid_mod(i,k)/p_uv_relax)
        end if
        if ( pmid_mod(i,k) < p_T_relax ) then
          tfac(i,k) = tfac(i,k) * max(0.01_r8, pmid_mod(i,k)/p_T_relax)
        end if
        if ( pmid_mod(i,k) < p_q_relax ) then
          qfac(i,k) = qfac(i,k) * max(0.01_r8, pmid_mod(i,k)/p_q_relax)
        end if
        if ( pmid_mod(i,k) < p_norelax ) then
          ufac(i,k) = 0._r8
          vfac(i,k) = 0._r8
          tfac(i,k) = 0._r8
          qfac(i,k) = 0._r8
        end if
      end do
    end do
  end if

  !Exclude the nudging tendency in boundary layer
  !Note: PBL height is in AGL  
  do k = pver-1,1,-1
     do i = 1,ncol
        if (abs(zm_mod(i,k)-pblh(i)) < (zi_mod(i,k)-zi_mod(i,k+1))*0.5_r8) kpblt(i) = k
     end do
  end do

  !--------------------------------------------------
  ! Strategy to exclude the nudging at near surface 
  !--------------------------------------------------
  do i = 1, ncol
    do k = pver, 1, -1
      select case (no_pbl_uv)
         case (0)
           wuprof(i,k) = 1.0_r8
           wvprof(i,k) = 1.0_r8
         case (1)
           wuprof(i,k) = 0.5_r8 * (1.0_r8 + tanh((real(kpblt(i))-real(k))/0.1_r8))
           wvprof(i,k) = 0.5_r8 * (1.0_r8 + tanh((real(kpblt(i))-real(k))/0.1_r8))
         case (2)
           wuprof(i,k) = 0.5_r8 * (1.0_r8 + tanh((zm_mod(i,k)-z_min)/(0.1_r8*z_min)))
           wvprof(i,k) = 0.5_r8 * (1.0_r8 + tanh((zm_mod(i,k)-z_min)/(0.1_r8*z_min)))
         case default
         call endrun('nudging_tend error: invalid option for no_pbl_uv nudging')
      end select

      select case (no_pbl_t)
         case (0)
           wtprof(i,k) = 1.0_r8
         case (1)
           wtprof(i,k) = 0.5_r8 * (1.0_r8 + tanh((real(kpblt(i))-real(k))/0.1_r8))
         case (2)
           wtprof(i,k) = 0.5_r8 * (1.0_r8 + tanh((zm_mod(i,k)-z_min)/(0.1_r8*z_min)))
         case default
         call endrun('nudging_tend error: invalid option for no_pbl_t nudging')
      end select

      select case (no_pbl_q)
         case (0)
           wqprof(i,k) = 1.0_r8
         case (1)
           wqprof(i,k) = 0.5_r8 * (1.0_r8 + tanh((real(kpblt(i))-real(k))/0.1_r8))
         case (2)
           wqprof(i,k) = 0.5_r8 * (1.0_r8 + tanh((zm_mod(i,k)-z_min)/(0.1_r8*z_min)))
         case default
           call endrun('nudging_tend error: invalid option for no_pbl_q nudging')
      end select
    end do
  end do

  ! apply the weighting on the nudging coefficient 
  do i = 1, ncol
    do k = pver, 1, -1
      ufac(i,k) = ufac(i,k) * wuprof(i,k)
      vfac(i,k) = vfac(i,k) * wvprof(i,k)
      tfac(i,k) = tfac(i,k) * wtprof(i,k)
      qfac(i,k) = qfac(i,k) * wqprof(i,k)
    end do
  end do

  return 
  end subroutine !update_nudge_prof

  !!!!!!!!!!!!!!!!!!!!!!
  !ps nudigng subroutine 
  !!!!!!!!!!!!!!!!!!!!!!
  subroutine ps_nudging(ncol, dtime, use_ps_adj, zi_obs, zm_obs, &
                        pint_obs, pint_mod, pmid_obs, pmid_mod, &
                        lnpint_obs, lnpint_mod, lnpmid_obs, lnpmid_mod, &
                        pdel_obs, pdel_mod, rpdel_obs, rpdel_mod, &
                        uref, vref, tref, qref, tobs, psref, &
                        phis_obs, phis_mod, psobs, psmod, psfac, psdt )

  use physconst,     only  : rga, cpair, gravit, rair, zvir, cappa
  use shr_vmath_mod, only  : shr_vmath_log
  use hycoef,        only  : hycoef_init, hyam, hybm, hyai, hybi, ps0
  use ppgrid,        only  : pver,pverp,pcols

  logical,  intent(in) :: use_ps_adj
  integer,  intent(in) :: ncol
  real(r8), intent(in) :: dtime

  real(r8), intent(in) :: phis_mod(pcols)
  real(r8), intent(in) :: phis_obs(pcols)
  real(r8), intent(in) :: tobs(pcols,pver)
  real(r8), intent(in) :: zi_obs(pcols,pverp)
  real(r8), intent(in) :: zm_obs(pcols,pver)

  real(r8), intent(in) :: psobs(pcols)
  real(r8), intent(in) :: pint_obs(pcols,pverp)
  real(r8), intent(in) :: pmid_obs(pcols,pver)
  real(r8), intent(in) :: pdel_obs(pcols,pver)
  real(r8), intent(in) :: rpdel_obs(pcols,pver)
  real(r8), intent(in) :: lnpint_obs(pcols,pverp)
  real(r8), intent(in) :: lnpmid_obs(pcols,pver)

  real(r8), intent(inout) :: psmod(pcols)
  real(r8), intent(inout) :: psfac(pcols)
  real(r8), intent(inout) :: pint_mod(pcols,pverp)
  real(r8), intent(inout) :: pmid_mod(pcols,pver)
  real(r8), intent(inout) :: pdel_mod(pcols,pver)
  real(r8), intent(inout) :: rpdel_mod(pcols,pver)
  real(r8), intent(inout) :: lnpint_mod(pcols,pverp)
  real(r8), intent(inout) :: lnpmid_mod(pcols,pver)

  real(r8), intent(inout) :: psdt(pcols)
  real(r8), intent(inout) :: psref(pcols)
  real(r8), intent(inout) :: uref(pcols,pver)
  real(r8), intent(inout) :: vref(pcols,pver)
  real(r8), intent(inout) :: tref(pcols,pver)
  real(r8), intent(inout) :: qref(pcols,pver)

  !local variables 
  !Parameters to adjust surface pressure 
  !Reference: !https://github.com/NOAA-GFDL/GFDL_atmos_cubed_sphere/blob/main/tools/fv_clim_nudge.F90
  real(r8) :: dtdz       = -6.5E-3_r8 ! Lapse rate for adjustment of PS
  real(r8) :: t_ref1     = 290.5_r8   ! reference temperature for adjustment of PS
  real(r8) :: t_ref2     = 255.0_r8   ! reference temperature for adjustment of PS
  real(r8) :: dz_thres   = 1.E-3_r8   ! threshold of topography difference for adjustment of PS 

  real(r8), parameter :: z_min      = 150._r8

  integer  :: i, k, m, kk
  real(r8) :: x
  real(r8) :: del_phis, tbot, pbot, tmp
  real(r8) :: tsurf, lapse, t0
  real(r8) :: tvobs(pcols,pver)
  real(r8) :: zobs(pcols,pverp)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  do i = 1,ncol
    psref(i) = psobs(i) 
  end do 

  !adjust the nudging coefficient 
  do i = 1, ncol
    if ( abs(psmod(i)-psref(i)) > 2.E2_r8 ) then
      psfac(i) = psfac(i) * 2.E2_r8 / abs(psmod(i)-psref(i))
    else
      psfac(i) = psfac(i) * 1.0_r8
    end if
  enddo

  !surface pressure adjustment due to topography difference 
  !adjust observed ps to model topography (if the pre-processing does not appy it)
  !Following the method developed by ECMWF  
  if ( use_ps_adj ) then

    do i = 1,ncol

      del_phis = phis_obs(i) - phis_mod(i)

      if(abs(del_phis) > dz_thres) then 
        !Tbot and Pbot are determined from the first model level that is at
        !least 150m above the surface
 
        kk = pver
        do k = pver-1,1,-1
          if (abs(zm_obs(i,k)-z_min) < (zi_obs(i,k)-zi_obs(i,k+1))*0.5_r8) kk = k
        end do 
 
        lapse = -dtdz
        k     = kk
        ! Define Tbot & Pbot
        tbot  = tobs(i,k)
        pbot  = pmid_obs(i,k)
        tmp   = lapse*(rair/gravit)*(psobs(i)/pbot - 1.0_r8)
        tsurf = tbot*(1.0_r8 + tmp)
        t0    = tsurf + lapse*phis_obs(i)/gravit
 
        if (t0 .gt. t_ref1 .and. tsurf .le. t_ref1) then
          lapse = (t_ref1 - tsurf)*gravit/phis_obs(i)
        else if (t0 .gt. t_ref1 .and. tsurf .gt. t_ref1) then
          lapse = 0._r8
          tsurf = (t_ref1 + tsurf)*0.5_r8
        end if
 
        if (tsurf .lt. t_ref2) then
          lapse = -dtdz
          tsurf = (t_ref2 + tsurf)*0.5_r8
        end if
  
        x   = lapse*del_phis/(gravit*tsurf)
        tmp = 1._r8 - x/2._r8 + x**2._r8/3._r8
        tmp = del_phis/(rair*tsurf)*tmp
        psref(i) = psobs(i)*exp(tmp)
      end if

    end do 

  end if 
 
  do i = 1,ncol
    psdt(i) = psref(i) - psmod(i)
    !limit the ps tendency 
    psdt(i) = sign ( min(10.E2, abs(psdt(i))), psdt(i) )
    psdt(i) = psdt(i) * psfac(i)
  end do

  !convert the values to conserve the geopotential height 
  do i=1,ncol
    do k=1,pver
      uref(i,k) = uref(i,k) * pdel_mod(i,k)
      vref(i,k) = vref(i,k) * pdel_mod(i,k)
      tref(i,k) = tref(i,k) * (lnpint_mod(i,k+1) - lnpint_mod(i,k))
      qref(i,k) = qref(i,k) * pdel_mod(i,k)
    enddo
  end do 

  !Update ps
  do i = 1,ncol
    psmod(i) = psmod(i) + psdt(i)*dtime
  end do 

  !Pdel etc.
  do i = 1,ncol
    do k = 1,pver
      pint_mod(i,k)     = hyai(k)*ps0+hybi(k)*psmod(i)
      lnpint_mod(i,k)   = log(pint_mod(i,k))
    end do 
    pint_mod(i,pverp)   = hyai(pverp)*ps0+hybi(pverp)*psmod(i)
    lnpint_mod(i,pverp) = log(pint_mod(i,pverp))
    do k = 1,pver
      pdel_mod(i,k)     = pint_mod(i,k+1) - pint_mod(i,k) 
      rpdel_mod(i,k)    = 1._r8/pdel_mod(i,k)
    end do
  end do 

 !Convert back to model quantities 
  do i = 1,ncol
    do k=1,pver
      uref(i,k) = uref(i,k) / pdel_mod(i,k)
      vref(i,k) = vref(i,k) / pdel_mod(i,k)
      tref(i,k) = tref(i,k) / (lnpint_mod(i,k+1) - lnpint_mod(i,k))
      qref(i,k) = qref(i,k) / pdel_mod(i,k)
    enddo
  end do

  return
  end subroutine !ps_nudging

  real(r8) function dist_latlon_2pts(rlat1,rlon1,rlat2,rlon2,PI)
   !
   ! dist_latlon_2pts: for the given lat and lon at 2 locations, set the 
   !                   great circle distance. 
   !===============================================================

   implicit none
   ! Arguments
   !--------------
   real(r8) :: rlat1,rlon1
   real(r8) :: rlat2,rlon2
   real(r8) :: PI

   ! Local values
   !----------------
   real(r8) :: half_PI, two_PI, deg2rad
   real(r8) :: dlon, rtemp

   half_PI = 0.5_r8 * PI
   two_PI  = 2.0_r8 * PI
   deg2rad = PI/180._r8 

   !calculate the great circle difference (assume radius of earth = 1)
   if ((rlat1 /= rmissing) .and. (rlon1 /= rmissing) .and. (rlat2 /= rmissing) .and. (rlon2 /= rmissing)) then
     dlon = rlon1 - rlon2 
     if (abs(rlat1) >= half_PI .or. abs(rlat2) >= half_PI .or. dlon == 0.0_r8) then
       dist_latlon_2pts = abs(rlat1 - rlat2)
     else
       rtemp = sin(rlat2) * sin(rlat1) + &
               cos(rlat2) * cos(rlat1) * cos(dlon)
       if (rtemp < -1.0_r8) then
         dist_latlon_2pts = PI
       else if (rtemp > 1.0_r8) then
         dist_latlon_2pts = 0.0_r8
       else
         dist_latlon_2pts = acos(rtemp)
       endif
     end if
   else 
     dist_latlon_2pts = fillvalue 
    !call endrun('dist_latlon_2pts:: lat/lon contains missing values')
   end if

   ! End Routine
   !------------
   return
  end function ! dist_latlon_2pts 

  subroutine latlon2se_interp_init(ngcols, rlon_d, rlat_d, &
                                   nrlon, nrlat, reclon, reclat, &
                                   close_ind, weight_qdl)
   !---------------------------------------------------------------------------
   ! This program computes weighting functions to map a variable on 
   ! (nlat,nlon) resolution to global SE grid 
   ! Author: Shixuan Zhang  -- August 2023
   !---------------------------------------------------------------------------
   use shr_kind_mod,    only : r8 => shr_kind_r8
   use shr_const_mod,   only : SHR_CONST_PI, SHR_CONST_REARTH
   use dycore,          only : dycore_is

   implicit none
   integer, parameter   :: ncorner = 4
   real(r8), parameter  :: re      = SHR_CONST_REARTH
   real(r8), parameter  :: PI      = SHR_CONST_PI
   real(r8), parameter  :: deg2rad = SHR_CONST_PI/180._r8

   integer,  intent(in) :: nrlat, nrlon, ngcols
   real(r8), intent(in) :: rlon_d(ngcols) !lon for model grid in radians 
   real(r8), intent(in) :: rlat_d(ngcols) !lat for model grid in radians  
   real(r8), intent(in) :: reclon(nrlon)
   real(r8), intent(in) :: reclat(nrlat)

   integer,  allocatable, intent(out) :: close_ind(:,:)
   real(r8), allocatable, intent(out) :: weight_qdl(:,:)

   real(r8), allocatable :: close_lat(:,:)
   real(r8), allocatable :: close_lon(:,:)
   real(r8), allocatable :: close_dst(:,:)
   real(r8), allocatable :: trlon(:), olon(:)
   real(r8), allocatable :: trlat(:), olat(:)

   ! local variables
   logical  :: neg_root
   integer  :: i, j, n, m, k
   integer  :: iorg, ntcs, nmin, nrecg
   integer  :: oc(1)
   integer  :: sh_corn(ncorner)
   real(r8) :: sh_rlat(ncorner), sh_rlon(ncorner), sh_rdst(ncorner)
   real(r8) :: rlat, rlon, dlat, dlon, ddlon, rtemp
   real(r8) :: maxdist, dmin, this_dist
   real(r8) :: dx,dy, p, q

   if (.not. dycore_is('UNSTRUCTURED')) then
      call endrun ('latlon2se_interp_init ERROR: only works for  SE dycore')
   end if

   !maxdistance for the search 
   maxdist = 1.2_r8*(21600.0_r8/ngcols)*deg2rad

   !initialize arrary for interpolation 
   allocate(close_lat(ngcols,5))
   allocate(close_lon(ngcols,5))
   allocate(close_dst(ngcols,5))
   allocate(weight_qdl(ngcols,5))
   allocate(close_ind(ngcols,5))
   close_lat(:,:)  = rmissing
   close_lon(:,:)  = rmissing
   close_dst(:,:)  = rmissing
   weight_qdl(:,:) = rmissing
   close_ind(:,:)  = imissing

   nrecg   = nrlat * nrlon
   allocate(trlat(nrecg))
   allocate(trlon(nrecg))
   do j = 1, nrlat
     do i = 1, nrlon
       m = (j-1)*nrlon + i
       trlat(m) = reclat(j)
       trlon(m) = reclon(i)
     end do
   end do

   !find the cell-center grid point closest to machine learning lat-lon grid
   do n = 1, ngcols
     !use a temporary array olat/olon so that we can modify values 
     allocate(olat(nrecg))
     allocate(olon(nrecg))
     olat  = trlat
     olon  = trlon
     rlat  = rlat_d(n)
     rlon  = rlon_d(n)
     !Find center and 4 corner values that are closed to target grid 
     !-----------------------------------------------------
     !locate the model grid closest to the lat-lon grid
     !corners --(4,4)--         (i+1,j+1)
     !corners   (i,j+1)       --(3,3)--
     !corners        --(xp,yp)-- 
     !corners --(1,1)--         (i+1,j) 
     !corners   (i,j)         --(2,2)--
     do k = 5,1,-1 ! 4 corners + 1 center location 
       if (k < 5 ) then !nearest 
         olat(:) = rmissing
         olon(:) = rmissing
         dlat    = close_lat(n,5) - rlat
         dlon    = close_lon(n,5) - rlon
         if (dlat <= 0._r8 .and. dlon <= 0._r8) then
           if (k == 1) then !(i,j), bottom left side of rlat,rlon
             olat = trlat
             olon = trlon
           end if
           if (k == 2) then !(i+1,j), bottom right side of rlat,rlon
             olat = merge(trlat,olat,trlat<min(close_lat(n,3),rlat).and.trlon>max(close_lon(n,5),rlon))
             olon = merge(trlon,olon,trlat<min(close_lat(n,3),rlat).and.trlon>max(close_lon(n,5),rlon))
           end if
           if (k == 3) then !(i+1,j+1), top right side of rlat,rlon
             olat = merge(trlat,olat,trlat>max(close_lat(n,5),rlat).and.trlon>max(close_lon(n,4),rlon))
             olon = merge(trlon,olon,trlat>max(close_lat(n,5),rlat).and.trlon>max(close_lon(n,4),rlon))
           end if
           if (k == 4) then !(i,j+1), top left side of rlat, rlon
             olat = merge(trlat,olat,trlat>max(close_lat(n,5),rlat).and.trlon<=rlon)
             olon = merge(trlon,olon,trlat>max(close_lat(n,5),rlat).and.trlon<=rlon)
           end if
         else if (dlat <= 0._r8 .and. dlon > 0._r8) then
           if (k == 1) then !(i,j)
             olat = merge(trlat,olat,trlat<min(close_lat(n,4),rlat).and.trlon<min(close_lon(n,2),rlon))
             olon = merge(trlon,olon,trlat<min(close_lat(n,4),rlat).and.trlon<min(close_lon(n,2),rlon))
           end if
           if (k == 2) then !(i+1,j)
             olat = trlat
             olon = trlon
           end if
           if (k == 3) then !(i+1,j+1)
             olat = merge(trlat,olat,trlat>max(close_lat(n,5),rlat).and.trlon>max(close_lon(n,4),rlon))
             olon = merge(trlon,olon,trlat>max(close_lat(n,5),rlat).and.trlon>max(close_lon(n,4),rlon))
           end if
           if (k == 4) then !(i,j+1)
             olat = merge(trlat,olat,trlat>max(close_lat(n,5),rlat).and.trlon<min(close_lon(n,5),rlon))
             olon = merge(trlon,olon,trlat>max(close_lat(n,5),rlat).and.trlon<min(close_lon(n,5),rlon))
           end if
         else if (dlat > 0._r8 .and. dlon > 0._r8) then
           if (k == 1) then !(i,j)
             olat = merge(trlat,olat,trlat<min(close_lat(n,4),rlat).and.trlon<min(close_lon(n,2),rlon))
             olon = merge(trlon,olon,trlat<min(close_lat(n,4),rlat).and.trlon<min(close_lon(n,2),rlon))
           end if
           if (k == 2) then !(i+1,j)
             olat = merge(trlat,olat,trlat<min(close_lat(n,3),rlat).and.trlon>max(close_lon(n,4),rlon))
             olon = merge(trlon,olon,trlat<min(close_lat(n,3),rlat).and.trlon>max(close_lon(n,4),rlon))
           end if
           if (k == 3) then !(i+1,j+1)
             olat = trlat
             olon = trlon
           end if
           if (k == 4) then !(i,j+1), top left side of rlat, rlon
             olat = merge(trlat,olat,trlon<min(close_lon(n,5),rlon).and.trlat>=rlat)
             olon = merge(trlon,olon,trlon<min(close_lon(n,5),rlon).and.trlat>=rlat)
           end if
         else if (dlat > 0._r8 .and. dlon <= 0._r8) then
           if (k == 1) then !(i,j)
             olat = merge(trlat,olat,trlat<min(close_lat(n,4),rlat).and.trlon<min(close_lon(n,2),rlon))
             olon = merge(trlon,olon,trlat<min(close_lat(n,4),rlat).and.trlon<min(close_lon(n,2),rlon))
           end if
           if (k == 2) then !(i+1,j)
             olat = merge(trlat,olat,trlat<min(close_lat(n,3),rlat).and.trlon>max(close_lon(n,4),rlon))
             olon = merge(trlon,olon,trlat<min(close_lat(n,3),rlat).and.trlon>max(close_lon(n,4),rlon))
           end if
           if (k == 3) then !(i+1,j+1)
             olat = merge(trlat,olat,trlon>max(close_lon(n,4),rlon).and.trlat>=rlat)
             olon = merge(trlon,olon,trlon>max(close_lon(n,4),rlon).and.trlat>=rlat)
           end if
           if (k == 4) then !(i,j+1)
             olat = trlat
             olon = trlon
           end if
         else
           call endrun('se2latlon_interp_init: ERROR in buidling the coordinates')
         end if
       end if

       !loop model grid to locate nearest point 
       ntcs  = 0
       nmin  = 0
       dmin  = fillvalue ! initialize a large value 
       do j = 1, nrlat
         do i = 1, nrlon
           m = (j-1)*nrlon + i
           this_dist = dist_latlon_2pts(olat(m),olon(m),rlat,rlon,PI)
           if ((this_dist /= fillvalue) .and. (this_dist < dmin)) then
             dmin = this_dist
             nmin = m
             ntcs = ntcs + 1
           end if
         end do
       end do
       if (ntcs > 0 ) then
         ! save info for later use 
         close_lat(n,k) = olat(nmin)
         close_lon(n,k) = olon(nmin)
         close_dst(n,k) = dmin
         close_ind(n,k) = nmin
       else if (dmin == fillvalue) then
         close_lat(n,k) = rmissing
         close_lon(n,k) = rmissing
         close_dst(n,k) = rmissing
         close_ind(n,k) = imissing
       else if (abs(rlon) > maxval(abs(reclon)) .or.abs(rlat) > maxval(abs(reclat)) ) then
         close_lat(n,k) = rmissing
         close_lon(n,k) = rmissing
         close_dst(n,k) = rmissing
         close_ind(n,k) = imissing      
       else
         write(iulog,*) "latlon2se_interp_init: Can't find enclosing quadrilatersl."
         write(iulog,*) "latlon2se_interp_init: lat,lon,dmin= ",olat(nmin),olon(nmin),dmin
         call endrun ('latlon2se_interp_init: error in finding nearest grid')
       endif
     end do 

     !loop corner to calculate interpolation weights
     if ((close_dst(n,5) > maxdist) .or. any(close_dst(n,:) == rmissing)) then 
       weight_qdl(n,:) = 0._r8
       close_ind(n,:)  = 1 
      !write(iulog,*) "latlon2se_interp_init: close_dst(n,1)= ", close_dst(n,1)
      !write(iulog,*) "latlon2se_interp_init: close_dst(n,2)= ", close_dst(n,2)
      !write(iulog,*) "latlon2se_interp_init: close_dst(n,3)= ", close_dst(n,3)
      !write(iulog,*) "latlon2se_interp_init: close_dst(n,4)= ", close_dst(n,4)
      !write(iulog,*) "latlon2se_interp_init: close_dst(n,5)= ", close_dst(n,5)
      !call endrun ('latlon2se_interp_init: error in finding interpolation points')
     else 
       !cshift is used to shift array to the left so 
       !corner(4,4) is the nearest point to target 
       oc   = minloc(close_dst(n,1:4),mask=(close_dst(n,1:4)==close_dst(n,5)))
       iorg = oc(1)
      !write(iulog,*) 'se2latlon_interp_init: iorg,close_dst=', iorg,close_dst(m,iorg)
       sh_corn(1:4) = cshift(close_ind(n,1:4), iorg)
       sh_rdst(1:4) = cshift(close_dst(n,1:4), iorg)
       sh_rlat(1:4) = cshift(close_lat(n,1:4), iorg)
       sh_rlon(1:4) = cshift(close_lon(n,1:4), iorg)
       !calculate the interpolation weights  
       call coord_ind_weight(sh_corn,sh_rlat,sh_rlon,sh_rdst,rlat,rlon, & !In 
                             close_lat(n,iorg),close_lon(n,iorg),close_dst(n,iorg), & !In 
                             p,q) !Out

       !calculate the final weight for each corner 
       do k = 1, 5
         weight_qdl(n,k) = 0._r8
         if(k==1) weight_qdl(n,k) = (1.0_r8 - q) * p
         if(k==2) weight_qdl(n,k) = q * p
         if(k==3) weight_qdl(n,k) = q * (1.0_r8 - p)
         if(k==4) weight_qdl(n,k) = (1.0_r8 - q) * (1.0_r8 - p)
       end do
       close_ind(n,1:4) = sh_corn(1:4)
       close_dst(n,1:4) = sh_rdst(1:4)
       close_lat(n,1:4) = sh_rlat(1:4)
       close_lon(n,1:4) = sh_rlon(1:4)
       !if(masterproc) then
       !  write(iulog,*) "latlon2se_interp_init: close_dst, close_lat, close_lon, rlat, rlon, weight_qdl= "
       !  write(iulog,*) close_dst(n,1),close_dst(n,2),close_dst(n,3),close_dst(n,4),close_dst(n,5)
       !  write(iulog,*) close_lat(n,1),close_lat(n,2),close_lat(n,3),close_lat(n,4),close_lat(n,5),rlat
       !  write(iulog,*) close_lon(n,1),close_lon(n,2),close_lon(n,3),close_lon(n,4),close_lon(n,5),rlon
       !  write(iulog,*) weight_qdl(n,1),weight_qdl(n,2),weight_qdl(n,3),weight_qdl(n,4),weight_qdl(n,5) 
       !end if
     end if 
     deallocate(olat)
     deallocate(olon)
   end do 

   deallocate(close_lat)
   deallocate(close_lon)
   deallocate(close_dst)

   !final sanity check 
   if (any(close_ind == imissing)) then 
     call endrun('latlon2se_interp_init: ERROR, interpolation index contains missing values')
   end if 
   if (any( weight_qdl == rmissing)) then
     call endrun('latlon2se_interp_init: ERROR, interpolation weight contains missing values')
   end if 
   !write(iulog,*) 'latlon2se_interp_init: interpolation initialization succeed!'

   return 
  end subroutine !latlon2se_interp_init

  subroutine se2latlon_interp_init(ngcols, rlon_d, rlat_d, & 
                                   nrlon, nrlat, reclon, reclat, & 
                                   close_ind, weight_qdl) 
   !---------------------------------------------------------------------------
   ! This program computes weighting functions to map a variable on 
   ! SE grid to (nlat,nlon) resolution
   ! Author: Shixuan Zhang  -- August 2023
   !---------------------------------------------------------------------------
   use shr_kind_mod,    only : r8 => shr_kind_r8
   use shr_const_mod,   only : SHR_CONST_PI, SHR_CONST_REARTH
   use dycore,          only : dycore_is

   implicit none
   integer, parameter   :: ncorner = 4 
   real(r8), parameter  :: re      = SHR_CONST_REARTH
   real(r8), parameter  :: PI      = SHR_CONST_PI
   real(r8), parameter  :: half_PI = 0.5_r8 * SHR_CONST_PI
   real(r8), parameter  :: two_PI  = 2.0_r8 * SHR_CONST_PI
   real(r8), parameter  :: deg2rad = SHR_CONST_PI/180._r8

   integer,  intent(in) :: nrlat, nrlon, ngcols
   real(r8), intent(in) :: rlon_d(ngcols) !lon for model grid in radians 
   real(r8), intent(in) :: rlat_d(ngcols) !lat for model grid in radians  
   real(r8), intent(in) :: reclon(nrlon)
   real(r8), intent(in) :: reclat(nrlat)

   integer,  allocatable, intent(out) :: close_ind(:,:)
   real(r8), allocatable, intent(out) :: weight_qdl(:,:)

   real(r8), allocatable :: close_lat(:,:)
   real(r8), allocatable :: close_lon(:,:)
   real(r8), allocatable :: close_dst(:,:)
   real(r8), allocatable :: olon(:)
   real(r8), allocatable :: olat(:)

   ! local variables
   logical  :: neg_root  
   integer  :: i, j, ifld, n, m, k, ii, jj 
   integer  :: iorg, ntcs, nmin, nrecg
   integer  :: oc(1)
   integer  :: sh_corn(ncorner)
   real(r8) :: sh_rlat(ncorner), sh_rlon(ncorner), sh_rdst(ncorner)
   real(r8) :: rlat, rlon, dlat, dlon, rtemp
   real(r8) :: maxdist, dmin, this_dist
   real(r8) :: p, q

   if (.not. dycore_is('UNSTRUCTURED')) then
      call endrun ('se2latlon_interp_init ERROR: only works for  SE dycore')
   end if 

   !maxdistance for the search 
   maxdist = 1.2_r8*(43200.0_r8/ngcols)*deg2rad
   nrecg   = nrlat * nrlon

   !initialize arrary for interpolation 
   allocate(close_lat(nrecg,5))
   allocate(close_lon(nrecg,5))
   allocate(close_dst(nrecg,5))
   allocate(close_ind(nrecg,5))
   allocate(weight_qdl(nrecg,5))
   close_lat(:,:)  = rmissing
   close_lon(:,:)  = rmissing
   close_dst(:,:)  = rmissing
   close_ind(:,:)  = imissing 
   weight_qdl(:,:) = rmissing

   !find the cell-center grid point closest to machine learning lat-lon grid
   do j = 1, nrlat
     do i = 1, nrlon
       !use a temporary array olat/olon so that we can modify values 
       allocate(olat(ngcols))
       allocate(olon(ngcols))
       olat = rlat_d
       olon = rlon_d
       rlon = reclon(i)
       rlat = reclat(j)
       m    = (j-1)*nrlon + i
       !Find center and 4 corner values that are closed to target grid 
       !-----------------------------------------------------
       !locate the model grid closest to the lat-lon grid
       !corners --(4,4)--         (i+1,j+1)
       !corners   (i,j+1)       --(3,3)--
       !corners        --(xp,yp)-- 
       !corners --(1,1)--         (i+1,j) 
       !corners   (i,j)         --(2,2)--
       do k = 5,1,-1 ! 4 corners + 1 center location 
         if (k < 5 ) then !nearest 
           olat(:) = rmissing 
           olon(:) = rmissing
           dlat    = close_lat(m,5) - rlat 
           dlon    = close_lon(m,5) - rlon
           if(dlat <= 0._r8 .and. dlon <= 0._r8) then 
             if (k == 1) then !(i,j), bottom left side of rlat,rlon
               olat = rlat_d 
               olon = rlon_d 
             end if
             if (k == 2) then !(i+1,j), bottom right side of rlat,rlon
               olat = merge(rlat_d,olat,rlat_d<min(close_lat(m,3),rlat).and.rlon_d>max(close_lon(m,5),rlon))
               olon = merge(rlon_d,olon,rlat_d<min(close_lat(m,3),rlat).and.rlon_d>max(close_lon(m,5),rlon))
             end if
             if (k == 3) then !(i+1,j+1), top right side of rlat,rlon
               olat = merge(rlat_d,olat,rlat_d>max(close_lat(m,5),rlat).and.rlon_d>max(close_lon(m,4),rlon)) 
               olon = merge(rlon_d,olon,rlat_d>max(close_lat(m,5),rlat).and.rlon_d>max(close_lon(m,4),rlon))
             end if
             if (k == 4) then !(i,j+1), top left side of rlat, rlon
               olat = merge(rlat_d,olat,rlat_d>max(close_lat(m,5),rlat).and.rlon_d<=rlon)
               olon = merge(rlon_d,olon,rlat_d>max(close_lat(m,5),rlat).and.rlon_d<=rlon)
             end if
           else if (dlat <= 0._r8 .and. dlon > 0._r8) then
             if (k == 1) then !(i,j)
               olat = merge(rlat_d,olat,rlat_d<min(close_lat(m,4),rlat).and.rlon_d<min(close_lon(m,2),rlon))
               olon = merge(rlon_d,olon,rlat_d<min(close_lat(m,4),rlat).and.rlon_d<min(close_lon(m,2),rlon))
             end if
             if (k == 2) then !(i+1,j)
               olat = rlat_d
               olon = rlon_d
             end if
             if (k == 3) then !(i+1,j+1)
               olat = merge(rlat_d,olat,rlat_d>max(close_lat(m,5),rlat).and.rlon_d>max(close_lon(m,4),rlon))
               olon = merge(rlon_d,olon,rlat_d>max(close_lat(m,5),rlat).and.rlon_d>max(close_lon(m,4),rlon))
             end if
             if (k == 4) then !(i,j+1)
               olat = merge(rlat_d,olat,rlat_d>max(close_lat(m,5),rlat).and.rlon_d<min(close_lon(m,5),rlon))
               olon = merge(rlon_d,olon,rlat_d>max(close_lat(m,5),rlat).and.rlon_d<min(close_lon(m,5),rlon))
             end if
           else if (dlat > 0._r8 .and. dlon > 0._r8) then
             if (k == 1) then !(i,j)
               olat = merge(rlat_d,olat,rlat_d<min(close_lat(m,4),rlat).and.rlon_d<min(close_lon(m,2),rlon))
               olon = merge(rlon_d,olon,rlat_d<min(close_lat(m,4),rlat).and.rlon_d<min(close_lon(m,2),rlon))
             end if
             if (k == 2) then !(i+1,j)
               olat = merge(rlat_d,olat,rlat_d<min(close_lat(m,3),rlat).and.rlon_d>max(close_lon(m,4),rlon))
               olon = merge(rlon_d,olon,rlat_d<min(close_lat(m,3),rlat).and.rlon_d>max(close_lon(m,4),rlon))
             end if
             if (k == 3) then !(i+1,j+1)
               olat = rlat_d
               olon = rlon_d
             end if
             if (k == 4) then !(i,j+1), top left side of rlat, rlon
               olat = merge(rlat_d,olat,rlon_d<min(close_lon(m,5),rlon).and.rlat_d>=rlat)
               olon = merge(rlon_d,olon,rlon_d<min(close_lon(m,5),rlon).and.rlat_d>=rlat)
             end if
           else if (dlat > 0._r8 .and. dlon <= 0._r8) then
             if (k == 1) then !(i,j)
               olat = merge(rlat_d,olat,rlat_d<min(close_lat(m,4),rlat).and.rlon_d<min(close_lon(m,2),rlon))
               olon = merge(rlon_d,olon,rlat_d<min(close_lat(m,4),rlat).and.rlon_d<min(close_lon(m,2),rlon))
             end if
             if (k == 2) then !(i+1,j)
               olat = merge(rlat_d,olat,rlat_d<min(close_lat(m,3),rlat).and.rlon_d>max(close_lon(m,4),rlon))
               olon = merge(rlon_d,olon,rlat_d<min(close_lat(m,3),rlat).and.rlon_d>max(close_lon(m,4),rlon))
             end if
             if (k == 3) then !(i+1,j+1)
               olat = merge(rlat_d,olat,rlon_d>max(close_lon(m,4),rlon).and.rlat_d>=rlat)
               olon = merge(rlon_d,olon,rlon_d>max(close_lon(m,4),rlon).and.rlat_d>=rlat)
             end if
             if (k == 4) then !(i,j+1)
               olat = rlat_d
               olon = rlon_d
             end if
           else 
             call endrun('se2latlon_interp_init: ERROR in buidling the coordinates')
           end if 
         end if 

         !loop model grid to locate nearest point 
         ntcs  = 0  
         nmin  = 0
         dmin  = fillvalue ! initialize a large value 
         do n = 1,ngcols 
           this_dist = dist_latlon_2pts(olat(n),olon(n),rlat,rlon,PI)
           if ((this_dist /= fillvalue) .and. (this_dist < dmin)) then 
             dmin = this_dist
             nmin = n
             ntcs = ntcs + 1
           endif
         end do
         if (ntcs > 0) then
           !save info for later use 
           close_lat(m,k) = olat(nmin)
           close_lon(m,k) = olon(nmin)
           close_dst(m,k) = dmin
           close_ind(m,k) = nmin
         else 
           write(iulog,*) "se2latlon_interp_init: Can't find enclosing quadrilatersl."
           write(iulog,*) "se2latlon_interp_init: lat,lon,dmin= ",olat(nmin),olon(nmin),dmin
           call endrun ('se2latlon_interp_init: error in finding nearest grid')
         endif
       end do 

       !loop corner to calculate interpolation weights
       if ((close_dst(m,5) > maxdist) .or. any(close_dst(m,:) == rmissing)) then
         write(iulog,*) "se2latlon_interp_init: close_dst(n,1)= ", close_dst(n,1)
         write(iulog,*) "se2latlon_interp_init: close_dst(n,2)= ", close_dst(n,2)
         write(iulog,*) "se2latlon_interp_init: close_dst(n,3)= ", close_dst(n,3)
         write(iulog,*) "se2latlon_interp_init: close_dst(n,4)= ", close_dst(n,4)
         write(iulog,*) "se2latlon_interp_init: close_dst(n,5)= ", close_dst(n,5)
         call endrun ('se2latlon_interp_init: error in finding interpolation coefficients!')
       else 
         !cshift is used to shift array to the left so 
         !corner(4,4) is the nearest point to target 
         oc   = minloc(close_dst(m,1:4),mask=(close_dst(m,1:4)==close_dst(m,5)))
         iorg = oc(1)
         sh_corn(1:4) = cshift(close_ind(m,1:4), iorg)
         sh_rdst(1:4) = cshift(close_dst(m,1:4), iorg)
         sh_rlat(1:4) = cshift(close_lat(m,1:4), iorg)
         sh_rlon(1:4) = cshift(close_lon(m,1:4), iorg)
         !calculate the interpolation weights  
         call coord_ind_weight(sh_corn,sh_rlat,sh_rlon,sh_rdst,rlat,rlon, & !In 
                               close_lat(m,iorg),close_lon(m,iorg),close_dst(m,iorg), & !In 
                               p,q) !Out

         !calculate the final weight for each corner 
         do k = 1, 5 
           weight_qdl(m,k) = 0._r8
           if(k==1) weight_qdl(m,k) = (1.0_r8 - q) * p 
           if(k==2) weight_qdl(m,k) = q * p 
           if(k==3) weight_qdl(m,k) = q * (1.0_r8 - p) 
           if(k==4) weight_qdl(m,k) = (1.0_r8 - q) * (1.0_r8 - p) 
         end do 
         close_ind(m,1:4) = sh_corn(1:4)
         close_dst(m,1:4) = sh_rdst(1:4)
         close_lat(m,1:4) = sh_rlat(1:4)
         close_lon(m,1:4) = sh_rlon(1:4)
        !write(iulog,*) "se2latlon_interp_init: close_dst,close_lat,close_lon,rlat,rlon,weight_qdl= "
        !write(iulog,*) close_dst(m,1),close_dst(m,2),close_dst(m,3),close_dst(m,4),close_dst(m,5)
        !write(iulog,*) close_lat(m,1),close_lat(m,2),close_lat(m,3),close_lat(m,4),close_lat(m,5),rlat
        !write(iulog,*) close_lon(m,1),close_lon(m,2),close_lon(m,3),close_lon(m,4),close_lon(m,5),rlon
        !write(iulog,*) weight_qdl(m,1),weight_qdl(m,2),weight_qdl(m,3),weight_qdl(m,4),weight_qdl(m,5)
       end if 
       deallocate(olat)
       deallocate(olon)
     end do
   end do 

   deallocate(close_lat)
   deallocate(close_lon)
   deallocate(close_dst)
    
   !final sanity check 
   if (any(close_ind == imissing)) then 
     call endrun('se2latlon_interp_init: ERROR, interpolation index contains missing values')
   end if 
   if (any( weight_qdl == rmissing)) then 
     call endrun('se2latlon_interp_init: ERROR, interpolation weight contains missing values')
   end if 
   !write(iulog,*) 'se2latlon_interp_init: interpolation initialization succeed!' 

   return 
  end subroutine !se2latlon_interp_init
  
  subroutine coord_ind_weight(sh_corn,sh_rlat,sh_rlon,sh_rdst,rlat,rlon, & 
                              close_lat,close_lon,close_dst,p,q)
   !---------------------------------------------------------------------------
   ! This program computes weighting functions to map a variable on 
   ! SE grid to (nlat,nlon) resolution
   ! Author: Shixuan Zhang  -- August 2023
   !---------------------------------------------------------------------------
   use shr_kind_mod,    only : r8 => shr_kind_r8
   use shr_const_mod,   only : SHR_CONST_PI, SHR_CONST_REARTH
   use dycore,          only : dycore_is
   implicit none
   integer, parameter    :: ncorner = 4
   real(r8), parameter   :: re      = SHR_CONST_REARTH
   real(r8), parameter   :: PI      = SHR_CONST_PI
   real(r8), parameter   :: half_PI = 0.5_r8 * SHR_CONST_PI
   real(r8), parameter   :: two_PI  = 2.0_r8 * SHR_CONST_PI
   real(r8), parameter   :: deg2rad = SHR_CONST_PI/180._r8

   integer,  intent(in)  :: sh_corn(ncorner)
   real(r8), intent(in)  :: sh_rlat(ncorner)
   real(r8), intent(in)  :: sh_rlon(ncorner)
   real(r8), intent(in)  :: sh_rdst(ncorner)
   real(r8), intent(in)  :: rlat
   real(r8), intent(in)  :: rlon
   real(r8), intent(in)  :: close_lat 
   real(r8), intent(in)  :: close_lon 
   real(r8), intent(in)  :: close_dst
   real(r8), intent(out) :: p, q 

   ! local variables
   logical  :: neg_root
   integer  :: i, j, ifld, n, m, k, ii, jj
   integer  :: iorg, ntcs, nmin, nrecg
   integer  :: oc(1)
   real(r8) :: dlat, dlon, rtemp, this_dist
   real(r8) :: xaa, xbb, xcc, as, bs, cs, disc
   real(r8) :: xaxbr, det, p1, p2, p_neg, q_neg
   real(r8) :: lon1c, lon2c, dist, angle, bearo, x_o, y_o
   real(r8) :: xa(3), bearg(3), xplr(3), yplr(3), xb(2)

   xa    = rmissing
   xb    = rmissing
   xaxbr = rmissing
   xplr  = rmissing
   yplr  = rmissing
   bearg = rmissing
   do ii = 3,1,-1  ! clockwise loop from 1 to 1
     if (half_PI - abs(sh_rlat(4)) < epsilon(sh_rlat(4))) then
       lon1c = 0.0_r8
     else
       lon1c = sh_rlon(4)
     endif
     if (half_PI - abs(sh_rlat(ii)) < epsilon(sh_rlat(ii))) then
       lon2c = 0.0_r8
     else
       lon2c = sh_rlon(ii)
     endif
     dlon = lon2c - lon1c
     dlon = mod(dlon,PI) - PI*int(dlon/PI)
     bearg(ii) = atan2(cos(sh_rlat(ii))*sin(dlon),  &
                       cos(sh_rlat(4))*sin(sh_rlat(ii)) &
                         - sin(sh_rlat(4))*cos(sh_rlat(ii))*cos(dlon))

     this_dist = dist_latlon_2pts(sh_rlat(4),sh_rlon(4),sh_rlat(ii),sh_rlon(ii),PI)
     if (this_dist /= fillvalue) then 
       angle    = bearg(ii) - bearg(3) !bearg(3) - bearg(ii)
       angle    = mod(angle,PI) - PI*int(angle/PI)
       xplr(ii) = this_dist * cos(angle)
       yplr(ii) = this_dist * sin(angle)
     else 
       write(iulog,*) "qdl_interp_init: corner lat/lon = ", sh_rlat(4), sh_rlon(4), sh_rlat(ii), sh_rlon(ii) 
       call endrun('qdl_interp_init: ERROR in computing distance between two corner points')
     end if 
   end do
   xaxbr = bearg(3)
   xa(3) = xplr(3)
   xa(2) = xplr(1)
   xa(1) = xplr(2) - xplr(1) - xplr(3)
   xb(2) = yplr(1)
   xb(1) = yplr(2) - yplr(1)

   !target grid location
   if ((half_PI - abs(close_lat)) < epsilon(close_lat)) then
     lon1c = 0.0_r8
   else
     lon1c = close_lon
   end if
   if ((half_PI - abs(rlat)) < epsilon(rlat)) then
     lon2c = 0.0_r8
   else
     lon2c = rlon
   end if
   dlon  = lon2c - lon1c
   dlon  = mod(dlon,PI) - PI*int(dlon/PI)
   bearo = atan2(cos(rlat)*sin(dlon),  &
                 cos(close_lat)*sin(rlat) &
                  - sin(close_lat)*cos(rlat)*cos(dlon))
   angle = bearo - xaxbr !xaxbr - bearo
   angle = mod(angle,PI) - PI*int(angle/PI)
   x_o   = close_dst * cos(angle)
   y_o   = close_dst * sin(angle)

   !solving the quadratic equation to obtain interpolation weight
   p     = rmissing
   q     = rmissing
   p1    = rmissing
   p2    = rmissing
   p_neg = rmissing
   q_neg = rmissing

   xaa   = xa(1) * xb(2) - xa(2) * xb(1)
   xbb   = xa(3) * xb(2) - xa(1) * y_o + xb(1) * x_o
   xcc   = -xa(3) * y_o
   as    = xaa / max(abs(xaa), abs(xbb), abs(xcc))
   bs    = xbb / max(abs(xaa), abs(xbb), abs(xcc))
   cs    = xcc / max(abs(xaa), abs(xbb), abs(xcc))

   if (abs(as) < epsilon(as)) then
     p1 = - cs / bs
   else
     disc = bs * bs - 4.0_r8 * as * cs
     if (disc >= 0.0_r8) then
       if(bs >= 0.0_r8) then
         p1 = (-bs - sqrt(disc)) / (2.0_r8 * as)
       else
         p1 = (-bs + sqrt(disc)) / (2.0_r8 * as)
       end if
       if (p1 == 0.0_r8) then
         p2 = 0.0_r8
       else
         p2 = cs / (as * p1)
       end if
     end if
   end if
   if (p1 == rmissing .and. p2 == rmissing) then
     call endrun('qdl_interp_init: ERROR in find interp p1, p2 coef.')
   end if

   !calculate the weight based on solution from quadratic equation
   if (xaa > 0.0_r8) then
     if (p1 >=0.0_r8 .and. p1 <= 1._r8) then
       p = p1
     else if (p2 >=0.0_r8 .and. p2 <= 1._r8) then
       p = p2
     else
       p = rmissing
     end if
   else if (p1 /= rmissing .and. p2 == rmissing ) then
     p = p1
   else
     p = p1
     if (xbb > 0.0_r8) then
       if (p > 1.0_r8) then
         p = p2
       else
         p_neg = p2
       end if
     else
       write(iulog,*) 'qdl_interp_init: xaa < 0 and xbb <= 0: no mapping is possible'
       write(iulog,*) 'qdl_interp_init: xaa,xbb,xcc,x_o,y_o,angle = ', xaa, xbb, xcc, x_o, y_o, angle
       call endrun('qdl_interp_init: ERROR in find interp p coef.')
     end if
   end if

   if (p < 0.0_r8 .or. p > 1.0_r8) then
     write(iulog,*) 'qdl_interp_init: target grid is out of bounds p,p1,p2 = [0,1] ',p,p1,p2
     call endrun('qdl_interp_init: ERROR in find interp p coef.')
   end if

   ! Use p to calculate the other unit square coordinate value, 'q'.
   det = xa(3) + xa(1) * p
   if (det /= 0.0_r8) then
     q = (x_o - xa(2)*p) / det
   else
     write(iulog,*) 'qdl_interp_init: xa(3), xa(1), p = ', xa(3), xa(1), p
     call endrun('qdl_interp_init: ERROR in find interp q coef.')
   end if

   ! Repeat for the -root, if it is a possibility.
   if (p_neg /= rmissing) then
     det = (xa(3) + xa(1)*p_neg)
     if (det /= 0.0_r8) then
       q_neg = (x_o - xa(2)*p_neg) / det
     else
       write(iulog,*) 'qdl_interp_init: xa(3), xa(1), p_neg = ', xa(3), xa(1), p_neg
       call endrun('qdl_interp_init: ERROR in find interp q_neg coef.')
     end if
   end if
   if (q < 0.0_r8 .or. q > 1.0_r8) then
     write(iulog,*) 'qdl_interp_init: out of range [0,1], q = ', q
     call endrun('qdl_interp_init: ERROR in find interp q coef.')
   end if

   ! select the right values in l and r.
   neg_root = p_neg >= 0.0_r8 .and. p_neg <= 1.0_r8 .and. &
              q_neg >= 0.0_r8 .and. q_neg <= 1.0_r8

   if (p >= 0.0_r8 .and. p <= 1.0_r8 .and. &
       q >= 0.0_r8 .and. q <= 1.0_r8 ) then

     ! Both roots yield a good mapping.
     if (neg_root) then
       write(iulog,*) 'qdl_interp_init: BOTH roots of the m quadratic &
                       yield usable mappings.  The +root is being used.'
     endif

   else if (neg_root) then

     ! The -root yields a good mapping.  Pass along the -root m and l.
     p = p_neg
     q = q_neg
     write(iulog,*) 'qdl_interp_init: The negative root of the m quadratic &
                     yields the only usable mapping.'

   end if

   !Diag output
   !write(iulog,*) 'qdl_interp_init: normalized coef as,bs,cs=', as,bs,cs
   !write(iulog,*) 'qdl_interp_init: xa(3)           =', xa(1),xa(2),xa(3) 
   !write(iulog,*) 'qdl_interp_init: xb(2)           =', xb(1),xb(2) 
   !write(iulog,*) 'qdl_interp_init: xaxbr           =', xaxbr 
   !write(iulog,*) 'qdl_interp_init: xaxbr           =', xaxbr 
   !write(iulog,*) 'qdl_interp_init: sh_ind          =', sh_corn(1),sh_corn(2),sh_corn(3),sh_corn(4) 
   !write(iulog,*) 'qdl_interp_init: sh_lat          =', sh_rlat(1),sh_rlat(2),sh_rlat(3),sh_rlat(4)
   !write(iulog,*) 'qdl_interp_init: sh_lon          =', sh_rlon(1),sh_rlon(2),sh_rlon(3),sh_rlon(4)
   !write(iulog,*) 'qdl_interp_init: sh_dst          =', sh_rdst(1),sh_rdst(2),sh_rdst(3),sh_rdst(4)
   !write(iulog,*) 'qdl_interp_init: close_lat       =', close_lat
   !write(iulog,*) 'qdl_interp_init: close_lon       =', close_lon
   !write(iulog,*) 'qdl_interp_init: close_dst       =', close_dst
   !write(iulog,*) 'qdl_interp_init: p,q,p_neg,q_neg =', p,q,p_neg,q_neg

   return 
  end subroutine !coord_ind_weight

end module nudging
