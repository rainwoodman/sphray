!> \file global.F90

!> \brief the module that handles global variables
!<
module global_mod
use myf03_mod 
use gadget_sphray_header_class, only: gadget_sphray_header_type
use particle_system_mod, only: particle_system_type
use oct_tree_mod, only: oct_tree_type
use atomic_rates_mod, only: atomic_rates_table_type
use atomic_rates_mod, only: atomic_rates_type
use ray_mod, only: src_ray_type
implicit none


!> snapshot information type. 
!==============================
type snap_info_type
   real(r8b)    :: TimeAt            !< t at snapshot
   real(r8b)    :: ScalefacAt        !< scale factor at snapshot
   real(r8b)    :: TimeToNext        !< TimeAt(n+1) - TimeAt(n)
   real(r8b)    :: RunTime           !< duration to trace snapshot
   real(r8b)    :: StartTime         !< t0 for tracing snapshot 
   integer(i8b) :: RaysFromSrcHeader !< rays listed in source header
   integer(i8b) :: SrcRays           !< rays to trace for snapshot
end type snap_info_type


!> run planning type.
!======================
type run_planning_type  
   type(snap_info_type), allocatable :: snap(:) !< snap information
   real(r8b), allocatable :: OutputTimes(:)     !< times to make outputs 
end type run_planning_type


! global variables
!=====================
type(particle_system_type) :: psys        !< particles + sources + box
type(oct_tree_type) :: tree               !< octree
type(src_ray_type), allocatable :: active_rays(:)   !< rays currently tracing
type(atomic_rates_table_type) :: rtable   !< rates read in from file
type(atomic_rates_type) :: isoT_k         !< static rates for iso-temperature run
type(atomic_rates_type) :: cmbT_k         !< static rates for cmb-temperature
type(atomic_rates_type) :: xHII_k         !< static rates for xHII-temperature 

type(run_planning_type) :: PLAN           !< run plan

type(gadget_sphray_header_type), allocatable :: saved_gheads(:,:) !< all headers (nsnaps,nfiles)


 
!> global variables type. 
!=========================
type global_variables_type


   ! particle and source sizes
   !-----------------------------------------------------------
   integer(i4b) :: bytesperpar  !< bytes per particle
   integer(i4b) :: bytespersrc  !< bytes per source
   real(r8b) :: MB              !< tracks memory consumption


   ! these are all set in get_planning_data in main_input.f90
   !-----------------------------------------------------------
   integer(i4b) :: Nsnaps !< snap count (particle=source)

   real(r8b) :: cgs_len   !< code length [cm/h]
   real(r8b) :: cgs_mass  !< code mass [g/h]
   real(r8b) :: cgs_vel   !< code velocity [cm/s] 

   real(r8b) :: cgs_time  !< code time [s/h]
   real(r8b) :: cgs_rho   !< code density [g/cm^3 h^2]
   real(r8b) :: cgs_prs   !< code pressure [dyne/cm^2 h^2]
   real(r8b) :: cgs_enrg  !< code energy [ergs/h]
   real(r8b) :: Lunit     !< code source luminosity unit [photons/s]

   real(r8b) :: LenFac_cm  !< code (input) length -> cm = cgs_len * a / h
   real(r8b) :: MassFac_g  !< code (input) mass -> g = cgs_mass / h
   real(r8b) :: TimeFac_s  !< code (input) time -> s = cgs_time / h

   real(r8b) :: OmegaM    !< matter / critical density z=0
   real(r8b) :: OmegaB    !< baryon / critical density z=0
   real(r8b) :: OmegaL    !< lambda / critical density z=0
   real(r8b) :: LittleH   !< Hubble parameter z=0 in units of 100 km/s/Mpc

   ! these are set in do_output_planning and do_ray_planning in intialize.f90
   !--------------------------------------------------------------------------
   integer(i4b) :: OutputIndx                     !< keeps track of outputs   
   real(r8b)    :: TotalSimTime                   !< total time to ray trace
   integer(i4b) :: NumTotOuts                     !< total outputs to do
 

   ! these should be reset each time a new snapshot is read in
   !------------------------------------------------------------
   integer(i8b) :: itime           !< integer time. measures ticks from start time

   real(r8b) :: start_time_code    !< starting time in code units
   real(r8b) :: start_time_s       !< starting time in seconds
   real(r8b) :: start_time_myr     !< starting time in Myr

   real(r8b) :: time_elapsed_code  !< elapsed time in code units
   real(r8b) :: time_elapsed_s     !< elapsed time in seconds
   real(r8b) :: time_elapsed_myr   !< elapsed time in Myr 

   real(r8b) :: dt_code            !< one tick in code units
   real(r8b) :: dt_s               !< one tick in seconds
   real(r8b) :: dt_myr             !< one tick in Myr

   real(r8b) :: BoxLwrs(3)         !< input coords [code] of lower x,y,z corner
   real(r8b) :: BoxUprs(3)         !< input coords [code] of upper x,y,z corner

   real(r8b) :: total_mass         !< summed mass of all particles in a snapshot
   real(r8b) :: total_lum          !< summed luminosity of sources in a snapshot
   real(r8b) :: total_atoms        !< sum of all atoms in computational volume
   real(r8b) :: total_photons      !< sum of all photons to be released
   real(r8b) :: Tcmb_cur           !< CMB temperature for the current snapshot

   real(r8b) :: sf_gamma_eos       !< index for polytropic equation of state for star forming gas.

   real(r8b) :: UVB_gammaHI_cloudy !< magnitude of UVB from cloudy ionization table
   
   ! these are updated continuosly while the code runs 
   ! and most are initialized in initialize.f90
   !----------------------------------------------------
   character(clen) :: ionfrac_file !< file where mini outputs are put 
   integer(i4b) :: ionlun          !< lun for ionfrac log file

   character(clen) :: raystat_file !< file where ray stats are put
   integer(i4b) :: raystatlun      !< lun for ray stat log file

   character(clen) :: pardata_file !< file with particle data summaries
   integer(i4b) :: pardatalun      !< lun for particle data log file
   
   character(clen) :: srcdata_file !< file with source data summaries
   integer(i4b) :: srcdatalun      !< lun for source data log file

   character(clen) :: rayplan_file  !< file for ray planning
   integer(i4b) :: rayplanlun      !< lun for ray planning file

   character(clen) :: outplan_file  !< file for out planning
   integer(i4b) :: outplanlun      !< lun for out planning file

   integer(i4b) :: CurSnapNum      !< current snapshot number

   real(r8b) :: nwionfrac          !< number weighted ionization fraction
   real(r8b) :: mwionfrac          !< mass weighted ionization fraction
   real(r8b) :: vwionfrac          !< volume weighted ionization fraction

   real(r8b) :: TotalSourceRaysCast    !< total rays traced from user def. sources
   real(r8b) :: TotalDiffuseRaysCast   !< total recombination rays traced
   real(r8b) :: IonizingPhotonsPerSec  !< ionizing photons emitted per second
   real(r8b) :: TotalPhotonsCast       !< total number of photons emitted
   real(r8b) :: TotalPhotonsAbsorbed   !< total number of photons absorbed
   real(r8b) :: PhotonsLeavingBox      !< total photons leaving the box
   real(r8b) :: TotalIonizations       !< total number of photoionizations
   real(r8b) :: TotalRecombinations    !< total number of recombinations
   
   real(r8b) :: PeakUpdates            !< max updates for a single particle
   real(r8b) :: AverageUpdatesPerPar   !< average number of updates per particle
   real(r8b) :: ParticleCrossings      !< number of ray / particle intersections
   real(r8b) :: TotalDerivativeCalls   !< times the solver being used has run
   integer(i8b) :: rayoops             !< number of times a ray hits a par hit by a future ray
   integer(i8b) :: totalhits !< a total number of hits
end type global_variables_type
 

type(global_variables_type) :: GV           !< global variables          





end module global_mod
