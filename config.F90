!> \file config.F90

!> \brief the module that handles the config file
!<
module config_mod
use myf03_mod
#ifdef useHDF5
use hdf5_wrapper
#endif
implicit none
private

public :: config_variables_type
public :: CV 
public :: read_config_file
public :: write_config_hdf5_lun
public :: dummy_check_config_variables
public :: config_info_to_file


! these variables are read in directly from the config file
!--------------------------------------------------------------
type config_variables_type


   integer(i4b)    :: Verbosity           !<  0=silent, 1=whisper, 2=talk, 3=debug

   logical         :: DoTestScenario      !<  set true if performing a test problem
   character(clen) :: TestScenario        !<  one of {iliev_test1, iliev_test2, iliev_test3, iliev_test4}

   logical         :: JustInit            !<  set true to stop after initialization
   logical         :: Comoving            !<  set true if values to be read are in comoving coords

   real(r8b)       :: IsoTemp             !<  if > zero all pars fixed @ IsoTemp (FixSnapTemp must be F)
   logical         :: FixSnapTemp         !<  if T, fix temp at snapshot values (IsoTemp must be <= 0)

   real(r8b)       :: EOStemp             !<  if non-negative, initialize EOS particles w/ T = EOStemp 
   real(r8b)       :: InitxHI             !<  if non-negative, all xHI initialized to this value
   logical         :: RayDepletion        !<  remove photons from rays as they travel?

   integer(i8b)    :: IntSeed             !<  seed for mersenne twister

   real(r8b)       :: StaticFieldSimTime  !<  sim time for single snapshot jobs
   character(clen) :: StaticSimTimeUnit   !<  one of {codetime,myr}


   integer(i4b)    :: InputType           !<  one of Gadget {1: Public 2: CosmoBH 3: OWLS/GIMIC 4: V.Bromm 5: Public HDF5}
   character(clen) :: SnapPath            !<  dir where particle snapshots are
   character(clen) :: SourcePath          !<  dir where source snapshots are


   character(clen) :: SpectraFile         !<  file containing spectra tables
   character(clen) :: b2cdFile            !<  file containing b2cd tables
   character(clen) :: AtomicRatesFile     !<  file containnig atomic rates tables

   character(clen) :: ParFileBase         !<  particle snapshot file base
   character(clen) :: SourceFileBase      !<  source snapshot file base

   integer(i4b)    :: StartSnapNum        !<  snapshot to start with
   integer(i4b)    :: EndSnapNum          !<  snapshot to end with

   integer(i4b)    :: ParFilesPerSnap     !<  files per particle snapshot
   integer(i4b)    :: SourceFilesPerSnap  !<  files per source snapshot


   character(clen) :: RayScheme           !<  one of {raynum, header}
   real(r8b)       :: ForcedRayNumber     !<  number of rays to trace if RayScheme = raynum

   logical         :: RayStats            !<  T = massive output file on ray statistics in raystats.dat
   integer(i4b)    :: BndryCond           !<  one of {-1:reflecting 0:vacuum 1:periodic}

   real(r8b)       :: RayPhotonTol        !<  fractional ray depletion to stop ray
   real(r8b)       :: MaxRayDist          !<  max ray distance in physical code units, negative=default
   
   logical         :: HydrogenCaseA       !<  T = use case A for Hydrogen Recombinations
   logical         :: HeliumCaseA         !<  T = use case A for Helium Recombinsations

   integer(i4b)    :: IonTempSolver       !<  one of {1:euler, 2:bdf}

   real(r8b)       :: Tfloor              !<  minimum allowed temperature
   real(r8b)       :: Tceiling            !<  maximum allowed temperature

   real(r8b)       :: xfloor              !<  minimum allowed ionization fraction
   real(r8b)       :: xceiling            !<  maximum allowed ionization fraction

   real(r8b)       :: NeBackground        !<  constant background electron number density from metals
   integer(i8b)    :: NraysUpdateNoHits   !<  update all pars not hit by a ray in last NraysUpdateNoHits

   real(r8b)       :: H_mf                !<  hydrogen mass fraction
   real(r8b)       :: He_mf               !<  helium mass fraction

   character(clen) :: OutputDir           !<  path to output directory
   character(clen) :: OutputFileBase      !<  output file base

   integer(i4b)    :: OutputType          !<  one of {1:Standard Binary Gadget 2:HDF5 Gadget}

   character(clen) :: OutputTiming        !<  one of {standard, forced}
   integer(i4b)    :: NumStdOuts          !<  if OutputTiming = "standard", # of outputs (maybe +1 initial)

   logical         :: DoInitialOutput     !<  produces output before any raytracing
   integer(i8b)    :: IonFracOutRays      !<  do mini output every IonFracOutRays src rays

   character(clen) :: ForcedOutFile       !<  file with forced output times
   character(clen) :: ForcedUnits         !<  one of {codetime, myr, mwionfrac, vwionfrac}

   integer(i4b)    :: PartPerCell         !<  minimum particles in a tree leaf

   character(clen) :: config_file         !< name of the config file
end type config_variables_type


type(config_variables_type) :: CV         !< config variables


contains


!> writes the config file 
!==============================
subroutine write_config_hdf5_lun(lun)
  integer :: lun

#ifdef useHDF5
  call hdf5_write_attribute( lun, 'Config/Verbosity',          CV%Verbosity )

  if (CV%DoTestScenario) then
     call hdf5_write_attribute( lun, 'Config/DoTestScenario',  1)           
  else
     call hdf5_write_attribute( lun, 'Config/DoTestScenario',  0)           
  endif
  
  call hdf5_write_attribute( lun, 'Config/TestScenario',       CV%TestScenario )
  
  if (CV%JustInit) then
     call hdf5_write_attribute( lun, 'Config/JustInit',        1)           
  else
     call hdf5_write_attribute( lun, 'Config/JustInit',        0)           
  endif
  
  if (CV%Comoving) then
     call hdf5_write_attribute( lun, 'Config/Comoving',        1)         
  else
     call hdf5_write_attribute( lun, 'Config/Comoving',        0)         
  endif
  
  call hdf5_write_attribute( lun, 'Config/IsoTemp',            CV%IsoTemp)            
  
  if (CV%FixSnapTemp) then
     call hdf5_write_attribute( lun, 'Config/FixSnapTemp',     1)         
  else
     call hdf5_write_attribute( lun, 'Config/FixSnapTemp',     0)         
  endif
  
  call hdf5_write_attribute( lun, 'Config/EOStemp',            CV%EOStemp)            
  
  call hdf5_write_attribute( lun, 'Config/InitxHI',            CV%InitxHI)            
  
  if (CV%RayDepletion) then
     call hdf5_write_attribute( lun, 'Config/RayDepletion',    1)         
  else
     call hdf5_write_attribute( lun, 'Config/RayDepletion',    0)         
  endif
  
  call hdf5_write_attribute( lun, 'Config/IntSeed',            CV%IntSeed)            

  call hdf5_write_attribute( lun, 'Config/StaticFieldSimTime', &
       CV%StaticFieldSimTime) 
  call hdf5_write_attribute( lun, 'Config/StaticSimTimeUnit',  trim(CV%StaticFieldSimTimeUnit)) 
  
  call hdf5_write_attribute( lun, 'Config/InputType',          CV%InputType)          
  call hdf5_write_attribute( lun, 'Config/SnapPath',           CV%SnapPath)           
  call hdf5_write_attribute( lun, 'Config/SourcePath',         CV%SourcePath)           
  
  call hdf5_write_attribute( lun, 'Config/SpectraFile',        CV%SpectraFile)           
  call hdf5_write_attribute( lun, 'Config/b2cdFile',           CV%b2cdFile)           
  call hdf5_write_attribute( lun, 'Config/AtomicRatesFile',    CV%AtomicRatesFile)    
  
  call hdf5_write_attribute( lun, 'Config/ParFileBase',        CV%ParFileBase)        
  call hdf5_write_attribute( lun, 'Config/SourceFileBase',     CV%SourceFileBase)        
  
  call hdf5_write_attribute( lun, 'Config/StartSnapNum',       CV%StartSnapNum)            
  call hdf5_write_attribute( lun, 'Config/EndSnapNum',         CV%EndSnapNum)            
  
  call hdf5_write_attribute( lun, 'Config/ParFilesPerSnap',    CV%ParFilesPerSnap)    
  call hdf5_write_attribute( lun, 'Config/SourceFilesPerSnap', CV%SourceFilesPerSnap)    
  
  call hdf5_write_attribute( lun, 'Config/RayScheme',          CV%RayScheme)
  call hdf5_write_attribute( lun, 'Config/ForcedRayNumber',    CV%ForcedRayNumber)
  
  if (CV%RayStats) then
     call hdf5_write_attribute( lun, 'Config/RayStats',        1)
  else
     call hdf5_write_attribute( lun, 'Config/RayStats',        0)
  endif
  
  call hdf5_write_attribute( lun, 'Config/BndryCond',          CV%BndryCond)          
  call hdf5_write_attribute( lun, 'Config/RayPhotonTol',       CV%RayPhotonTol)         
  
  if (CV%HydrogenCaseA) then
     call hdf5_write_attribute( lun, 'Config/HydrogenCaseA',   1)
  else
     call hdf5_write_attribute( lun, 'Config/HydrogenCaseA',   0)
  endif
  
  if (CV%HeliumCaseA) then
     call hdf5_write_attribute( lun, 'Config/HeliumCaseA',     1)
  else
     call hdf5_write_attribute( lun, 'Config/HeliumCaseA',     0)
  endif
  
  call hdf5_write_attribute( lun, 'Config/IonTempSolver',      CV%IonTempSolver)         
  
  call hdf5_write_attribute( lun, 'Config/Tfloor',             CV%Tfloor)             
  call hdf5_write_attribute( lun, 'Config/Tceiling',           CV%Tceiling)           
  
  call hdf5_write_attribute( lun, 'Config/xfloor',             CV%xfloor)             
  call hdf5_write_attribute( lun, 'Config/xceiling',           CV%xceiling)           
  
  call hdf5_write_attribute( lun, 'Config/NeBackground',       CV%NeBackground)       
  
  call hdf5_write_attribute( lun, 'Config/NraysUpdateNoHits',  CV%NraysUpdateNoHits)                 
  
  call hdf5_write_attribute( lun, 'Config/H_mf',               CV%H_mf)               
  call hdf5_write_attribute( lun, 'Config/He_mf',              CV%He_mf)              
  
  call hdf5_write_attribute( lun, 'Config/OutputDir',          CV%OutputDir)          
  call hdf5_write_attribute( lun, 'Config/OutputFileBase',     CV%OutputFileBase)     
  call hdf5_write_attribute( lun, 'Config/OutputType',         CV%OutputType)         
  
  call hdf5_write_attribute( lun, 'Config/OutputTiming',       CV%OutputTiming)         
  call hdf5_write_attribute( lun, 'Config/NumStdOuts',         CV%NumStdOuts)         
  
  if (CV%DoInitialOutput) then
     call hdf5_write_attribute( lun, 'Config/DoInitialOutput', 1)         
  else
     call hdf5_write_attribute( lun, 'Config/DoInitialOutput', 0)         
  endif
  
  call hdf5_write_attribute( lun, 'Config/IonFracOutRays',     CV%IonFracOutRays)         
  call hdf5_write_attribute( lun, 'Config/ForcedOutFile',      CV%ForcedOutFile)         
  call hdf5_write_attribute( lun, 'Config/ForcedUnits',        CV%ForcedUnits)         
  
  call hdf5_write_attribute( lun, 'Config/PartPerCell',        CV%PartPerCell)        
#endif

  
end subroutine write_config_hdf5_lun


!> reads the config file 
!==============================
subroutine read_config_file(config_file)

  character(clen), intent(in) :: config_file  !< file to read config vars from
  character(clen) :: keyword
  character(clen) :: str
  logical :: file_exists
  integer :: verb

  character(clen), parameter :: myname="read_config_file"
  logical, parameter :: crash = .true.

    write(str,'(A,A)') 'using configuration file: ', trim(config_file)
    verb = 0
    call mywrite(str,verb) 

    inquire( file=config_file, exist=file_exists )
    if (.not. file_exists) then
       call myerr("config file does not exist",myname,crash)
    end if


    keyword = "Verbosity:"
    call scanfile(config_file,keyword,CV%Verbosity)
    myf03_verbosity = CV%verbosity

    keyword = "DoTestScenario:"
    call scanfile(config_file,keyword,CV%DoTestScenario)

    keyword = "TestScenario:"
    call scanfile(config_file,keyword,CV%TestScenario)
!-----------------------

    keyword = "JustInit:"
    call scanfile(config_file,keyword,CV%JustInit)

    keyword = "Comoving:"
    call scanfile(config_file,keyword,CV%Comoving)

    keyword = "IsoTemp:"
    call scanfile(config_file,keyword,CV%IsoTemp)

    keyword = "FixSnapTemp:"
    call scanfile(config_file,keyword,CV%FixSnapTemp)

    keyword = "EOStemp:"
    call scanfile(config_file,keyword,CV%EOStemp)

    keyword = "InitxHI:"
    call scanfile(config_file,keyword,CV%InitxHI)

    keyword = "RayDepletion:"
    call scanfile(config_file,keyword,CV%RayDepletion)

    keyword = "IntSeed:"
    call scanfile(config_file,keyword,CV%IntSeed)

    keyword = "StaticFieldSimTime:"
    call scanfile(config_file,keyword,CV%StaticFieldSimTime)

    keyword = "StaticSimTimeUnit:"
    call scanfile(config_file,keyword,CV%StaticSimTimeUnit)    

    !   input snapshot information
    !------------------------------
    keyword = "InputType:"
    call scanfile(config_file,keyword,CV%InputType)

    keyword = "SnapPath:"
    call scanfile(config_file,keyword,CV%SnapPath)

    keyword = "SourcePath:"
    call scanfile(config_file,keyword,CV%SourcePath)

    keyword = "SpectraFile:"
    call scanfile(config_file,keyword,CV%SpectraFile)

    keyword = "b2cdFile:"
    call scanfile(config_file,keyword,CV%b2cdFile)

    keyword = "AtomicRatesFile:"
    call scanfile(config_file,keyword,CV%AtomicRatesFile)

    keyword = "ParFileBase:"
    call scanfile(config_file,keyword,CV%ParFileBase)

    keyword = "SourceFileBase:"
    call scanfile(config_file,keyword,CV%SourceFileBase)

    keyword = "StartSnapNum:"
    call scanfile(config_file,keyword,CV%StartSnapNum)

    keyword = "EndSnapNum:"
    call scanfile(config_file,keyword,CV%EndSnapNum)

    keyword = "ParFilesPerSnap:"
    call scanfile(config_file,keyword,CV%ParFilesPerSnap)

    keyword = "SourceFilesPerSnap:"
    call scanfile(config_file,keyword,CV%SourceFilesPerSnap)


    !   ray tracing
    !----------------------------
    keyword = "RayScheme:"
    call scanfile(config_file,keyword,CV%RayScheme)

    keyword = "ForcedRayNumber:"
    call scanfile(config_file,keyword,CV%ForcedRayNumber)
 
    keyword = "RayStats:"
    call scanfile(config_file,keyword,CV%RayStats)

    keyword = "BndryCond:"
    call scanfile(config_file,keyword,CV%BndryCond)

    keyword = "RayPhotonTol:"
    call scanfile(config_file,keyword,CV%RayPhotonTol)


    !   ion/temp solving
    !----------------------------
    keyword = "HydrogenCaseA:"
    call scanfile(config_file,keyword,CV%HydrogenCaseA)

    keyword = "HeliumCaseA:"
    call scanfile(config_file,keyword,CV%HeliumCaseA)

    keyword = "IonTempSolver:"
    call scanfile(config_file,keyword,CV%IonTempSolver)

    keyword = "Tfloor:"
    call scanfile(config_file,keyword,CV%Tfloor)

    keyword = "Tceiling:"
    call scanfile(config_file,keyword,CV%Tceiling)

    keyword = "xfloor:"
    call scanfile(config_file,keyword,CV%xfloor)

    keyword = "xceiling:"
    call scanfile(config_file,keyword,CV%xceiling)

    keyword = "NeBackground:"
    call scanfile(config_file,keyword,CV%NeBackground)

    keyword = "NraysUpdateNoHits:"
    call scanfile(config_file,keyword,CV%NraysUpdateNoHits)

    keyword = "H_mf:"
    call scanfile(config_file,keyword,CV%H_mf)

    keyword = "He_mf:"
    call scanfile(config_file,keyword,CV%He_mf)

    !   output
    !-------------
    keyword = "OutputDir:"
    call scanfile(config_file,keyword,CV%OutputDir)

    keyword = "OutputFileBase:"
    call scanfile(config_file,keyword,CV%OutputFileBase)

    keyword = "OutputType:"
    call scanfile(config_file,keyword,CV%OutputType)

    keyword = "OutputTiming:"
    call scanfile(config_file,keyword,CV%OutputTiming)

    keyword = "NumStdOuts:"
    call scanfile(config_file,keyword,CV%NumStdOuts)

    keyword = "DoInitialOutput:"
    call scanfile(config_file,keyword,CV%DoInitialOutput)

    keyword = "IonFracOutRays:"
    call scanfile(config_file,keyword,CV%IonFracOutRays)

    keyword = "ForcedOutFile:"
    call scanfile(config_file,keyword,CV%ForcedOutFile)

    keyword = "ForcedUnits:"
    call scanfile(config_file,keyword,CV%ForcedUnits)
 
!--------------------
    keyword = "PartPerCell:"
    call scanfile(config_file,keyword,CV%PartPerCell)


    CV%config_file = config_file

    call dummy_check_config_variables()
    call config_info_to_file()

end subroutine read_config_file


!> run dummy checks on config variables
!========================================
subroutine dummy_check_config_variables()

  character(clen) :: Cwarning, Mwarning
  logical :: config_good
  logical :: charmatch

  Cwarning = "please edit " // trim(CV%config_file)
  Mwarning = "please edit Makefile"
  config_good = .true. 

  if ( CV%InputType /= 1 .and. CV%InputType /= 2 .and. &
       CV%InputType /= 3 .and. CV%InputType /= 4) then
     write(*,*) "Input Type ", CV%InputType, " not recognized"
     write(*,*) "must be 1 (Gadget2 Public Standard), &
          2 (Gadget CosmoBH), 3 (Gadget OWLS/GIMIC HDF5), &
          4 (Gadget V. Bromm), or 5 (Gadget Public HDF5)"
     config_good = .false. 
  end if

  if (CV%OutputType /= 1 .and. CV%OutputType /= 2) then
     write(*,*) "Output Type ", CV%OutputType, " not recognized"
     write(*,*) "must be 1 (Standard Binary Gadget) or 2 (HDF5)"
     config_good = .false. 
  end if


  if (CV%Tfloor < 0.0 .or. CV%Tceiling < 0.0) then
     write(*,*) "Tfloor and Tceiling must be greater than or equal to 0.0"
     config_good = .false. 
  end if

  if (CV%Tfloor > 1.0e9 .or. CV%Tceiling > 1.0e9) then
     write(*,*) "Tfloor and Tceiling must be less than or equal to 1.0e9"
     config_good = .false. 
  end if

  if (CV%Tfloor > CV%Tceiling) then
     write(*,*) "Tceiling must be greater than Tfloor"
     config_good = .false. 
  end if

  if (CV%IsoTemp > 0.0) then
  
     if (CV%IsoTemp < CV%Tfloor) then
        write(*,*) "IsoTemp cannot be set lower than Tfloor"
        config_good = .false. 
     end if
     
     if (CV%IsoTemp > CV%Tceiling) then
        write(*,*) "IsoTemp cannot be set higher than Tceiling"
        config_good = .false. 
     end if

     if (CV%FixSnapTemp) then 
        write(*,*) "Only one of the following can be true:"
        write(*,*) "FixSnapTemp = T (then IsoTemp must be < 0)"
        write(*,*) "IsoTemp > 0     (then FixSnapTemp must = F)"
        config_good = .false.
     endif

  end if

  if (CV%StartSnapNum < 0 .or. CV%EndSnapNum < 0) then
     write(*,*) "Starting and Ending snapshot numbers must be > 0"
     config_good = .false. 
  end if

  if (CV%StartSnapNum > CV%EndSnapNum) then
     write(*,*) "Starting snapshot number cannot be > Ending snapshot number"
     config_good = .false. 
  end if


  if ( CV%IonTempSolver /= 1 .and. CV%IonTempSolver /= 2) then
     config_good = .false.
     write(*,*) "IonTempSolver: ", CV%IonTempSolver, " not recognized"
     write(*,*) "must be '1'=euler or '2'=bdf "
  end if


  charmatch = .false.
  if ( trim(CV%StaticSimTimeUnit) == "codetime" ) charmatch = .true.
  if ( trim(CV%StaticSimTimeUnit) == "myr"      ) charmatch = .true.
  if (.not. charmatch) then
     write(*,*) "StaticSimTimeUnit: ", trim(CV%StaticSimTimeUnit), " not recognized"
     write(*,*) "must be 'codetime' or 'myr' "
     config_good = .false.
  end if

  charmatch = .false.
  if ( trim(CV%RayScheme) == "raynum" ) charmatch = .true.
  if ( trim(CV%RayScheme) == "header"   ) charmatch = .true.
  if (.not. charmatch) then
     write(*,*) "RayScheme: ", trim(CV%RayScheme), " not recognized"
     write(*,*) "must be 'raynum' or 'header' "
     config_good = .false.
  end if

  charmatch = .false.
  if ( trim(CV%OutputTiming) == "standard" ) charmatch = .true.
  if ( trim(CV%OutputTiming) == "forced"   ) charmatch = .true.
  if (.not. charmatch) then
     write(*,*) "OutputTiming: ", trim(CV%OutputTiming), " not recognized"
     write(*,*) "must be 'standard' or 'forced' "
     config_good = .false.
  end if

  charmatch = .false.
  if ( trim(CV%ForcedUnits) == "codetime" ) charmatch = .true.
  if ( trim(CV%ForcedUnits) == "myr"      ) charmatch = .true.
  if ( trim(CV%ForcedUnits) == "mwionfrac") charmatch = .true.
  if (.not. charmatch) then
     write(*,*) "ForcedUnits: ", trim(CV%ForcedUnits), " not recognized"
     write(*,*) "must be 'codetime', 'myr', or 'mwionfrac' "
     config_good = .false.
  end if

#ifdef incHe
  if (CV%He_mf == 0.0) then
     write(*,*) "You have defined the incHe macro in the Makefile, but"
     write(*,*) "the Helium mass fraction is set to 0.0 in the config file."
     write(*,*) "Please comment out incHe in the Makefile or set the Helium"
     write(*,*) "mass fraction to something greater than zero in the config file."
     config_good = .false.
  end if
#else
  if (CV%He_mf /= 0.0) then
     write(*,*) "You have not defined the incHe macro in the Makefile, but"
     write(*,*) "the Helium mass fraction is non-zero in the config file."
     write(*,*) "Please uncomment the incHe line in the Makefile or set the"
     write(*,*) "Helium mass fraction to zero in the config file."
     config_good = .false.
  end if
#endif




  if (CV%DoTestScenario) then  
     charmatch = .false.

     if ( trim(CV%TestScenario) == "iliev_test1" ) then
        charmatch = .true.
        if (CV%IsoTemp /= 1.0d4) then
           config_good = .false.
           write(*,*) "iliev test 1 must have IsoTemp = 1.0e4"
        end if
     end if

     if ( trim(CV%TestScenario) == "iliev_test1He" ) then
        charmatch = .true.
        if (CV%IsoTemp /= 1.0d4) then
           config_good = .false.
           write(*,*) "iliev test 1 He must have IsoTemp = 1.0e4"
        end if
     end if

     if ( trim(CV%TestScenario) == "iliev_test2" ) then
        charmatch = .true.
        if (CV%IsoTemp > 0.0) then
           config_good = .false.
           write(*,*) "iliev test 2 must have IsoTemp <= 0.0"
        end if
     end if

     if ( trim(CV%TestScenario) == "iliev_test3" ) then
        charmatch = .true.
        if (CV%IsoTemp > 0.0) then
           config_good = .false.
           write(*,*) "iliev test 3 must have IsoTemp <= 0.0"
        end if
     end if

     if ( trim(CV%TestScenario) == "iliev_test4" ) then
        charmatch = .true.
        if (CV%IsoTemp > 0.0) then
           config_good = .false.
           write(*,*) "iliev test 4 must have IsoTemp <= 0.0"
        end if
     end if


     if (.not. charmatch) then
        write(*,*) "TestScenario: ", trim(CV%TestScenario), " not recognized"
        write(*,*) "must be 'iliev_test1(He)', 'iliev_test2', 'iliev_test3' or 'iliev_test4' "
        config_good = .false.
     end if
  end if

  if (.not. config_good) then
     write(*,*) Cwarning
     stop
  end if


end subroutine dummy_check_config_variables




!> writes the configuration file information to the output directory
!====================================================================
subroutine config_info_to_file()

  character(200) :: config_log_file
  integer(i4b) :: lun

  config_log_file = trim(CV%OutputDir) // "/config_values_used.log"
  
  call open_formatted_file_w(config_log_file,lun)

  105 format(T2,A,I3.3)
  111 format(T2,A,I10)
  120 format(T2,A,ES10.3)

  write(lun,*)"============================================"
  write(lun,*)"The SPHRAY configuration file has been read "
  write(lun,*)"The SPHRAY configuration file variables are "
  write(lun,*)"============================================"
  write(lun,*)"Verbosity " , CV%Verbosity
  write(lun,*)
  write(lun,*)"Do a test scenario? " , CV%DoTestScenario
  write(lun,*)
  write(lun,*)"Which test scenario? " , trim(CV%TestScenario)
  write(lun,*)
  write(lun,*)"Just initialize? " , CV%JustInit
  write(lun,*)
  write(lun,*)"Input is in comoving coords? " , CV%Comoving
  write(lun,*)
  write(lun,*)"Iso temperature (if > 0.0 fixed single temperature): ", CV%IsoTemp
  write(lun,*) 
  write(lun,*)"Fix temperature at snapshot values?: ", CV%FixSnapTemp
  write(lun,*) 
  write(lun,*)"EOS temperature (negative = snapshot temperature): ", CV%EOStemp
  write(lun,*) 
  write(lun,*)"Initial xHI (negative = snapshot or collisional equil.): ", CV%InitxHI
  write(lun,*) 
  write(lun,*)"Rremove photons from ray as it travels?: ", CV%RayDepletion
  write(lun,*) 
  write(lun,*)"Integer Seed for RNG ", CV%IntSeed
  write(lun,*) 
  write(lun,*)"Static field simulation time: ", CV%StaticFieldSimTime
  write(lun,*) 
  write(lun,*)"Static field time unit: ", trim(CV%StaticSimTimeUnit)
  write(lun,*) 

  write(lun,*)

  write(lun,*)"Input Type (1=Gadget Public, 2=Gadget Cooling, 3=Gadget HDF5, 4=Gadget Bromm)", CV%InputType
  write(lun,*)
  write(lun,*)"Path to Snapshot File(s):"
  write(lun,*)trim(CV%SnapPath)
  write(lun,*)
  write(lun,*)"Path to Source File(s):"
  write(lun,*)trim(CV%SourcePath)
  write(lun,*) 
  write(lun,*)"Path to the Spectra file:"
  write(lun,*)trim(CV%SpectraFile)
  write(lun,*) 
  write(lun,*)"Path to the impact parameter -> column depth file:"
  write(lun,*)trim(CV%b2cdFile)
  write(lun,*) 
  write(lun,*)"Path to the atomic rates file:"
  write(lun,*)trim(CV%AtomicRatesFile)

  write(lun,*)
  write(lun,*)"Particle file base: ", trim(CV%ParFileBase)
  write(lun,*)"Source file base:   ", trim(CV%SourceFileBase)
  write(lun,*)
  write(lun,"(A,I3.3)") "Starting snapshot number:", CV%StartSnapNum
  write(lun,"(A,I3.3)") "Ending snapshot number:  ", CV%EndSnapNum
  write(lun,"(A,I3)") "Par Files per snapshot:    ", CV%ParFilesPerSnap
  write(lun,"(A,I10)") "Source Files per snapshot:", CV%SourceFilesPerSnap

  write(lun,*)
  write(lun,*)
  write(lun,*)  "Ray Scheme : ", trim(CV%RayScheme)
  write(lun,120) "Forced Ray Number : ", real(CV%ForcedRayNumber)

  write(lun,*) "Report ray statistics in raystats.dat?", CV%RayStats

  if(CV%BndryCond==-1) then
  write(lun,*)  "Boundry Conditions : ", "reflecting"
  else if(CV%BndryCond==0) then
  write(lun,*)  "Boundry Conditions : ", "vacuum"
  else if(CV%BndryCond==1) then  
  write(lun,*)  "Boundry Conditions : ", "periodic"
  end if

  write(lun,*)  "Ray Photon Tol     : ", CV%RayPhotonTol
  write(lun,*)  "Maximum Distance to trace a ray [physical code units], negative = default", CV%MaxRayDist

  write(lun,*)  "Use Case A recombination rates for Hydrogen? :", CV%HydrogenCaseA
  write(lun,*)  "Use Case A recombination rates for Helium?   :", CV%HeliumCaseA

  write(lun,*)  "Ionization and temperature solver :", CV%IonTempSolver

  write(lun,*)  "Temperature floor   : ", CV%Tfloor
  write(lun,*)  "Temperature ceiling : ", CV%Tceiling
  write(lun,*)  "Ionization floor    : ", CV%xfloor
  write(lun,*)  "Ionization ceiling  : ", CV%xceiling


  write(lun,*)  "ne background      : ", CV%NeBackground
  write(lun,*)  "Rays between all particle update:  ", CV%NraysUpdateNoHits

  write(lun,*) 
  write(lun,*) 
  write(lun,*)  "Hydrogen Mass Fraction: ", CV%H_mf
  write(lun,*)  "Helium Mass Fraction: ", CV%He_mf
  
  write(lun,*)
  write(lun,*)
  write(lun,*)  "Output Dir         : ", trim(CV%OutputDir)
  write(lun,*)  "Output File Base   : ", trim(CV%OutputFileBase)
  write(lun,*)  "Output Type (1=Std. Gadget, 2=HDF5 Gadget) : ", CV%OutputType

  write(lun,*)  "Output timing plan : ", trim(CV%OutputTiming)

  write(lun,*)  "Number Std Outs    : ", CV%NumStdOuts
  write(lun,*)  "Do Initial Output  : ", CV%DoInitialOutput

  write(lun,*)  "Ion Frac Out Rays  : ", CV%IonFracOutRays

  write(lun,*)  "ForcedOutFile      : ", trim(CV%ForcedOutFile)
  write(lun,*)  "ForcedUnits        : ", trim(CV%ForcedUnits)

  write(lun,*)  "Particles Per Tree Cell: ", CV%PartPerCell


  write(lun,*)"====================================="
  write(lun,*)"   End SPHRAY Configuration Output   "
  write(lun,*)"====================================="
  write(lun,*)
  write(lun,*)

  if (CV%FixSnapTemp) then
     write(lun,*) "***********************************************************"
     write(lun,*) "you are running a constant temperature simulation."
     write(lun,*) "the temperatures are fixed at the readin snapshot values"
     write(lun,*) "***********************************************************"
  else
     if (CV%IsoTemp > 0.0) then
        write(lun,*) "***********************************************************"
        write(lun,*) "you are running a constant temperature simulation."
        write(lun,*) "the temperature is fixed at T (K) = ", CV%IsoTemp
        write(lun,*) "***********************************************************"
     end if
  end if

 
  close(lun)

end subroutine config_info_to_file


end module config_mod
