!> \file ion_temperature_update.F90

!> \brief Module that is the starting point for ionization
!! and temperature updates.  It calls a runge kutta or a 
!! backwards difference formula (depending on the Makefile
!! flags)  
!<
module ion_temperature_update
use myf03_mod
use ray_mod
use raylist_mod
use ionpar_mod
use physical_constants_mod
use particle_system_mod, only: particle_system_type
use particle_system_mod, only: particle_type
use particle_system_mod, only: box_type
use oct_tree_mod, only: oct_tree_type
use spectra_mod, only: rn2freq
use euler_mod, only: eulerint, recombeulerint
use bdf_mod, only: bdfint
use atomic_rates_mod, only: get_atomic_rates
use global_mod, only: GV, saved_gheads, rtable
use global_mod, only: active_rays

use config_mod, only: CV
implicit none

private
public :: update_raylist

 
contains


!> sets the logicals that can be determined from global variables
!!-----------------------------------------------------------------
subroutine set_bools( He, caseA, isoT, fixT )
  logical, intent(out) :: He
  logical, intent(out) :: caseA(2)
  logical, intent(out) :: isoT
  logical, intent(out) :: fixT

#ifdef incHe
  He = .true.
#else
  He = .false.
#endif

  caseA = .false.
  if (CV%HydrogenCaseA) caseA(1) = .true.
  if (CV%HeliumCaseA)   caseA(2) = .true.

  if (CV%IsoTemp > 0.0) then
     isoT = .true.
  else
     isoT = .false.
  end if

  if (CV%FixSnapTemp) then
     fixT = .true.
  else
     fixT = .false.
  end if

end subroutine set_bools



!> updates the particles intersected by a ray 
!!-----------------------------------------------------------------
subroutine update_raylist(raylist, pars, box)

  type(raylist_type), intent(inout) :: raylist !< ray/particle intersections
  type(particle_type), intent(inout) :: pars(:)  !< particle system
  type(box_type), intent(in) :: box  !< particle system

  type(particle_type) :: par
  type(ionpart_type) :: ipar
  type(intersection_type) :: intersection
  type(src_ray_type) :: ray
  integer(i8b) :: impact  
  integer(i8b) :: scalls  ! number of calls to solver
  logical :: photo
  logical :: He
  logical :: caseA(2)
  logical :: isoT
  logical :: fixT


  ! set booleans
  !-------------------------------------------------------
  call set_bools( He, caseA, isoT, fixT )
  photo = .true.

  
  ! loop through the ray particle intersections
  !-------------------------------------------------------
  ray = active_rays(raylist%rayn)
  impact_loop: do impact = 1,raylist%nnb

     GV%ParticleCrossings = GV%ParticleCrossings + 1     
     intersection = raylist%intersections(impact)
     par = pars(intersection%pindx)

     

     ! check we dont have double intersections when we shouldn't
     !-------------------------------------------------------------
     if (ray%srcray) then
        if (box%tbound(1)==1 .and. ray%emit_time == par%lasthit) then
           ! here we have periodic BCs and a particle has been hit
           ! twice by the same ray so we stop tracing 
           GV%PhotonsLeavingBox = GV%PhotonsLeavingBox + ray%pcnt
           raylist%lastnnb = impact-1
           exit
        else if (box%tbound(1)==0 .and. ray%emit_time == par%lasthit) then
           ! here we have transmissive BCs and a particle has been hit
           ! twice by the same ray so something is wrong
           write(*,*) "transmissive BCs and particle hit twice in one ray!"
           write(*,*) "impact  = ", impact
           ipar%strtag = "check boundary conditions"
           call ionpar2screen(ipar)
           stop
        end if
     end if
     if (ray%emit_time < par%lasthit) then
         stop "ray emit_time < par%lasthit, causaulity broken, stopping"
     endif
     call initialize_ionpar(ipar, intersection, He)


!     write(*,*) "d,dl:", raylist%intersection(impact)%d, ipar%dl
!     write(*,"(A,4F12.6)") "pos,nH: ", ipar%pos, ipar%nH
!     write(*,*) "inside: ", ipar%inside
!     write(*,*) 

     if (ray%srcray) then
        if (CV%IonTempSolver==1) then
           call eulerint(ipar,scalls,photo,caseA,He,isoT,fixT)
           ipar%strtag = "on_eulerint_output"
        else if (CV%IonTempSolver==2) then
           call bdfint(ipar,scalls,photo,caseA,He,isoT,fixT)
           ipar%strtag = "on_bdfint_output"
        end if
     else
        call recombeulerint(ipar,scalls)
        ipar%strtag = "on_recombeulerint_output"
     end if
     call check_x(ipar)

     raylist%lastnnb = impact

     GV%TotalDerivativeCalls = GV%TotalDerivativeCalls + scalls
     if (scalls .GT. GV%PeakUpdates) GV%PeakUpdates = scalls


     !  put the updated particle data into the particle system
     !===========================================================
     call ionpar2par(ipar,par)
     if (par%T < CV%Tfloor) par%T = CV%Tfloor

#ifdef outGammaHI
     par%gammaHI = pars(ipar%index)%gammaHI + ipar%gammaHI * ipar%dt_s
     par%time    = pars(ipar%index)%time + ipar%dt_s
#endif

     if (ray%srcray) then
        par%lasthit = ray%emit_time 
     end if

     pars(intersection%pindx) = par


     !  use the solution to set some global variables
     !=================================================
     GV%TotalPhotonsAbsorbed = GV%TotalPhotonsAbsorbed + ipar%pdeps
     GV%TotalIonizations     = GV%TotalIonizations + &
                               (ipar%xHII - ipar%xHII_in) * ipar%Hcnt
#ifdef incHe
     GV%TotalIonizations     = GV%TotalIonizations + &
                               (ipar%xHeII - ipar%xHeII_in) * ipar%Hecnt     
     GV%TotalIonizations     = GV%TotalIonizations + &
                               (ipar%xHeIII - ipar%xHeIII_in) * ipar%Hecnt     
#endif


     ! if the particle satisfies the rec ray tol put it on the recomb list
     !=====================================================================
#ifdef incHrec
     ! this code is never going to work as pars(interseection%pindx) has already
     ! been tainted by par, several lines above!
     stop 'read the code here'
     if (ray%srcray) then
        if (.not. pars(intersection%pindx)%OnRecList) then
           if (pars(intersection%pindx)%xHIIrc > CV%RecRayTol) then
              pars(intersection%pindx)%OnRecList = .true.
              GV%recpt = GV%recpt + 1
              reclist(GV%recpt) = intersection%pindx
           end if
        end if
     end if
#endif
     
     
     
     !  determine if we move to the next particle and track some ray stats
     !=====================================================================

     if (CV%RayDepletion) then
        ray%pcnt = ray%pcnt - ipar%pdeps
     endif
     
     
     ! if photons are exhausted
     !-------------------------------
     if (ray%pini > 0.0) then
        if (ray%pcnt / ray%pini < CV%RayPhotonTol) then
           exit
        end if
     end if

     active_rays(raylist%rayn) = ray
     
     ! if vacuum BCs and exiting box
     !-------------------------------
     if(box%tbound(1)==0) then
        if(impact==raylist%nnb) then
           GV%PhotonsLeavingBox = GV%PhotonsLeavingBox + ray%pcnt
        end if
     end if

     
  end do impact_loop

!  stop "finished impact loop"

end subroutine update_raylist










end module ion_temperature_update
