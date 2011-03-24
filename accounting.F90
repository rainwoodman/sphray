!>
!>
!<


module accounting_mod
use myf03_mod 
implicit none
type timer_type
   real(r8b):: clock
   real(r8b):: timer
   logical(i4b):: running
   character(len=20) :: tag
endtype timer_type

type accounting_variables_type
   integer(i8b) :: TotalSourceRaysCast    !< total rays traced from user def. sources
   integer(i8b) :: TotalDiffuseRaysCast   !< total recombination rays traced
   real(r8b) :: IonizingPhotonsPerSec  !< ionizing photons emitted per second
   real(r8b) :: TotalPhotonsCast       !< total number of photons emitted
   real(r8b) :: TotalPhotonsAbsorbed   !< total number of photons absorbed
   real(r8b) :: PhotonsLeavingBox      !< total photons leaving the box
   real(r8b) :: TotalIonizations       !< total number of photoionizations
   real(r8b) :: TotalRecombinations    !< total number of recombinations
   integer(i8b) :: PeakUpdates            !< max updates for a single particle
   integer(i8b) :: ParticlesCrossed        !< total number of updated particles(with dup)
   integer(i8b) :: ParticleCrossings      !< number of ray / particle intersections, updated
   integer(i8b) :: TotalDerivativeCalls   !< times the solver being used has run
   type(timer_type) :: timers(10)
end type accounting_variables_type

type(accounting_variables_type) :: AV       !< global accounting variables

contains
subroutine reduce_accounting_variables(ja)
  type(accounting_variables_type), intent(in)  :: ja
  AV%TotalSourceRaysCast = AV%TotalSourceRaysCast + ja%TotalSourceRaysCast
  AV%TotalDiffuseRaysCast = AV%TotalDiffuseRaysCast + ja%TotalDiffuseRaysCast
  AV%IonizingPhotonsPerSec = AV%IonizingPhotonsPerSec + ja%IonizingPhotonsPerSec
  AV%TotalPhotonsCast = AV%TotalPhotonsCast + ja%TotalPhotonsCast
  AV%TotalPhotonsAbsorbed = AV%TotalPhotonsAbsorbed + ja%TotalPhotonsAbsorbed
  AV%PhotonsLeavingBox = AV%PhotonsLeavingBox + ja%PhotonsLeavingBox
  AV%TotalIonizations = AV%TotalIonizations + ja%TotalIonizations
  AV%TotalRecombinations = AV%TotalRecombinations + ja%TotalRecombinations
  AV%PeakUpdates = max(AV%PeakUpdates , ja%PeakUpdates)
  AV%ParticleCrossings = AV%ParticleCrossings + ja%ParticleCrossings
  AV%TotalDerivativeCalls = AV%TotalDerivativeCalls + ja%TotalDerivativeCalls
end subroutine reduce_accounting_variables
subroutine clear_accounting_variables(ja)
  type(accounting_variables_type), intent(inout)  :: ja
  ja%TotalSourceRaysCast = 0
  ja%TotalDiffuseRaysCast = 0
  ja%IonizingPhotonsPerSec = 0
  ja%TotalPhotonsCast = 0
  ja%TotalPhotonsAbsorbed = 0
  ja%PhotonsLeavingBox = 0
  ja%TotalIonizations = 0
  ja%TotalRecombinations = 0
  ja%PeakUpdates = 0
  ja%ParticleCrossings = 0
  ja%TotalDerivativeCalls = 0
end subroutine clear_accounting_variables

subroutine timer_init(av, timerid, n, running)
  integer(i4b),intent(in) :: timerid
  type(accounting_variables_type), intent(inout) :: av
  character(*), intent(in) :: n
  logical(i4b), intent(in) :: running
  real(r8b) OMP_GET_WTIME
  av%timers(timerid)%clock = OMP_GET_WTIME()
  av%timers(timerid)%timer = 0
  av%timers(timerid)%tag = n
  av%timers(timerid)%running = running
endsubroutine timer_init
subroutine timer_pause(av, timerid)
  integer(i4b),intent(in) :: timerid
  type(accounting_variables_type), intent(inout) :: av
  real(r8b) OMP_GET_WTIME
  av%timers(timerid)%timer = av%timers(timerid)%timer + OMP_GET_WTIME() - av%timers(timerid)%clock
  av%timers(timerid)%running = .False.
endsubroutine timer_pause
subroutine timer_resume(av, timerid)
  integer(i4b),intent(in) :: timerid
  type(accounting_variables_type), intent(inout) :: av
  real(r8b) OMP_GET_WTIME
  av%timers(timerid)%clock = OMP_GET_WTIME()
  av%timers(timerid)%running = .True.
endsubroutine timer_resume
subroutine timer_print(av, timerid)
  integer(i4b),intent(in) :: timerid
  type(accounting_variables_type), intent(inout) :: av
  real(r8b) OMP_GET_WTIME
  logical(i4b) :: was_running
  was_running = av%timers(timerid)%running
  if (was_running) then
    call timer_pause(av, timerid)
  endif
    print *, av%timers(timerid)%tag, av%timers(timerid)%timer
  if (was_running) then
    call timer_resume(av, timerid)
  endif
endsubroutine timer_print
end module accounting_mod
