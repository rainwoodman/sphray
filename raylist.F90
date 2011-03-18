!> \file raylist.F90

!> \brief The raylist module 
!!
!<

module raylist_mod
use myf03_mod
use ray_mod
use particle_system_mod, only: particle_system_type
use particle_system_mod, only: particle_type
use particle_system_mod, only: transformation_type
use oct_tree_mod, only: oct_tree_type
use global_mod, only: active_rays, psys

implicit none

 private
 public :: intersection_type
 public :: trace_ray
 public :: raylist_type
 public :: prepare_raysearch
 public :: kill_raylist

 real, parameter :: zero = 0.0d0
 real, parameter :: half = 0.5d0

 integer,parameter :: MAX_RAYLIST_LENGTH = 1000000  !< default maximum


!> holds a particle index, an impact parameter, and a distance along a ray
!! usefull for the ionization and temperature updating
! -----------------------------------------------------------------------
   type intersection_type
      integer(i8b) :: pindx   !< particle index
      integer(i4b) :: rayn    !< ray index within activa_rays
      real :: b          !< impact parameter
      real :: d          !< distance along ray
      real :: dl         !< path length
   end type intersection_type
 
!> grouping of all things ray + impacts
!--------------------------------------- 
   type raylist_type
      integer :: nnb              !< number of intersections
      type(intersection_type), allocatable :: intersections(:) !< ray/par 
   end type raylist_type



contains


 
!> set intersection values
!-----------------------------------------
subroutine set_intersection(intersection, curay, rayn, pindx)
   type(intersection_type), intent(out) :: intersection !< the intersection
   type(src_ray_type) :: curay !< the ray
   integer(i4b) :: rayn         !< the ray
   integer(i8b)  :: pindx            !< the particle index
   real(r8b) :: d
   
     intersection%pindx = pindx
     intersection%b = src_ray_dist2pt(curay, psys%par(pindx)%pos, d)
     intersection%d = d
     intersection%rayn = rayn
 end subroutine set_intersection

 
!> initialize raylist variables and set search images
!------------------------------------------------------
 subroutine prepare_raysearch(psys, raylist)

   type(particle_system_type) psys !< particle system
   type(raylist_type) raylist      !< ray list
  
   raylist%nnb            = 0

   allocate(raylist%intersections(MAX_RAYLIST_LENGTH))

 end subroutine prepare_raysearch


!> check all of the search images for intersection
!---------------------------------------------------
 subroutine fullsearch(psys, searchtree, rayn, raylist)

   integer(i4b), intent(in) :: rayn !< active_rays(rayn) is to be traced
   type(particle_system_type) :: psys  !< particle system
   type(oct_tree_type) :: searchtree   !< oct-tree to search
   type(raylist_type) raylist          !< raylist
   type(src_ray_type) :: curay !< transformed ray
   integer(i4b) :: searchimage
   integer(i8b) :: searchcell  !< index of cell being searched
   
   do searchimage = lbound(psys%box%trafo, 1), ubound(psys%box%trafo, 1)
      searchcell = 1
      call src_ray_transform( active_rays(rayn), curay, psys%box%trafo(searchimage) )
      call raysearch(psys, searchtree, curay, rayn, raylist, searchcell)
! this is strange, if raylist%searchcell != 0 something must be wrong,
! shall we print an error instead?
      if (searchcell /= 0) return 
   enddo

 end subroutine fullsearch


!> checks a single search image for intersection
!-----------------------------------------------
 subroutine raysearch(psys, searchtree, curay, rayn, raylist, searchcell)

   type(src_ray_type), intent(in) :: curay !< transformed ray
   integer(i4b), intent(in) :: rayn    !< activerays[rayn] will be traced
   type(particle_system_type), intent(in) :: psys         !< particle system
   type(oct_tree_type), intent(in), target :: searchtree  !< oct-tree to search
   type(raylist_type) :: raylist                          !< raylist
   integer(i8b), intent(inout) :: searchcell  !< index of cell being searched


   type(oct_tree_type), pointer :: tree      !< pointer to tree
   integer(i8b) :: this, daughter, next      !< octree indices
   integer(i8b) :: par_in_cell               !< number of particles in current search cell
   logical :: par_hit                        !< ray / particle intersection test
   logical :: cell_hit                       !< ray / cell intersection test
   logical :: long                           !< past max distance? 
   integer(i8b) :: i, orderindx

   tree => searchtree

   ! curay%start = ray%start * fac + shift
   next = searchcell

   do while (next /= 0)

      this     = next
      daughter = tree%cell(this)%daughter
      next     = tree%cell(this)%next
           
      ! if we've reached a leaf
      !----------------------------
      if (daughter == 0) then

         ! return if we go over max intersections
         !--------------------------------------------
         par_in_cell = tree%cell(next)%start - tree%cell(this)%start
         if (raylist%nnb + par_in_cell > size(raylist%intersections, 1)) then             
            write(*,*) ' *** reached max intersections *** '
            searchcell = this
            return            
         endif
         
         ! add intersected particles to list
         !----------------------------------------------
         do i = tree%cell(this)%start, tree%cell(next)%start - 1
            orderindx = tree%partorder(i)
            par_hit = src_ray_part_intersection( curay, psys%par(orderindx) ) 
            if (par_hit) then
               raylist%nnb = raylist%nnb + 1
               call set_intersection(raylist%intersections(raylist%nnb), curay, rayn, orderindx)
            endif
         enddo


      ! if we need to descend further
      !--------------------------------
      else

         cell_hit = src_ray_cell_intersection( curay, tree%cell(this) )
         if ( cell_hit ) next = daughter  

      endif

   enddo
   
   searchcell = next
   
 end subroutine raysearch


!> kill a raylist
!---------------------------------
 subroutine kill_raylist(raylist)

   type(raylist_type) :: raylist !< the raylist to kill
  
     raylist%nnb=0
     deallocate(raylist%intersections)

 end subroutine kill_raylist



!> error handling
!---------------------------------
 subroutine raylistError(string,i)
   character(*) :: string  !< error message
   integer, optional :: i  !< error number

     print*,' Error detected:'
  
     if(present(i)) then
        print*,string,i
     else
        print*,string
     endif
  
     stop

 end subroutine raylistError

!> sorts a raylist according to the distance along the ray with the particles 
!! closest to the origin of the ray first in the list.   The corresponding 
!! changes are made to pindx and b
!---------------------------------------------------------------------------
 subroutine sort3_raylist(raylist)
 use m_mrgrnk, only: mrgrnk
 implicit none

   type(raylist_type) :: raylist     !< raylist to sort

   integer(i8b) :: N 
   integer(i8b) :: indexx(raylist%nnb)
   real(r8b) :: darr(raylist%nnb)

   N = raylist%nnb
   darr(1:N)=raylist%intersections(1:N)%d
   call mrgrnk(darr,indexx)

   raylist%intersections(1:N)%d=raylist%intersections(indexx(1:N))%d  
   raylist%intersections(1:N)%b=raylist%intersections(indexx(1:N))%b 
   raylist%intersections(1:N)%pindx=raylist%intersections(indexx(1:N))%pindx  
   raylist%intersections(1:N)%rayn=raylist%intersections(indexx(1:N))%rayn

 end subroutine sort3_raylist





!> given a ray creates a raylist with intersections
!------------------------------------------------------
 subroutine trace_ray(rayn, raylist, psys, searchtree) 
   integer(i4b), intent(in) :: rayn   !< the ray number to be traced (active_rays[rayn])
   type(raylist_type), intent(inout) :: raylist    !< the returned raylist
   type(particle_system_type), intent(in) :: psys  !< the particle system
   type(oct_tree_type), intent(in) :: searchtree   !< the oct-tree to search

   logical :: wantsort

   call fullsearch(psys, searchtree, rayn, raylist)

   call sort3_raylist(raylist)

 end subroutine trace_ray




end module raylist_mod
