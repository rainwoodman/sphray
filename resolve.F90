! the word resolution was inspired by the United Nations. - Yu
! this module contains the subroutines to plan for the parallized
! executation of temperature and ionization updating.


module resolve_mod
  
  use myf03_mod
  use raylist_mod, only: intersection_type, raylist_type
  type resolution_type
    integer(i8b),allocatable :: encoded_impacts(:)  !< tid * max_nnb + impact
    integer(i8b):: remaining_nnb  !< nnb that remains in encoded_impacts
    integer(i8b):: max_nnb
    integer(i8b),allocatable :: good(:)
    integer(i8b):: good_nnb
  end type resolution_type
contains

  function race(a, b) result(r)
     type (intersection_type), intent(in) :: a, b
     logical(i4b) :: r
     if(a%pindx == b%pindx .or. a%rayn == b%rayn) then 
        r = .True.
     else 
        r = .False.
     endif
  end function race
  function cconflict(a, b) result(r)
     type (intersection_type), intent(in) :: a, b
     logical(i4b) :: r
     r = a%t > b%t .and. race(a,b)
  end function cconflict

  subroutine prepare_resolution(resolution, raylists)
     type(raylist_type), intent(in), allocatable :: raylists(:)
     type(resolution_type), intent(inout) :: resolution
     integer(i8b) :: total_nnb, max_nnb, j
     integer(i4b) :: tid
     total_nnb = 0
     max_nnb = 0
     j = 1
     do tid = 0, size(raylists, 1) - 1
       total_nnb = total_nnb + raylists(tid)%nnb
       max_nnb = max(max_nnb, raylists(tid)%nnb)
       print *, size(raylists, 1), tid, total_nnb, max_nnb
     enddo
     allocate(resolution%encoded_impacts(total_nnb))
     resolution%max_nnb = max_nnb
     do tid = 0, size(raylists, 1)- 1
       do impact = 1, raylists(tid)%nnb
         resolution%encoded_impacts(j) = tid * (max_nnb + 1)+ impact
         j = j + 1
       enddo
     enddo
     allocate(resolution%good(total_nnb))
     resolution%good_nnb = 0
     resolution%remaining_nnb = total_nnb
  end subroutine prepare_resolution 

  subroutine resolve_more(resolution, raylists)
     type(raylist_type), intent(in), allocatable :: raylists(:)
     type(resolution_type), intent(inout) :: resolution
     integer(i8b) :: bad(size(resolution%encoded_impacts, 1))
     integer(i8b) :: good_tail
     integer(i8b) :: bad_tail
     integer(i8b) :: pool_head
     type (intersection_type) :: a, b, c
     logical(i4b) :: good_candidate
     integer(i8b) :: candidate_index
     integer(i8b) :: i
     bad_tail = 0
     good_tail = 0
     pool_head = 1
     do while (pool_head <= resolution%remaining_nnb)
       good_candidate = .True.
       candidate_index = pool_head
       ! pop up the head from pool, to c
       call resolution_get_intersection(resolution, raylists, pool_head, c)
       pool_head = pool_head + 1
       do i = bad_tail, 1, -1
         call resolution_get_intersection(resolution, raylists, bad(i), a)
         if (cconflict(c, a)) then
           good_candidate = .False.
           exit
         endif
       enddo
       do i = good_tail, 1, -1
         call resolution_get_intersection(resolution, raylists, resolution%good(i), a)
         if (race(a, c)) then 
           good_candidate = .False.
           exit
         endif
       enddo
       if (good_candidate) then
          good_tail = good_tail + 1
          resolution%good(good_tail) = candidate_index
       else
          bad_tail = bad_tail + 1
          bad(bad_tail) = candidate_index
       endif
       
       if(good_tail > 1000) then
!         print *, 'stop planning, good/bad = ', good_tail, bad_tail
!         do i = 1, bad_tail, 1
!            call resolution_get_intersection(resolution, raylists, bad(i), a)
!            print *, 'b', a%t, a%rayn, a%pindx
!         enddo
!         do i = 1, good_tail, 1
!            call resolution_get_intersection(resolution, raylists, resolution%good(i), a)
!            print *, 'g', a%t, a%rayn, a%pindx
!         enddo
         exit
       endif
     enddo
     ! convert the index to encoded_impacts
     do i = 1, good_tail, 1
       resolution%good(i) = resolution%encoded_impacts(resolution%good(i))
     enddo
     do i = 1, bad_tail, 1
       ! use bad as the intermediate array to avoid conflicts
       bad(i) = resolution%encoded_impacts(bad(i))
     enddo
     !  put the bad queue to the begining of the pool
     do i = 1, bad_tail, 1
       resolution%encoded_impacts(i) = bad(i)
     enddo
     !early termination, thus move the remaining of the pool forward
     do i = pool_head, resolution%remaining_nnb, 1
       ! bad_tail + 1 always <= pool_head, no worry here. compiler happy? at least icc can vectorize this
       resolution%encoded_impacts(bad_tail + 1 + i - pool_head) = resolution%encoded_impacts(i)
     enddo
     resolution%remaining_nnb = bad_tail + resolution%remaining_nnb - pool_head + 1
     resolution%good_nnb = good_tail
  end subroutine resolve_more

  subroutine resolution_get_intersection(resolution, raylists, j, intersection)
     type(raylist_type), intent(in),allocatable :: raylists(:)
     type(resolution_type), intent(in) :: resolution
     type (intersection_type), intent(out) :: intersection
     integer(i8b), intent(in):: j
     tid = resolution%encoded_impacts(j) / (resolution%max_nnb + 1)
     impact = mod(resolution%encoded_impacts(j),  resolution%max_nnb + 1)
     intersection = raylists(tid)%intersections(impact)
  endsubroutine resolution_get_intersection

  subroutine resolution_get_resolved_intersection(resolution, raylists, j, intersection)
     type(raylist_type), intent(in),allocatable :: raylists(:)
     type(resolution_type), intent(in) :: resolution
     type (intersection_type), intent(out) :: intersection
     integer(i8b), intent(in):: j
     tid = resolution%good(j) / (resolution%max_nnb + 1)
     impact = mod(resolution%good(j),  resolution%max_nnb + 1)
     intersection = raylists(tid)%intersections(impact)
  endsubroutine resolution_get_resolved_intersection
  subroutine kill_resolution(resolution)
    type(resolution_type), intent(inout):: resolution
    deallocate(resolution%encoded_impacts)
    deallocate(resolution%good)
    resolution%max_nnb = 0
  endsubroutine kill_resolution
end module
