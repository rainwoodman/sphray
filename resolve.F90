! the word resolution was inspired by the United Nations. - Yu
! this module contains the subroutines to plan for the parallized
! executation of temperature and ionization updating.


module resolve_mod
  
  use myf03_mod
  use raylist_mod, only: intersection_type, raylist_type
  type resolution_type
    integer(i8b),allocatable :: encoded_impacts(:)  !< tid * max_nnb + impact
    integer(i8b):: max_nnb
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
  end subroutine prepare_resolution 

  subroutine resolve(resolution, raylists, chunk)
     type(raylist_type), intent(in), allocatable :: raylists(:)
     type(resolution_type), intent(inout) :: resolution
     integer(i4b):: chunk
     integer(i4b):: NTRD, tid
     integer(i8b):: chunk_start
     type(intersection_type) :: a, b, c

  end subroutine resolve

  subroutine resolution_get_intersection(resolution, raylists, j, intersection)
     type(raylist_type), intent(in),allocatable :: raylists(:)
     type(resolution_type), intent(in) :: resolution
     type (intersection_type), intent(out) :: intersection
     integer(i8b), intent(in):: j
     tid = resolution%encoded_impacts(j) / (resolution%max_nnb + 1)
     impact = mod(resolution%encoded_impacts(j),  resolution%max_nnb + 1)
     intersection = raylists(tid)%intersections(impact)
  endsubroutine resolution_get_intersection
  subroutine kill_resolution(resolution)
    type(resolution_type), intent(inout):: resolution
    deallocate(resolution%encoded_impacts)
    resolution%max_nnb = 0
  endsubroutine kill_resolution
end module
