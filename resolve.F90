! the word resolution was inspired by the United Nations. - Yu
! this module contains the subroutines to plan for the parallized
! executation of temperature and ionization updating.


module resolve_mod
  
  use myf03_mod
  use raylist_mod, only: intersection_type, raylist_type
  type resolution_type
    integer(i8b),allocatable :: encoded_impacts(:)  !< tid * max_nnb + impact
    integer(i8b):: remaining_nnb  !< nnb that remains in encoded_impacts
    integer(i8b):: secret !< used to encode the raylist id and impact id
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
     integer(i8b) :: total_nnb,  j
     integer(i4b) :: rlid
     total_nnb = 0
     do rlid = 1, size(raylists, 1)
       total_nnb = total_nnb + raylists(rlid)%nnb
     enddo
     resolution%secret = size(raylists, 1) + 1
     allocate(resolution%good(total_nnb))
     resolution%good_nnb = 0
     resolution%remaining_nnb = total_nnb
  end subroutine prepare_resolution 

  subroutine resolve_more(resolution, raylists)
     type(raylist_type), intent(in), allocatable :: raylists(:)
     type(resolution_type), intent(inout) :: resolution
     integer(i8b) :: good_tail
     integer(i8b) :: pool_head(size(raylists, 1))
     type (intersection_type) :: a, c
     logical(i4b) :: good_candidate
     integer(i8b) :: i, j, k
     good_tail = 0
     pool_head = 1

     j = 0
     do i = 1, size(raylists, 1)
       do k = 1, raylists(i)%nnb
         j = j + 1
         resolution%good(j) = k* resolution%secret + i
       enddo
     enddo
     resolution%good_nnb = j
     resolution%remaining_nnb = 0
     return 
     do i = 1, size(pool_head), 1
        if (pool_head(i) > raylists(i)%nnb) cycle
        c = raylists(i)%intersections(pool_head(i))
        good_candidate = .True.
        do j = 1, good_tail, 1
           call resolution_get_resolved_intersection(resolution, raylists, j, a)
           if (a%pindx  == c%pindx ) then 
              good_candidate = .False.
              exit
           endif
        enddo
        if (good_candidate) then
           outer: do j = 1, i - 1, 1
             do k = pool_head(j), raylists(j)%nnb, 1
                if (raylists(j)%intersections(k)%pindx == c%pindx) then
                   good_candidate = .False.
                   exit outer
                endif
             enddo
           enddo outer
        endif
        if (good_candidate) then
          good_tail = good_tail + 1
          resolution%good(good_tail) = resolution%secret * pool_head(i) + i
          pool_head(i) = pool_head(i) + 1
        endif
     enddo

     resolution%good_nnb = good_tail
     resolution%remaining_nnb = resolution%remaining_nnb - good_tail
  end subroutine resolve_more

  subroutine resolution_get_resolved_intersection(resolution, raylists, j, intersection)
     type(raylist_type), intent(in),allocatable :: raylists(:)
     type(resolution_type), intent(in) :: resolution
     type (intersection_type), intent(out) :: intersection
     integer(i8b), intent(in):: j
     impact = resolution%good(j) / (resolution%secret)
     tid = mod(resolution%good(j),  resolution%secret)
     intersection = raylists(tid)%intersections(impact)
  endsubroutine resolution_get_resolved_intersection
  subroutine kill_resolution(resolution)
    type(resolution_type), intent(inout):: resolution
    deallocate(resolution%good)
    resolution%secret = 0
  endsubroutine kill_resolution
end module
