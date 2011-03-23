! the word resolution was inspired by the United Nations. - Yu
! this module contains the subroutines to plan for the parallized
! executation of temperature and ionization updating.


module resolve_mod
  
  use myf03_mod
  use m_mrgrnk, only: mrgrnk
  use raylist_mod, only: intersection_type, raylist_type
  implicit none
  type resolution_type
    integer(i8b):: remaining_nnb  !< nnb that remains in encoded_impacts
    integer(i8b):: secret !< used to encode the raylist id and impact id
    integer(i8b), allocatable:: encoded(:)
    integer(i8b), allocatable:: plink(:)
    integer(i8b),allocatable :: good(:)
    integer(i8b):: good_nnb
    integer(i8b),allocatable :: pool_head(:)
  end type resolution_type
contains

  function encode(resolution, rayln, impact) result(r)
    type(resolution_type), intent(in):: resolution
    integer(i4b), intent(in) :: rayln
    integer(i8b), intent(in) :: impact
    integer(i8b) :: r
    r = impact * resolution%secret + rayln
  end function encode
  subroutine decode(resolution, j, rayln, impact)
    type(resolution_type), intent(in):: resolution
    integer(i4b), intent(out) :: rayln
    integer(i8b), intent(out) :: impact
    integer(i8b), intent(in) :: j
    rayln = mod (j, resolution%secret)
    impact = j / resolution%secret
  end subroutine decode

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
     integer(i8b),allocatable :: pindx(:), indexx(:)
     integer(i8b) :: total_nnb,  j, this, last, first, impact
     integer(i4b) :: rayln
     total_nnb = 0
     do rayln = 1, size(raylists, 1)
       total_nnb = total_nnb + raylists(rayln)%nnb
     enddo
     resolution%secret = size(raylists, 1) + 1
     allocate(resolution%good(total_nnb))
     allocate(resolution%encoded(total_nnb))
     allocate(resolution%plink(total_nnb))
     ! create the link list of intersections with same particles
     allocate(pindx(total_nnb))
     allocate(indexx(total_nnb))
     j = 0
     do rayln = 1, size(raylists, 1)
       do impact = 1, raylists(rayln)%nnb
         j = j + 1
         resolution%encoded(j) = encode(resolution, rayln, impact)
         pindx(j) = raylists(rayln)%intersections(impact)%pindx
       enddo
     enddo

     call mrgrnk(pindx, indexx)

     print *, 'sorted'
     first = indexx(1)
     last =indexx(1)
     do j = 2, total_nnb, 1
       this = indexx(j)
       if (pindx(this) /= pindx(last)) then
         resolution%plink(last) = first
         first = this
         last = this
       else
         resolution%plink(last) = this
         last = this
       endif
     enddo
     resolution%plink(last) = first
     print *, 'linked'
     if (.False.) then
     !$OMP PARALLEL DO PRIVATE(first, this)
     do j = 1, total_nnb, 1
       first = resolution%plink(j)
       this = resolution%plink(j)
       last = 0
       do while (.True.)
         if(pindx(j) /= pindx(this)) then
             stop 'failed pindx list check'
         endif
         this = resolution%plink(this)
         last = last + 1
         if(this == first) exit
       enddo
     enddo
     !$OMP END PARALLEL DO
     endif
     deallocate(pindx)
     deallocate(indexx)

     print *, 'checked'
     allocate(resolution%pool_head(size(raylists, 1)))
     resolution%good_nnb = 0
     resolution%remaining_nnb = total_nnb
     resolution%pool_head(:) = 1
  end subroutine prepare_resolution 

  subroutine resolve_all(resolution, raylists)
     type(raylist_type), intent(in), allocatable :: raylists(:)
     type(resolution_type), intent(inout) :: resolution
     integer(i8b) ::  j, k
     integer(i4b) :: i
     j = 0
     do i = 1, size(raylists, 1)
       do k = 1, raylists(i)%nnb
         j = j + 1
         resolution%good(j) = encode(resolution, i, k)
       enddo
     enddo
     resolution%good_nnb = j
     resolution%remaining_nnb = 0
  end subroutine resolve_all

  subroutine resolve_more(resolution, raylists)
     type(raylist_type), intent(in), allocatable :: raylists(:)
     type(resolution_type), intent(inout) :: resolution
     integer(i8b) :: good_tail
     type (intersection_type) :: a, c
     logical(i4b) :: good_candidate
     integer(i4b) :: i
     integer(i8b) :: j, k
     good_tail = 0
     do i = 1, size(resolution%pool_head), 1
        if (resolution%pool_head(i) > raylists(i)%nnb) then 
            cycle
        endif
        c = raylists(i)%intersections(resolution%pool_head(i))
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
             do k = resolution%pool_head(j), raylists(j)%nnb, 1
                if (raylists(j)%intersections(k)%pindx == c%pindx) then
                   good_candidate = .False.
                   exit outer
                endif
             enddo
           enddo outer
        endif
        if (good_candidate) then
          good_tail = good_tail + 1
          resolution%good(good_tail) = encode(resolution, i, resolution%pool_head(i))
          resolution%pool_head(i) = resolution%pool_head(i) + 1
        else
          exit
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
     integer(i8b):: impact
     integer(i4b):: rayln
     call decode(resolution, resolution%good(j), rayln, impact)
     intersection = raylists(rayln)%intersections(impact)
  endsubroutine resolution_get_resolved_intersection
  subroutine kill_resolution(resolution)
    type(resolution_type), intent(inout):: resolution
    deallocate(resolution%good)
    deallocate(resolution%pool_head)
    resolution%secret = 0
  endsubroutine kill_resolution
end module
