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
    integer(i8b),allocatable :: pool_cur(:)
    integer(i8b),allocatable :: pool_tail(:)
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
     allocate(resolution%pool_head(size(raylists, 1)))
     allocate(resolution%pool_cur(size(raylists, 1)))
     allocate(resolution%pool_tail(size(raylists, 1)))
     total_nnb = 0
     do rayln = 1, size(raylists, 1)
       total_nnb = total_nnb + raylists(rayln)%nnb
     enddo
     
     resolution%pool_head(1) = 1
     resolution%pool_tail(1) = raylists(1)%nnb
     do rayln = 2, size(raylists, 1)
       resolution%pool_head(rayln) = resolution%pool_head(rayln - 1) + raylists(rayln - 1)%nnb
       resolution%pool_tail(rayln) = resolution%pool_tail(rayln - 1) + raylists(rayln)%nnb
     enddo
     resolution%pool_cur = resolution%pool_head

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
     resolution%good_nnb = 0
     resolution%remaining_nnb = total_nnb
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
     integer(i4b) :: rayln, raylm
     integer(i8b) :: j, k, jmpact, this, first
     good_tail = 0
     !do rayln = 1, size(raylists, 1)
     !  print *, resolution%pool_head(rayln), resolution%pool_cur(rayln), resolution%pool_tail(rayln)
     !enddo
     do rayln = 1, size(resolution%pool_cur), 1
        if (resolution%pool_cur(rayln) > resolution%pool_tail(rayln)) then 
            cycle
        endif
        c = raylists(rayln)%intersections(resolution%pool_cur(rayln) - resolution%pool_head(rayln) + 1)
        good_candidate = .True.
        do j = 1, good_tail, 1
           call resolution_get_resolved_intersection(resolution, raylists, j, a)
           if (a%pindx  == c%pindx ) then 
              good_candidate = .False.
              exit
           endif
        enddo
        if (good_candidate) then
           !print *, 'candidate', c%pindx, c%rayn, c%t
           first = resolution%pool_cur(rayln)
           this = resolution%plink(first)
           do while (this /= first)
              if(this > first) then
                this = resolution%plink(this)
                cycle
              endif
              call decode(resolution, resolution%encoded(this), raylm, jmpact)
              a = raylists(raylm)%intersections(jmpact)
              !print *, '  compare', a%pindx, a%rayn, a%t
              if(raylm >= rayln) then
                 stop "can't happen. an earlier intersection must be in an earlier ray"
              endif
              !print *, 'this', this, 'cur', resolution%pool_cur(raylm)
              if(this >= resolution%pool_cur(raylm)) then
              !  print *, 'rejected'
                 good_candidate = .False.
                 exit
              endif
              this = resolution%plink(this)
           enddo
        endif
        if (good_candidate) then
          good_tail = good_tail + 1
          resolution%good(good_tail) = encode(resolution, rayln, resolution%pool_cur(rayln) - resolution%pool_head(rayln) + 1)
          resolution%pool_cur(rayln) = resolution%pool_cur(rayln) + 1
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
    deallocate(resolution%plink)
    deallocate(resolution%encoded)
    deallocate(resolution%pool_head)
    deallocate(resolution%pool_cur)
    deallocate(resolution%pool_tail)
    resolution%secret = 0
  endsubroutine kill_resolution
end module
