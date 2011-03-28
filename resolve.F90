! the word resolution was inspired by the United Nations. - Yu
! this module contains the subroutines to plan for the parallized
! executation of temperature and ionization updating.


module resolve_mod
  
  use myf03_mod
!  use rbtree_mod, only: rbtree_type, rbtree_prepare, rbtree_kill, rbtree_eject, rbtree_return, rbtree_insert, rbtree_print
  use m_mrgrnk, only: mrgrnk
  use raylist_mod, only: intersection_type, raylist_type
  use particle_system_mod, only: particle_type
  use accounting_mod
  implicit none
 integer(i8b), parameter :: MAX_GOOD_LENGTH = 10000
  type resolution_type
    integer(i8b):: remaining_nnb  !< nnb that remains in encoded_impacts
    integer(i8b):: secret !< used to encode the raylist id and impact id
    integer(i8b), allocatable:: encoded(:)
    integer(i8b), allocatable:: plink_r(:)  !< circular link on the same particle, reversed time order
    integer(i8b), allocatable:: plink(:)  !< circular link on the same particle, in time order
    integer(i8b) :: good(MAX_GOOD_LENGTH)
    integer(i8b) :: good_pindx(MAX_GOOD_LENGTH)
    logical(i1b), allocatable :: par_in_good(:)
    integer(i8b):: good_nnb
    integer(i8b),allocatable :: pool_head(:)
    integer(i8b),allocatable :: pool_cur(:)
    integer(i8b),allocatable :: pool_tail(:)
    integer(i8b),allocatable :: priority(:)
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

  subroutine prepare_resolution(resolution, raylists, par)
     type(raylist_type), intent(in), allocatable :: raylists(:)
     type(particle_type), intent(in), allocatable :: par(:)
     type(resolution_type), intent(inout) :: resolution
     integer(i8b),allocatable :: pindx(:), indexx(:)
     integer(i8b) :: total_nnb,  j, this, last, first, impact
     integer(i4b) :: rayln
     allocate(resolution%par_in_good(size(par, 1)))
     do j = 1, size(par, 1), 1
        resolution%par_in_good(j) = .False.
     enddo
     allocate(resolution%pool_head(size(raylists, 1)))
     allocate(resolution%pool_cur(size(raylists, 1)))
     allocate(resolution%pool_tail(size(raylists, 1)))
     allocate(resolution%priority(size(raylists, 1)))
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

     do rayln = 1, size(raylists, 1), 1
        resolution%priority(rayln) = resolution%pool_cur(rayln) - resolution%pool_tail(rayln)
     enddo

     resolution%secret = size(raylists, 1) + 1
     allocate(resolution%encoded(total_nnb))
     allocate(resolution%plink_r(total_nnb))
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

     ! first categroy by the particle ids
     ! note that mgrrnk is stable, thus the sorted array
     ! preserves the time ordering
     call mrgrnk(pindx, indexx)

     ! create the reversed link
     first = indexx(total_nnb)
     last =indexx(total_nnb)
     do j = total_nnb - 1, 1, -1
       this = indexx(j)
       if (pindx(this) /= pindx(last)) then
         resolution%plink_r(last) = first
         first = this
         last = this
       else
         if (this >= last) then
           stop "sort subroutine unstable!"
         endif
         resolution%plink_r(last) = this
         last = this
       endif
     enddo
     resolution%plink_r(last) = first

     ! create the forward link
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
     print *, 'particle links created'

     do j = 1, total_nnb, 1
       if (resolution%plink_r(j) == 0) then
          stop "plink_r failed"
       endif
     enddo
    
     do j = 1, total_nnb, 1
       if (resolution%plink(j) == 0) then
          stop "plink failed"
       endif
     enddo
     deallocate(pindx)
     deallocate(indexx)

     resolution%good_nnb = 0
     resolution%remaining_nnb = total_nnb
     total_nnb = size(raylists, 1)

  end subroutine prepare_resolution 

  subroutine resolve_all(resolution, raylists)
     type(raylist_type), intent(in), allocatable :: raylists(:)
     type(resolution_type), intent(inout) :: resolution
     integer(i8b) ::  j, k
     integer(i4b) :: i
     j = resolution%remaining_nnb 
     j = size(resolution%encoded, 1) - resolution%remaining_nnb + 1
     do k = 1, MAX_GOOD_LENGTH
         resolution%good(k) = k + j - 1
         if(k + j -1 == size(resolution%encoded, 1) .or. k == MAX_GOOD_LENGTH) then 
            exit
         endif
     enddo
     resolution%good_nnb = k
     resolution%remaining_nnb = resolution%remaining_nnb - k
  end subroutine resolve_all

  subroutine resolve_more(resolution, raylists)
     type(raylist_type), intent(in), allocatable :: raylists(:)
     type(resolution_type), intent(inout) :: resolution
     integer(i8b) :: good_tail
     type (intersection_type) :: a, c
     logical(i4b) :: good_candidate
     integer(i4b) :: rayln, raylm
     integer(i8b) :: j, k, jmpact, this, first, next, prev
     integer(i8b) :: key, value
     integer(i8b) :: counter
     integer(i8b) :: c_pindx
     integer(i8b) :: handled
     integer(i8b) :: rayindex(size(raylists, 1))
     good_tail = 0
     !do rayln = 1, size(raylists, 1)
     !  print *, resolution%pool_head(rayln), resolution%pool_cur(rayln), resolution%pool_tail(rayln)
     !enddo
     !call mrgrnk(resolution%priority, rayindex)
     
     !print *, rayindex(1), resolution%priority(rayindex(1))
     handled = 1
     do rayln = 1, size(raylists, 1), 1
        if (resolution%pool_cur(rayln) > resolution%pool_tail(rayln)) then 
            cycle
        endif
        c_pindx = raylists(rayln)%intersections(resolution%pool_cur(rayln) - resolution%pool_head(rayln) + 1)%pindx
        good_candidate = .True.
        call timer_resume(AV, 3)
        if (resolution%par_in_good(c_pindx)) then
              good_candidate = .False.
        endif
        call timer_pause(AV, 3)
        call timer_resume(AV, 4)
        if (good_candidate) then
           counter = 0
           !print *, 'candidate', c%pindx, c%rayn, c%t
           first = resolution%pool_cur(rayln)
           this = resolution%plink_r(first)
           do while (this /= first)
              if(this > first) then
                 exit ! because the link list is reverse sorted
                      ! if 'this' > 'first', first is already the
                      ! earliest impact on the particle
                      ! we can move it wherever we want
              endif
              call decode(resolution, resolution%encoded(this), raylm, jmpact)
              !a = raylists(raylm)%intersections(jmpact)
              !print *, '  compare', a%pindx, a%rayn, a%t
              if(raylm >= rayln) then
                 stop "can't happen. an earlier intersection must be in an earlier ray"
              endif
              !print *, 'this', this, 'cur', resolution%pool_cur(raylm)
              if(this >= resolution%pool_cur(raylm)) then
              !  print *, 'rejected'
!                 resolution%priority(raylm) = resolution%priority(raylm) + rayindex(rayln)
                 good_candidate = .False.
                 exit
              endif
              this = resolution%plink_r(this)
              counter = counter + 1
           enddo
        endif
        call timer_pause(AV, 4)
        if (good_candidate) then
          good_tail = good_tail + 1
          this = resolution%pool_cur(rayln)
          next = resolution%plink(this)
          prev = resolution%plink_r(this)
          resolution%plink(prev) = next
          resolution%plink_r(next) = prev
          resolution%plink(this) = this
          resolution%plink_r(this) = this
          resolution%good(good_tail) = this
          resolution%par_in_good(c_pindx) = .True.
          resolution%good_pindx(good_tail) = c_pindx
          resolution%pool_cur(rayln) = resolution%pool_cur(rayln) + 1
          resolution%priority(rayln) = resolution%pool_cur(rayln) - resolution%pool_tail(rayln)
!          print *, 'taken', rayln, '1', resolution%pool_cur(1) - resolution%pool_tail(1)
        endif
        if (good_tail == MAX_GOOD_LENGTH) then 
          exit
        endif
     enddo

     do j = 1, good_tail, 1
       resolution%par_in_good(resolution%good_pindx(j)) = .False.
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
     call decode(resolution, resolution%encoded(resolution%good(j)), rayln, impact)
     intersection = raylists(rayln)%intersections(impact)
  endsubroutine resolution_get_resolved_intersection
  subroutine kill_resolution(resolution)
    type(resolution_type), intent(inout):: resolution
    deallocate(resolution%plink_r)
    deallocate(resolution%priority)
    deallocate(resolution%plink)
    deallocate(resolution%encoded)
    deallocate(resolution%pool_head)
    deallocate(resolution%pool_cur)
    deallocate(resolution%pool_tail)
    deallocate(resolution%par_in_good)
    resolution%secret = 0
  endsubroutine kill_resolution
end module
