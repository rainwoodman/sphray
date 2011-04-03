! the word resolution was inspired by the United Nations. - Yu
! this module contains the subroutines to plan for the parallized
! executation of temperature and ionization updating.


module resolve_mod
  
  use myf03_mod
!  use rbtree_mod, only: rbtree_type, rbtree_prepare, rbtree_kill, rbtree_eject, rbtree_return, rbtree_insert, rbtree_print
  use m_mrgrnk, only: mrgrnk
  use ray_mod, only: src_ray_type
  use raylist_mod, only: intersection_type, raylist_type
  use particle_system_mod, only: particle_type
  use accounting_mod
  implicit none
 integer(i8b), parameter :: MAX_GOOD_LENGTH = 10000
  type resolution_type
    integer(i8b):: remaining_nnb  !< nnb that remains in encoded_impacts
    integer(i8b):: secret!< encoding rayln and impact, fortran has no pointer
    integer(i8b), allocatable:: encoded(:)
    integer(i8b), allocatable:: plink_r(:)  !< circular link on the same particle, reversed time order
    integer(i8b), allocatable:: plink(:)  !< circular link on the same particle, in time order
    integer(i8b),allocatable :: rhead(:)
    integer(i8b), allocatable:: rlink(:)  !< link on the same ray, in time order
    integer(i8b) :: good(MAX_GOOD_LENGTH)
    integer(i8b) :: good_pindx(MAX_GOOD_LENGTH)
    logical(i1b), allocatable :: par_in_good(:)
    integer(i8b):: good_nnb
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
  subroutine build_list(key, indexx, link, head)
     integer(i8b),allocatable, intent(inout) :: key(:), indexx(:), link(:)
     integer(i8b),allocatable, optional, intent(inout):: head(:)
     integer(i8b) :: first, last, total_nnb, j, this
     total_nnb = size(key, 1)
     first = indexx(1)
     last =indexx(1)
     if (present(head)) then
       head(:) = 0
       head(key(first)) = first
     endif
     do j = 2, total_nnb, 1
       this = indexx(j)
       if (key(this) /= key(last)) then
         if (present(head)) then
         head(key(this)) = this
         endif
         link(last) = 0
         first = this
         last = this
       else
         link(last) = this
         last = this
       endif
     enddo
     link(last) = 0
     print *, 'forward link created'
    
  endsubroutine build_list

  subroutine build_circular_list(key, indexx, link, link_r, head, tail)
     integer(i8b),allocatable, intent(inout) :: key(:), indexx(:), link(:), link_r(:)
     integer(i8b),allocatable, optional, intent(inout):: head(:), tail(:)
     integer(i8b) :: first, last, total_nnb, j, this
     total_nnb = size(key, 1)
     ! create the reversed link
     first = indexx(total_nnb)
     last =indexx(total_nnb)
     if (present(tail)) then
       tail(key(last)) = last
     endif
     do j = total_nnb - 1, 1, -1
       this = indexx(j)
       if (key(this) /= key(last)) then
         if (present(tail)) then
           tail(key(this)) = this
         endif
         link_r(last) = first
         first = this
         last = this
       else
         if (this >= last) then
           stop "sort subroutine unstable!"
         endif
         link_r(last) = this
         last = this
       endif
     enddo
     link_r(last) = first

     ! create the forward link
     first = indexx(1)
     last =indexx(1)
     if (present(head)) then
       head(key(first)) = first
     endif
     do j = 2, total_nnb, 1
       this = indexx(j)
       if (key(this) /= key(last)) then
         if (present(head)) then
         head(key(this)) = this
         endif
         link(last) = first
         first = this
         last = this
       else
         link(last) = this
         last = this
       endif
     enddo
     link(last) = first
     print *, 'particle links created'

     do j = 1, total_nnb, 1
       if (link_r(j) == 0) then
          stop "plink_r failed"
       endif
     enddo
    
     do j = 1, total_nnb, 1
       if (link(j) == 0) then
          stop "plink failed"
       endif
     enddo
  endsubroutine build_circular_list
  subroutine prepare_resolution(resolution, raylists, ray, par)
     type(raylist_type), intent(in), allocatable :: raylists(:)
     type(particle_type), intent(in), allocatable :: par(:)
     type(src_ray_type), intent(in), allocatable :: ray(:)
     type(resolution_type), intent(inout) :: resolution
     real(r8b), allocatable :: t(:)
     integer(i8b),allocatable :: pindx(:), indexx(:), rayn(:)
     integer(i8b) :: total_nnb,  j, this, last, first, impact
     integer(i4b) :: rayln
     allocate(resolution%par_in_good(size(par, 1)))
     do j = 1, size(par, 1), 1
        resolution%par_in_good(j) = .False.
     enddo

     total_nnb = 0
     do rayln = 1, size(raylists, 1)
       total_nnb = total_nnb + raylists(rayln)%nnb
     enddo
     
     resolution%secret = size(raylists, 1) + 1
     allocate(resolution%encoded(total_nnb))
     allocate(resolution%plink_r(total_nnb))
     allocate(resolution%plink(total_nnb))
     allocate(resolution%rlink(total_nnb))
     allocate(resolution%rhead(size(ray, 1)))
     ! create the link list of intersections with same particles
     allocate(pindx(total_nnb))
     allocate(rayn(total_nnb))
     allocate(indexx(total_nnb))
     allocate(t(total_nnb))
     j = 0
     do rayln = 1, size(raylists, 1)
       do impact = 1, raylists(rayln)%nnb
         j = j + 1
         resolution%encoded(j) = encode(resolution, rayln, impact)
         t(j) = raylists(rayln)%intersections(impact)%t
       enddo
     enddo

     call mrgrnk(t, indexx)
     resolution%encoded = resolution%encoded(indexx)
     t = t(indexx)

     do j= 1, total_nnb, 1
       call decode(resolution, resolution%encoded(j), rayln, impact)
       pindx(j) = raylists(rayln)%intersections(impact)%pindx
       rayn(j) = raylists(rayln)%intersections(impact)%rayn
     enddo

     ! first categroy by the particle ids
     ! note that mgrrnk is stable, thus the sorted array
     ! preserves the time ordering
     call mrgrnk(pindx, indexx)
     call build_circular_list(pindx, indexx, resolution%plink, resolution%plink_r)

     call mrgrnk(rayn, indexx)
     call build_list(rayn, indexx, resolution%rlink, resolution%rhead)
     deallocate(indexx)
     deallocate(pindx)
     deallocate(rayn)
     deallocate(t)

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
     integer(i4b) :: rayln, raylm, rayn
     integer(i8b) :: j, k, jmpact, this, first, next, prev, impact
     integer(i8b) :: key, value
     integer(i8b) :: counter
     integer(i8b) :: c_pindx
     integer(i8b) :: handled
     integer(i8b) :: rayindex(size(raylists, 1))
     good_tail = 0

     do rayn = 1, size(resolution%rhead, 1), 1
        if (resolution%rhead(rayn) == 0 ) then 
            cycle
        endif
        this = resolution%rhead(rayn)
        call decode(resolution, resolution%encoded(this), rayln, impact)
        c_pindx = raylists(rayln)%intersections(impact)%pindx
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
           prev = resolution%plink_r(this)
           if (prev < this) then
             ! if this is not the first update in the particle link list
             good_candidate = .False.
           endif
        endif
        call timer_pause(AV, 4)
        if (good_candidate) then
          good_tail = good_tail + 1
          ! remove this from the particle link list
          next = resolution%plink(this)
          prev = resolution%plink_r(this)
          resolution%plink(prev) = next
          resolution%plink_r(next) = prev
          resolution%plink(this) = this
          resolution%plink_r(this) = this


          ! rmeove this from the ray link list
          next = resolution%rlink(this)
          resolution%rhead(rayn) = next

          resolution%good(good_tail) = this
          resolution%par_in_good(c_pindx) = .True.
          resolution%good_pindx(good_tail) = c_pindx
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
    deallocate(resolution%rhead)
    deallocate(resolution%rlink)
    deallocate(resolution%plink_r)
    deallocate(resolution%plink)
    deallocate(resolution%encoded)
    deallocate(resolution%par_in_good)
    resolution%secret = 0
  endsubroutine kill_resolution
end module
