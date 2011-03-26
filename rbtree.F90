!< Note a red black tree yet. currently only a binary tree is implemented
!< I strongly feel the binary tree will degrade to a link list very quickly.
!< but a r/b tree is too complicated to implement.

module  rbtree_mod
use myf03_mod
integer(i8b), parameter :: MAX_LENGTH = 1000
type rbtree_type
   integer(i8b), allocatable :: value(:)  !< the task ids
   integer(i8b), allocatable :: key(:)  !< the hungriness of the tasks
   integer(i1b), allocatable :: color(:) !< preserved for r/b tree
   integer(i8b), allocatable :: left(:)  !< left child
   integer(i8b), allocatable :: right(:) !< right child
   integer(i8b), allocatable :: parent(:) !< parent
   integer(i8b) :: leftmost     !< leftmost element, the one to be ejected
   integer(i8b) :: ejected      !< the ejected task, to be returned by rbtree_return
   integer(i8b) :: length       !< number of used storage elements
   integer(i8b) :: root         !< root of the tree
endtype rbtree_type

contains
  
  subroutine rbtree_prepare(rbtree, max_length)
    type(rbtree_type), intent(inout) :: rbtree
    integer(i8b), intent(in) :: max_length
    allocate(rbtree%value(max_length))
    allocate(rbtree%key(max_length))
    allocate(rbtree%color(max_length))
    allocate(rbtree%left(max_length))
    allocate(rbtree%right(max_length))
    allocate(rbtree%parent(max_length))
    rbtree%value = 0
    rbtree%color = 0
    rbtree%left = 0
    rbtree%right = 0
    rbtree%left = 0
    rbtree%leftmost = 0
    rbtree%root = 0
  endsubroutine rbtree_prepare

  subroutine rbtree_insert(rbtree, key, value)
    type(rbtree_type), intent(inout) :: rbtree
    integer(i8b), intent(in) :: key, value
    call rbtree_insert0(rbtree, key, value, 0)
  endsubroutine rbtree_insert

  subroutine rbtree_insert0(rbtree, key, value, id)
    type(rbtree_type), intent(inout) :: rbtree
    integer(i8b), intent(in) :: key, value
    integer(i8b), intent(in) :: id
    integer(i8b) :: parent, this, new

     if(id == 0) then
       new = rbtree%length + 1
       rbtree%length = new
     else 
       new = id
     endif

     rbtree%value(new) = value
     rbtree%key(new) = key
     rbtree%left(new) = 0
     rbtree%right(new) = 0
     rbtree%color(new) = 1

     if (rbtree%root == 0) then
       rbtree%root = new
       rbtree%leftmost = new
       return
     endif

     this = rbtree%root

     do while(this /= 0)
       parent = this
       if ( key < rbtree%key(this) ) then
          this = rbtree%left(this)
          if (this == 0) then
             rbtree%left(parent) = new
             rbtree%parent(new) = parent
             if(rbtree%leftmost == parent) then
               rbtree%leftmost = new
             endif
             exit
          endif
       else
          this = rbtree%right(this)
          if (this == 0) then
             rbtree%right(parent) = new
             rbtree%parent(new) = parent
             exit
          endif
       endif
     enddo
     
  endsubroutine rbtree_insert0


  subroutine rbtree_print(rbtree)
    type(rbtree_type), intent(inout) :: rbtree
    integer(i8b) :: parent, this, new
    if (rbtree%length == 0) then
      print *, 'empty tree'
      return
    endif
    print *, 'leftmostid', rbtree%leftmost, 'rootid', rbtree%root
    this = rbtree%root
    do while(this /= 0)
       if(rbtree%left(this) == 0) then
         print *, 'node', this, rbtree%key(this)
         this = rbtree%right(this)
         cycle
       endif
       next = rbtree%left(this)
       do while(rbtree%right(next) /= 0 .and. rbtree%right(next) /= this)
         next = rbtree%right(next)
       enddo
       if(rbtree%right(next) == 0) then
          rbtree%right(next) = this
          this = rbtree%left(this)
          cycle
       else
          rbtree%right(next) = 0
         print *, 'node', this, rbtree%key(this)
          this = rbtree%right(this)
          cycle
       endif
    enddo
  endsubroutine

  subroutine rbtree_eject(rbtree, key, value)
    !< select the most hungry task in the tree,
    !< returns the hungriness in 'key', and the value stored in value

    type(rbtree_type), intent(inout) :: rbtree
    integer(i8b),intent(out) :: key
    integer(i8b),intent(out) :: value
    integer(i8b) :: parent, ejected, new
    ejected = rbtree%leftmost
    print *, 'ejected', ejected
    if(ejected == rbtree%root) then
      new = rbtree%right(ejected)
      rbtree%root = new
      rbtree%leftmost = new
      rbtree%parent(new) = 0
    else
      parent = rbtree%parent(ejected)
      new = rbtree%right(ejected)
      rbtree%left(parent) = new
      rbtree%parent(new) = parent
      if(new == 0) then
       ! if already a leaf ejected, the parent is the new leftmost
        rbtree%leftmost = parent
      else 
        rbtree%leftmost = new
      endif
    endif
    rbtree%left(ejected) = 0
    rbtree%right(ejected) = 0
    rbtree%parent(ejected) = 0
    rbtree%color(ejected) = 0
    rbtree%ejected = ejected
    key = rbtree%key(ejected)
    value = rbtree%value(ejected)
  endsubroutine rbtree_eject

  subroutine rbtree_return(rbtree, key, value)
    !< return the returned hungry task from rbtree_eject back 
    !< with a new hungriness in 'key', and new value stored in value
    type(rbtree_type), intent(inout) :: rbtree
    integer(i8b),intent(in) :: key
    integer(i8b),intent(in) :: value
    call rbtree_insert0(rbtree, key, value, rbtree%ejected)
  endsubroutine rbtree_return

  subroutine rbtree_kill(rbtree)
    type(rbtree_type), intent(inout) :: rbtree
    deallocate(rbtree%value)
    deallocate(rbtree%key)
    deallocate(rbtree%color)
    deallocate(rbtree%left)
    deallocate(rbtree%right)
    deallocate(rbtree%parent)
  endsubroutine


endmodule rbtree_mod

program rbtree_test
use rbtree_mod
  type(rbtree_type) :: rbtree
  integer(i8b):: key, value
  integer(i4b):: i
  call rbtree_prepare(rbtree, 1000)
  call rbtree_insert(rbtree, 5, 1)
  call rbtree_insert(rbtree, 1, 2)
  call rbtree_insert(rbtree, 0, 3)
  call rbtree_insert(rbtree, 2, 4)
  call rbtree_insert(rbtree, 3, 5)
  call rbtree_print(rbtree)
  do i = 1, 5
  call rbtree_eject(rbtree, key, value)
  call rbtree_print(rbtree)
  call rbtree_return(rbtree, key+ 10, value)
  call rbtree_print(rbtree)
  enddo
end program
