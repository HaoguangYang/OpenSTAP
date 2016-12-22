module vector_class
    use node_class
    implicit none
private
    
    type vector
        integer :: index_
        integer :: length_ = 1
        integer :: length_real_ = 0
        integer, allocatable:: nodes_(:)
    end type vector

    interface set
        module procedure set_this
    end interface

    interface delete
        module procedure delete_this
    end interface
    
    interface add
        module procedure add_with_sort_this
    end interface
    
    public vector, add, set, delete

    contains
    subroutine set_this(this, index)
    type(vector) :: this
    integer :: index
    this%index_ = index
    allocate(this%nodes_(this%length_))
    end subroutine
    
    subroutine delete_this(this)
    type(vector) :: this
    deallocate(this%nodes_)
    end subroutine    
    
    ! add with out search
    subroutine add_with_sort_this(this, node0)
    type(vector):: this
    integer:: node0
    integer :: i,j
    if(this%index_ .gt. node0) return
    
    if( this%length_real_ .EQ. 0) then
        this%nodes_(1) = node0
        this%length_real_ = this%length_real_ + 1
    else
        do i = 1,this%length_real_
            if(node0 .lt. this%nodes_(i)) then
                do j = this%length_real_,i,-1
                    this%nodes_(j+1) = this%nodes_(j)
                end do
                this%nodes_(i) = node0
                this%length_real_ = this%length_real_+1
                return
            else if(node0 .eq. this%nodes_(i)) then
                return
            end if
        end do
        this%length_real_ = this%length_real_+1
        this%nodes_(this%length_real_) = node0
    end if
    end subroutine add_with_sort_this
    

end module vector_class