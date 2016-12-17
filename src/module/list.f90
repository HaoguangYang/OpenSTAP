module list_class
    use node_class
    implicit none
private
    
    type list
        integer :: row_    = 0
        integer :: column_ = 0
        integer :: length_ = 0
        integer :: sign_ = 0
        type(node), pointer :: head_ => null()
        type(node), pointer :: tail_ => null()
    end type list
    
    interface set
        module procedure set_this
    end interface
    
    interface length
        module procedure length_this
    end interface
    
    interface search
        module procedure search_this
    end interface
    
    interface add
        module procedure add_this
    end interface
    
    interface delete
        module procedure delete_this
    end interface
    
    public node, list, length, add, search, add_with_search, delete, set

    contains
    !new a node
    subroutine new_this(this, index)
        type(node) , pointer:: this
        integer    :: index
        integer(4) :: err
        allocate(this, stat = err)
        this%index_ = index
    end subroutine new_this
    
    ! set the row and colum for the list
    subroutine set_this(this, row, column)
        type(list) :: this
        integer :: row
        integer :: column
        this%row_ = row
        this%column_ = column
    end subroutine set_this
    
    ! return length of the list
    integer function length_this(this)
        type(list) :: this
        length_this = this%length_
    end function length_this

    ! search for p_node0
    logical function search_this(this, p_node0)
        type(list), intent(in) :: this
        type(node), pointer    :: p_node0
        type(node), pointer    :: p_node
        search_this = .false.
        p_node => this%head_
        do while( associated(p_node) )
            if(p_node .EQ. p_node0) then
                search_this = .true.
                exit
            end if
            p_node => p_node%next_
        end do
    end  function search_this
    
    ! add with out search
    subroutine add_this(this, p_node0)
    type(list):: this
    type(node), pointer:: p_node0
        if(associated(this%head_)) then
            this%tail_%next_ => p_node0
            p_node0%prev_ => this%tail_
            this%tail_ => p_node0
        else
            this%head_ => p_node0
            this%tail_ => p_node0
        end if
        this%length_ = this%length_ + 1
    end subroutine add_this
    
    subroutine add_with_search(this, p_node0)
    type(list):: this
    type(node), pointer :: p_node0
        if(.not. search(this, p_node0)) then
            if(associated(this%head_)) then
                this%tail_%next_ => p_node0
                this%tail_ => p_node0
            else
                this%head_ => p_node0
                this%tail_ => p_node0
            end if
            this%length_ = this%length_ + 1
        end if
    end subroutine add_with_search
    
    subroutine delete_this(this, p_node0)
        type(list) :: this
        type(node), pointer    :: p_node0
        type(node), pointer    :: p_node
        integer :: err
        p_node => this%head_
        do while( associated(p_node) )
            if(p_node .EQ. p_node0) then ! 找到了要删除的节点
                if(.not. associated(p_node%prev_)) then
                    ! 删除头节点
                    if(associated(p_node%next_)) then
                        this%head_ => p_node%next_
                        nullify(p_node%next_%prev_)
                    else
                        ! 原链表中只有一个节点
                        nullify(this%head_)
                        nullify(this%tail_)
                    end if
                else ! 不是头节点
                    if(associated(p_node%next_)) then
                        p_node%next_%prev_ => p_node%prev_
                        p_node%prev_%next_ => p_node%next_
                    else ! 删除尾节点
                        nullify(p_node%prev_%next_)
                        this%tail_ => p_node%prev_
                    end if
                end if
                this%length_ = this%length_ - 1
                
                deallocate(p_node, stat = err)
                if(err .ne. 0) stop 'Fail to delete this node'
                exit
            end if
            p_node => p_node%next_
        end do
    end subroutine delete_this
end module list_class