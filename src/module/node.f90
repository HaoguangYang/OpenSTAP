! ����ʵ�ֵ�node����ֻ����һ��index,����û����
module node_class
private
    type node
        integer :: index_
        type(node), pointer :: next_ => null()
        type(node), pointer :: prev_ => null()
    end type node
    
    
    interface new
        module procedure new_this
    end interface
    
    interface operator( .EQ. )
        module procedure equal_this
    end interface
    
public node, new, operator(.EQ.)
    contains
    ! �½�ָ��
    subroutine new_this(this, index)
    type(node), pointer:: this
    integer :: index
    integer(4) :: err
        allocate(this, stat = err)
        this%index_  = index
    end subroutine new_this

    ! �����½�operator .EQ.
    logical function equal_this(p_node1, p_node2)
    type(node), intent(in) :: p_node1, p_node2
        equal_this = p_node1%index_ .EQ. p_node2%index_
    end function equal_this
end module node_class