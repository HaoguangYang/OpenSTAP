module BandwidthOpt

private
    type node
        integer :: index
        type(node), pointer :: next => null()
    end type
    
    type list
        integer :: length = 0
        type(node), pointer :: head => null()
        type(node), pointer :: tail => null()
    end type
    
    ! ���û���ظ�ֵ���ͼ���
    subroutine search(queue, node0)
    subroutine add(queue, node0)
    type(list), intent(out):: queue
    type(node):: node0
    temp => queue%head
    if(temp .EQ. null()) then
        queue%head
    do while(temp .NE. null())
        if(temp%index .EQ. p_node%index) then
            return
        end if
        temp = temp => next
    end do
    list