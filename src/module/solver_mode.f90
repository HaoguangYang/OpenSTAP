subroutine bdopt(ID)
    USE LIST_CLASS
    USE NODE_CLASS
    USE GLOBALS, ONLY : IIN, IOUT, NEQ, NUMNP, NWK
    
    implicit none
    TYPE(LIST) :: LISTS(NEQ)
    TYPE(NODE), POINTER :: P_NODE
    INTEGER :: ID(6, NUMNP), ID_NEW
    INTEGER :: NUME ! 全部总共的节点个数
    INTEGER :: N    ! 节点中的函数
    INTEGER :: TEMP_NODE(8) ! 用来暂存
    INTEGER :: TEMP_ID(48) ! 用来暂存节点信息
    INTEGER :: I, J, K ! 循环变量
    INTEGER :: temp_length, TEMP_INDEX
    INTEGER :: FREE_DOF ! 单元中自由的自由度
    character (len = 10) :: cFmt
    DO I = 1,NUMNP
        DO J = 1,6
            IF(ID(J,I) .NE. 0) THEN
                CALL SET(lists(ID(J,I)), I, J)
            END IF
        END DO
    END DO
    ! input phase
     READ (IIN,"(I5)") NUME ! 总element数
     DO I = 1,NUME
         READ (IIN,"(I5)") N ! 这个element对应的节点数
         WRITE (IOUT, "('(',I5,'I5)')") N
         READ (IIN, "(8I5)") (TEMP_NODE(J), J = 1,N)
         FREE_DOF = 0
         DO J = 1,N
             DO K = 1,6
                 IF (ID(K, TEMP_NODE(J)) .NE. 0) THEN
                     FREE_DOF = FREE_DOF + 1
                     TEMP_ID(FREE_DOF) = ID(K, TEMP_NODE(J))
                 END IF
             END DO
         END DO
         DO J = 1, FREE_DOF
             DO K = 1, FREE_DOF
                 CALL NEW(P_NODE, TEMP_ID(J))
                 call ADD_WITH_SEARCH(lists(TEMP_ID(K)), P_NODE)
             END DO
         END DO
     END DO
        ! CM algorithm
        ! 先找到最小的节点
        ID_NEW = 1
        temp_length = LISTS(1)%length_
        TEMP_INDEX = 1
        DO I = 1,NEQ
             IF(temp_length > lists(I)%length_) THEN
                temp_length = lists(I)%length_
                TEMP_INDEX = I
            END IF
        END DO
        ID(lists(TEMP_INDEX)%column_, lists(TEMP_INDEX)%row_) = ID_NEW
        ID_NEW = ID_NEW + 1
        ! 后续的搜索循环
        p_node => lists(temp_index)%head_
        do while(associated(p_node))
            lists(p_node%index_)%sign_ = 1
            p_node => p_node%next_
        end do
        lists(temp_index)%sign_ = -1
        do while( id_new .le. neq)
            temp_index = 0
            temp_length = 100000 ! 应该要放最大整数，但是忘了是多少了
            do i = 1, neq
              if(lists(i)%sign_ .eq. 1) then
                  if(temp_length > lists(i)%length_) then
                    temp_length = lists(i)%length_
                    temp_index  = i
                  end if
             end if
           end do
           id(lists(temp_index)%column_, lists(temp_index)%row_) = id_new
           id_new = id_new+1
           p_node => lists(temp_index)%head_
           do while(associated(p_node))
              if(lists(p_node%index_)%sign_ .ne. -1) then
                  lists(p_node%index_)%sign_ = 1
              end if
              p_node => p_node%next_
           end do
           lists(temp_index)%sign_ = -1
        end do
        ! Write equation numbers
        WRITE (IOUT,"(//,' EQUATION NUMBERS',//,'   NODE',9X,  &
                     'DEGREES OF FREEDOM',/,'  NUMBER',/,  &
                     '     N',13X,'X    Y    Z   RX   RY   RZ',/,(1X,I5,9X,6I5))") (N,(ID(I,N),I=1,6),N=1,NUMNP)
     !  清空list，释放内存
     do i = 1, neq
         call delete_all(lists(i))
     end do
end subroutine bdopt
    
subroutine pardiso_input(ID)
    USE LIST_CLASS
    USE NODE_CLASS
    USE GLOBALS, ONLY : IIN, IOUT, NEQ, NUMNP, NWK, pardisodoor, TIM
    USE MEMALLOCATE
    
    implicit none
    TYPE(LIST) :: LISTS(NEQ)
    TYPE(NODE), POINTER :: P_NODE
    INTEGER :: ID(6, NUMNP), ID_NEW
    INTEGER :: NUME ! 全部总共的节点个数
    INTEGER :: N    ! 节点中的函数
    INTEGER :: TEMP_NODE(8) ! 用来暂存
    INTEGER :: TEMP_ID(48) ! 用来暂存节点信息
    INTEGER :: I, J, K ! 循环变量
    INTEGER :: temp_length, TEMP_INDEX
    INTEGER :: FREE_DOF ! 单元中自由的自由度
    ! input phase
     READ (IIN,"(I5)") NUME ! 总element数
     ! 为了节省空间，这里用sign_表示序号好了
     DO i = 1,NEQ
         lists(i)%sign_ = i
     end do
     DO I = 1,NUME
         READ (IIN,"(I5)") N ! 这个element对应的节点数
         READ (IIN, "(<N>I5)") (TEMP_NODE(J), J = 1,N)
         FREE_DOF = 0
         DO J = 1,N
             DO K = 1,6
                 IF (ID(K, TEMP_NODE(J)) .NE. 0) THEN
                     FREE_DOF = FREE_DOF + 1
                     TEMP_ID(FREE_DOF) = ID(K, TEMP_NODE(J))
                 END IF
             END DO
         END DO
         DO K = 1, FREE_DOF
             DO J = 1, FREE_DOF
                 CALL NEW(P_NODE, TEMP_ID(J))
                 if( (K.eq.1) .and. (I.eq.55) ) then
                     ID_new = K
                 end if
                 call ADD_WITH_SORT(lists(TEMP_ID(K)), P_NODE)
             END DO
         END DO
     END DO
     CALL SECOND (TIM(2))
     ! 注意这里分配了rowIndex
     CALL MEMALLOC(2,"rwInd",NEQ+1,1)                   !RowIndex
     CALL assign_rowIndex(lists, IA(NP(2)))
     CALL MEMALLOC(3,"STFF ",NWK,ITWO)
     CALL MEMALLOC(4,"R    ",NEQ,ITWO)
     CALL MEMALLOC(5,"colmn",NWK,1)
     CALL assign_columns(lists, IA(NP(5)))
     do i = 1, neq
         call delete_all(lists(i))
     end do
    end subroutine pardiso_input

subroutine assign_rowIndex(lists, rowIndex)
    use GLOBALS, only : neq, nwk
    use list_class
    use node_class
    
    implicit none
    type(list) :: lists(neq)
    integer :: rowIndex(neq+1)
    type(node), pointer :: p_node
    integer :: i, j
    
    nwk = 0
    do i=1, neq
        rowIndex(i) = 1 + nwk
        nwk = nwk + lists(i)%length_
    end do
    rowIndex(neq+1) = nwk+1
    end subroutine assign_rowIndex
    
subroutine assign_columns(lists, columns)
    use GLOBALS, only : neq, nwk
    use list_class
    use node_class
    
    implicit none
    type(list) :: lists(neq)
    integer :: columns(nwk)
    type(node), pointer :: p_node
    integer :: i, j
    
    j = 1
    do i = 1, neq
        p_node => lists(i)%head_
        do while(associated(p_node))
            columns(j) = p_node%index_
            p_node => p_node%next_
            j = j + 1
            if( j == nwk+1) then
                j = j-1
                exit
            end if
        end do
    end do
    end subroutine assign_columns
    