SUBROUTINE SOLVERMODE(ID)
    USE LIST_CLASS
    USE NODE_CLASS
    USE GLOBALS, ONLY : IIN, NEQ, NUMNP
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
     temp_length = LENGTH(LISTS(1))
     TEMP_INDEX = 1
     DO I = 1,NEQ
         IF(temp_length > length(lists(I))) THEN
            temp_length = length(lists(I))
            TEM_INDEX = I
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
                 if(temp_length > length(lists(i))) then
                    temp_length = length(lists(i))
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
END SUBROUTINE SOLVERMODE