'''
···································
·                                 ·
·          INP文件读取            ·
·                                 ·
·          Qi He                  ·
·                                 ·
·          Tsinghua University    ·
·                                 ·
·          2016                   ·
···································
'''



#----------输入文件位置------------
import sys,os
path = sys.path[0]
filename = input('Please input your filename:\n')
location = path[2:].replace('\\','/') + '/' + filename + '.inp'
#location = '/Users/Administrator/Desktop/INP/8H_8.inp'
read = open(location,'r')


#----------定义变量----------------
#临时变量
i = 0
j = 0
k = 0
l = 0
#读取文件
readline = 1    #读入变量（成功读入则继续while语句）
line = []       #文件内容


#-----------读入文件到line-----------
while readline:
    readline=read.readline()
    line.append(readline[:-1])  #删除最后的回车符号
    i = i+1
read.close
line_len = len(line)            #INP文件全长度


#-----------读取特征关键词位置-------

#关键词位置信息
keylocation = [0 for x in range(line_len)] 

for i in range(0,line_len-1):
    if line[i].find('*Part')!=-1:           #'Part'--1
        keylocation[i] = 1
    if line[i].find('*Node')!=-1:           #'Node'--2
        keylocation[i] = 2
    if line[i].find('*Element')!=-1:        #'Element'--3
        keylocation[i] = 3
    #继续添加关键字---------
#print (keylocation)

j = 0
for i in range(0,line_len-1):
    if keylocation[i] == 3:
        j = 1
        continue
    if j == 1:
        if line[i][:1] == '*':
            keylocation[i-1] = 4            #每个element结束位置--4
            j = 0
        

keynum = max(keylocation)    #关键字个数
keyloc = [[] for x in range(keynum)]          
for i in range(0,line_len-1):
    for j in range(0,keynum):
        if keylocation[i] == (j+1):
            keyloc[j].append(i)
            
partnum = len(keyloc[0])     #PART个数
#print(keyloc)    
    
        
#-----------读取NODE信息-------------
def readnode(beg,end):
    
    #节点信息
    node_num = []   #NODE信息-node序号
    node_xyz = []   #NODE信息-节点xyz坐标
    #临时变量
    i = 0
    j = 0
    
    for i in range(beg,end+1):
    
        nodeline = line[i][line[i].find(',')+1:]        #'          -5.,          -5.,          10.' in 8H
        node_num.append(j+1)
        node_xyz.append([])                             #为了能添加内容所初始化
        
        while nodeline.find(',') != -1 :
            loc = nodeline.find(',')                    #找出','的位置
            node_xyz[j].append(float(nodeline[:loc]))
            nodeline = nodeline[loc+1:]
        node_xyz[j].append(float(nodeline[:]))
        j = j + 1
    
    return node_num,node_xyz
    
    
#-----------读取ELEMENT信息-------------
def readelement(beg,end):
    
    #单元信息
    element_num = []    #ELEMENT信息-element序号
    element_name = []   #ELEMENT信息-element名称
    element_node = []   #ELEMENT信息-单元节点号
    #临时变量
    i = 0
    j = 0
    
    for i in range(beg,end+1):
    
        elementline = line[i][line[i].find(',')+1:]    #' 10, 11, 14, 13,  1,  2,  5,  4' in 8H
        element_num.append(j+1)
        element_node.append([])                        #为了能添加内容所初始化
        
        while elementline.find(',') != -1 :
            loc = elementline.find(',')                #找出','的位置
            element_node[j].append(int(elementline[:loc]))
            elementline = elementline[loc+1:]
        element_node[j].append(int(elementline[:]))
        j = j + 1
    
    return element_num,element_node
    

#-----------读取PART名称-------------
def readpartname(beg):
    
    partname = line[beg][line[beg].find('name')+5:]
    return partname


#-----------读取ELEMENT名称----------
def readelementname(beg):
    
    elementname = line[beg][line[beg].find('type')+5:]
    return elementname
    
    
#-----------读取PART信息-------------
#部件信息
part = [[] for x in range(partnum)]       #包含对于不同部件的所有信息

for i in range(0,partnum):
    part[i].append(readpartname(keyloc[0][i]))
    part[i].append(readnode(keyloc[1][i]+1,keyloc[2][i]-1))
    part[i].append(readelementname(keyloc[2][i]))
    part[i].append(readelement(keyloc[2][i]+1,keyloc[3][i]))

    

#-----------输出PART信息-------------
for i in range(len(part)):
    print('···········第' + str(i+1) + '单元···········')
    print('-------------部件名称---------------')
    print(part[i][0])
    print('\n')     #换行
    print('-------------节点编号---------------')
    print(part[i][1][0])
    print('\n')
    print('-------------节点坐标---------------')
    print(part[i][1][1])
    print('\n')
    print('-------------单元类型---------------')
    print(part[i][2])
    print('\n')
    print('-------------单元编号---------------')
    print(part[i][3][0])
    print('\n')
    print('-------------单元节点---------------')
    print(part[i][3][1])
    print('\n')
    
    print('·····························')
    print('\n'*5)
    
    
    
    
    
    
    
    
    
    
    
    
    
    
