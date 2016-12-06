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
#location = '/Users/Administrator/Desktop/INP/4Q_4.inp'
read = open(location,'r')


#----------定义变量----------------
#临时变量
i = 0
j = 0
#读取文件
readline = 1    #读入变量（成功读入则继续while语句）
line = []       #文件内容
#节点信息
node_num = []   #NODE信息-node序号
node_xyz = []   #NODE信息-节点xyz坐标
#单元信息
element_num = []    #ELEMENT信息-element序号
element_name = []   #ELEMENT信息-element名称
element_node = []   #ELEMENT信息-单元节点号
#部件信息
part = []       #包含对于不同部件的所有信息


#-----------读入文件到line-----------
while readline:
    readline=read.readline()
    line.append(readline[:-1])  #删除最后的回车符号
    i = i+1
read.close
line_len = len(line)            #INP文件全长度


#-----------读取特征关键词位置-------

#关键词位置信息
keyloc = [0 for x in range(line_len)] 

for j in range(0,line_len-1):
    if line[j].find('*Part')!=-1:
        keyloc[j] = 1
    if line[j].find('*Node')!=-1:
        keyloc[j] = 2
    if line[j].find('*Element')!=-1:
        keyloc[j] = 3
    #继续添加关键字---------
print (keyloc)
    
    

'''
#-----------读取PART信息-------------
i = 0
k = 0
l = 0
for j in range(0,line_len-1):
    if line[j].find('*Part')!=-1:
        #l = l + 1
        part.append([])
        part[l].append(line[j][line[j].find('name')+5:])
'''

        
#-----------读取NODE信息-------------
i = 0
k = 0
for j in range(0,line_len-1):
    if line[j].find('*Node')!=-1:
        i = j + 1
        
        while line[i][0] == ' ':
            nodeline = line[i][line[i].find(',')+1:]    #'          -5.,          -5.,          10.' in 8H
            node_num.append(k+1)
            node_xyz.append([])                         #为了能添加内容所初始化
            
            while nodeline.find(',') != -1 :
                loc = nodeline.find(',')                #找出','的位置
                node_xyz[k].append(float(nodeline[:loc]))
                nodeline = nodeline[loc+1:]
            node_xyz[k].append(float(nodeline[:]))
            
            i = i + 1
            k = k + 1

print('-----节点对应序号-----')            
print(node_num)
print('-----节点坐标-----')   
print(node_xyz)

#-----------读取ELEMENT信息-------------
i = 0
k = 0
for j in range(0,line_len-1):
    if line[j].find('*Element')!=-1:
        i = j + 1
        element_name = line[j][line[j].find('type')+5:] #读取单元名称
        
        while line[i][0] != '*':
            elementline = line[i][line[i].find(',')+1:]    #' 10, 11, 14, 13,  1,  2,  5,  4' in 8H
            element_num.append(k+1)
            element_node.append([])                         #为了能添加内容所初始化
            
            while elementline.find(',') != -1 :
                loc = elementline.find(',')                #找出','的位置
                element_node[k].append(int(elementline[:loc]))
                elementline = elementline[loc+1:]
            element_node[k].append(int(elementline[:]))
            
            i = i + 1
            k = k + 1 

print('-----单元类型-----')               
print(element_name) 
print('-----单元序号-----')             
print(element_num)
print('-----单元节点号----')  
print(element_node)
