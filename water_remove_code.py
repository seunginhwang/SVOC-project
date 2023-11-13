#-*- coding: utf-8 -*-
"""
Created on Fri Sep  2 15:44:08 2022

@author: iiz20
"""

#!/bin/python3

with open('sol_4000_cp.gro' , 'r') as lis:
    line = lis.read().splitlines()

#rem : 지울 z축 좌표 목록
#axis : 띄쓰마다 split
#line : 줄마다 split
 
line0 = line[0]
line1 = line[1]
rem = []

del line[0]
del line[0]

line_num = len(line) #21246 (index 0~21245)

n = 14 #atom number of molecule
nmol = 100 #number of molecule 


####################################################################


for i in range(2, line_num): #2 ~ line_num-1
    line_split = line[i].split()
    line_splt_num = len(line_split) - 1
    
    if float(line_split[line_splt_num]) <= 2.5:
        rem.append(line_split[line_splt_num])



rem = set(rem)
rem_sort = sorted(rem)
rem = list(rem_sort)

rem_num = len(rem)  #7209 (index 0~7208)   

#rem : 지울 z축 좌표 목록
#axis : 띄쓰마다 split
#line : 줄마다 split 

#####################################################################

line_list = list(line)
del_count = 0

for i in range( 7844 + n*nmol + 1 , line_num): #index 9245 ~ 21245
    split_line = line_list[i].split()
    sp_num_index = len(split_line)-1
    for j in range(rem_num): #index 0~7208  
        if rem[j] == split_line[sp_num_index]:
            line_list[i] = "removed"
            del_count += 1
        
print('number of removed water molecules = ', del_count)
line1 = int(line1) - 3*del_count

begin = []
begin.append(line0)
begin.append(line1)

while 'removed' in line_list:
    line_list.remove('removed') 

fin = begin + line_list
fin_num = len(fin)

fin[0] = str(fin[0])
fin[1] = str(fin[1])


for i in range(fin_num):
    fin[i] = fin[i] + '\n'


with open('removed.gro', 'w') as f:
    for line in fin:
        f.write(line)

print('number of final water molecules = ', (line_num-1 - (7845 + n*nmol))/3 - del_count)



















