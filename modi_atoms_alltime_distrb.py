# -*- coding: utf-8 -*-
"""
Created on Thu Sep 29 16:48:03 2022
@author: SI H
"""

import numpy as np
import MDAnalysis as mda
import pandas as pd

def distance(p1, p2):
    res_dis = np.sqrt(np.sum(np.square(p1 - p2)))
    return res_dis 

frame_number = 6001
t_inter = 50
row_number = (frame_number-1)/t_inter + 1
num_tpr = 1  #tpr의 마지막 number (갯수 x)

def distribution(universe):  
    row = 0
    timestep = 0
    result = np.zeros((int(row_number), 7))
    
    for ts in universe.trajectory[0:frame_number:t_inter]:  
        cutoff_res = 4
        coor_B = []
        coor_C = []
        coor_A = []
        
#############################################################
        #atom name은 itp file이나 pdb file 참고해 기입# 
        
        #ATOM A :
        for_UNK = universe.select_atoms("name C0B")
        po = for_UNK.positions
        for i in range(200):
            coor_A.append(list(po[i]))
        
        #ATOM B :
        for_UNK = universe.select_atoms("name C01")
        po = for_UNK.positions
        for i in range(200):
            coor_B.append(list(po[i]))
        
        #ATOM C :
        for_UNK = universe.select_atoms("name C0A")
        po = for_UNK.positions
        for i in range(200):
            coor_C.append(list(po[i]))
            
#############################################################
        
        B_B_d = []
        B_C_d = []
        B_A_d = []
        C_C_d = []
        C_A_d = []
        A_A_d = []     
        
        A_A_distrb = 0
        B_A_distrb = 0
        C_A_distrb = 0
        B_B_distrb = 0
        B_C_distrb = 0
        C_C_distrb = 0        

        for i in range(len(coor_B)):    
            point_1 = np.array(coor_B[i])
            for j in range(i+1, len(coor_B)): #B-B
                point_2 = np.array(coor_B[j])
                res_dis = distance(point_1, point_2)
                if 0 < res_dis <= cutoff_res:
                    B_B_distrb += 1
                    B_B_d.append(res_dis)
            for j in range(i+1, len(coor_C)): #B-C
                point_2 = np.array(coor_C[j])
                res_dis = distance(point_1, point_2)
                if 0 < res_dis <=cutoff_res:
                    B_C_distrb += 1
                    B_C_d.append(res_dis)
            for j in range(i+1, len(coor_A)): #B-A
                point_2 = np.array(coor_A[j])
                res_dis = distance(point_1, point_2)
                if 0 < res_dis <=cutoff_res:
                    B_A_distrb += 1
                    B_A_d.append(res_dis)
                    
        for i in range(len(coor_C)):    
            point_1 = np.array(coor_C[i])
            for j in range(i+1, len(coor_C)): #C-C
                point_2 = np.array(coor_C[j])
                res_dis = distance(point_1, point_2)
                if 0 < res_dis <= cutoff_res:
                    C_C_distrb += 1
                    C_C_d.append(res_dis)     
            for j in range(i+1, len(coor_A)): #C-A
                point_2 = np.array(coor_A[j])
                res_dis = distance(point_1, point_2)
                if 0 < res_dis <= cutoff_res:
                    C_A_distrb += 1
                    C_A_d.append(res_dis)
                    
        for i in range(len(coor_A)):
            point_1 = np.array(coor_A[i])
            for j in range(i+1, len(coor_A)): #A-A
                point_2 = np.array(coor_A[j])
                res_dis = distance(point_1, point_2)
                if 0 < res_dis <= cutoff_res:
                    A_A_distrb += 1
                    A_A_d.append(res_dis)
                
        print(universe.trajectory.time)
        #print("A_A =", A_A_distrb)
        #print("A_B =", B_A_distrb)
        #print("A_C =", C_A_distrb)
        #print("B_B =", B_B_distrb)
        #print("B_C =", B_C_distrb)
        #print("C_C =", C_C_distrb)
        
        ##md_n 데이터 기록
        result[row][0] = timestep
        result[row][1] = A_A_distrb
        result[row][2] = B_A_distrb
        result[row][3] = C_A_distrb
        result[row][4] = B_B_distrb
        result[row][5] = B_C_distrb
        result[row][6] = C_C_distrb

        row+=1
        timestep += t_inter
        
    return result


######################################################################################################

final1 = np.zeros((int(row_number), num_tpr+1))
final2 = np.zeros((int(row_number), num_tpr+1))
final3 = np.zeros((int(row_number), num_tpr+1))
final4 = np.zeros((int(row_number), num_tpr+1))
final5 = np.zeros((int(row_number), num_tpr+1))
final6 = np.zeros((int(row_number), num_tpr+1))
final_result = np.zeros((int(row_number), 7))
timestamp = np.zeros((int(row_number)))

##분석 시작
for i in range(0, num_tpr+1):
    result = distribution(mda.Universe(f"md_{i}" + '.tpr', f"300ddvp_md_{i}" + '.xtc'))
    for j in range(int(row_number)):
        timestamp[j] = result[j][0] 
        final1[j][i] = result[j][1]
        final2[j][i] = result[j][2]
        final3[j][i] = result[j][3]
        final4[j][i] = result[j][4]
        final5[j][i] = result[j][5]
        final6[j][i] = result[j][6]
    print("finish : md", i)


avg1 = final1.mean(axis=1)
avg2 = final2.mean(axis=1)
avg3 = final3.mean(axis=1)
avg4 = final4.mean(axis=1)
avg5 = final5.mean(axis=1)
avg6 = final6.mean(axis=1)

for i in range (int(row_number)):
    final_result[i][0] = timestamp[i]
    final_result[i][1] = avg1[i]
    final_result[i][2] = avg2[i]
    final_result[i][3] = avg3[i]
    final_result[i][4] = avg4[i]
    final_result[i][5] = avg5[i]
    final_result[i][6] = avg6[i]

last = final_result[100:,1:].mean(axis=0)


##최종 result 내보내기
df = pd.DataFrame(final_result)
df.to_csv('all_of_atom_distrb.csv', index=False, header=False)

df.to_csv('AA_t.csv', columns = [0, 1], index=False, header=False)
df.to_csv('AB_t.csv', columns = [0, 2], index=False, header=False)
df.to_csv('AC_t.csv', columns = [0, 3], index=False, header=False)
df.to_csv('BB_t.csv', columns = [0, 4], index=False, header=False)
df.to_csv('BC_t.csv', columns = [0, 5], index=False, header=False)
df.to_csv('CC_t.csv', columns = [0, 6], index=False, header=False)

df = pd.DataFrame(last)
df.to_csv('last_1ns_distrb.csv', header=False)



