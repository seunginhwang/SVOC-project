# -9- coding: utf-8 -*-
"""
Created on Thu Sep  1 15:52:55 2022

@author: iiz20
"""

import numpy as np
import MDAnalysis as mda
import pandas as pd


nmol = 10    #num of SOL
nmol1 = 100   #num of HQN
nmol2 = 100   #num of DDVP

frame_number = 6001
row_number = (frame_number-1)/50 + 1

num_tpr = 4

def density(universe): 
    z_arr = list(range(nmol+nmol1+nmol2))
    z_arr0 = list(range(nmol))
    z_arr1 = list(range(nmol1))
    z_arr2 = list(range(nmol2))
    SOL = universe.select_atoms("resname SOL")
    HQN = universe.select_atoms("resname HQN")
    UNK = universe.select_atoms("resname UNK")
    time_all = np.zeros((34,21))
    time_SOL = np.zeros((34,21))
    time_HQN = np.zeros((34,21))
    time_DDVP = np.zeros((34,21))
    t_inter = 50
    column = 0    
    for ts in universe.trajectory[5000:frame_number+1:t_inter]:
        result_all = [0 for i in range(34)]
        result_SOL = [0 for i in range(34)]
        result_HQN = [0 for i in range(34)]
        result_DDVP = [0 for i in range(34)]
        
        co = [0]
        pos = []
        row_index = list(range(34))
        
        for i in range(nmol):           
            stt_n = 3*i
            fin_n = 3*(i+1)-1
            for_SOL = SOL[stt_n:fin_n]
            
            co = for_SOL.center_of_mass()
            pos.append(list(co))
            
            z = co[2]
            z_arr0[i] = z
            z_arr[i] = z
            # z_arr = UNK209~UNK408 z좌표만 받아놓은 리스트 0~199
            
        for i in range(nmol1):
            stt_n1 = 14*i
            fin_n1 = 14*(i+1)-1
            for_HQN = HQN[stt_n1:fin_n1]
            
            co1 = for_HQN.center_of_mass()
            pos.append(list(co1))

            z1 = co1[2]
            z_arr1[i] = z1
            z_arr[i+nmol] = z1
            
        for i in range(nmol2):
            stt_n2 = 18*i
            fin_n2 = 18*(i+1)-1
            for_UNK = UNK[stt_n2:fin_n2]
            
            co2 = for_UNK.center_of_mass()
            pos.append(list(co2))

            z2 = co2[2]
            z_arr2[i] = z2
            z_arr[i+nmol+nmol1] = z2        
        
             
        for i in range(0, 34): #0~17
            z_dist = i
            row_index[i] = 5*z_dist
        
        '''
        전체        z_arr = list(range(nmol+nmol1+nmol2))
        물          z_arr0 = list(range(nmol))
        하이드로퀴논  z_arr1 = list(range(nmol1))
        디클로르보스  z_arr2 = list(range(nmol2))  
        '''
        
        for j in range(nmol+nmol1+nmol2):
            num = z_arr[j]
            for k in range(0,34):
                z_dist = k
                z_index = z_dist*5
                if z_dist*5 < num <= z_dist*5 + 5:
                    index = row_index.index(z_index)
                    result_all[index] += 1
                    break
            
        for j in range(nmol):
            num = z_arr0[j]
            for k in range(0,34):
                z_dist = k
                z_index = z_dist*5
                if z_dist*5 < num <= z_dist*5 + 5:
                    index = row_index.index(z_index)
                    result_SOL[index] += 1
                    break
            
        for j in range(nmol1):
            num = z_arr1[j]
            for k in range(0,34):
                z_dist = k
                z_index = z_dist*5
                if z_dist*5 < num <= z_dist*5 + 5:
                    index = row_index.index(z_index)
                    result_HQN[index] += 1
                    break
            
        for j in range(nmol2):
            num = z_arr2[j]
            for k in range(0,34):
                z_dist = k
                z_index = z_dist*5
                if z_dist*5 < num <= z_dist*5 + 5:
                    index = row_index.index(z_index)
                    result_DDVP[index] += 1
                    break
                
        print("sum =", sum(result_all), 
              "/ SOL_sum =", sum(result_SOL), 
              "/ HQN_sum =", sum(result_HQN),
              "/ DDVP_sum =", sum(result_DDVP))
 

        for i in range(34):
            time_all[i][column] = result_all[i]
            time_HQN[i][column] = result_HQN[i]
            time_DDVP[i][column] = result_DDVP[i]
            
            if nmol != 0:
                time_SOL[i][column] = result_SOL[i]
        
        column += 1
        print(universe.trajectory.time)
    
    final = np.zeros((34, 5))

    for i in range (0, 34):
        final[i][0] = i/2

    mean = time_all.mean(axis=1)
    mean0 = time_SOL.mean(axis=1)
    mean1 = time_HQN.mean(axis=1)
    mean2 = time_DDVP.mean(axis=1)


    for j in range(34):
        final[j][1] = mean[j]
        final[j][2] = mean0[j]   
        final[j][3] = mean1[j]
        final[j][4] = mean2[j]
        
    #0열 전체 평균 밀도
    #1열 SOL의 평균 밀도 
    #2열 HQN의 평균 밀도 
    #3열 DDVP의 평균 밀도 
    return final


#final 
#0열은 z축 좌표 목록
#1열은 좌표당 전체 분자 밀도
#2열은 좌표당 SOL 밀도
#3열은 좌표당 HQN 밀도
#4열은 좌표당 DDVP 밀도

final1 = np.zeros((34, num_tpr+1))
final2 = np.zeros((34, num_tpr+1))
final3 = np.zeros((34, num_tpr+1))
final4 = np.zeros((34, num_tpr+1))
final_result = np.zeros((34, 5))
zstamp = np.zeros((34))

##분석 시작
for i in range(0, num_tpr+1):
    result = density(mda.Universe(f"md_{i}" + '.tpr', f"280mix_md_{i}" + '.xtc'))
    for j in range(34):
        zstamp[j] = result[j][0] 
        final1[j][i] = result[j][1]
        final2[j][i] = result[j][2]
        final3[j][i] = result[j][3]
        final4[j][i] = result[j][4]
    print("finish : md", i)

avg1 = final1.mean(axis=1)
avg2 = final2.mean(axis=1)
avg3 = final3.mean(axis=1)
avg4 = final4.mean(axis=1)

for i in range (34):
    final_result[i][0] = zstamp[i]
    final_result[i][1] = avg1[i]
    final_result[i][2] = avg2[i]
    final_result[i][3] = avg3[i]
    final_result[i][4] = avg4[i]


df = pd.DataFrame(final_result)
df.to_csv('mix_density_avg.csv', index=False, header=False)

df.to_csv('all_den.csv', columns = [0, 1], index=False, header=False)
df.to_csv('sol_den.csv', columns = [0, 2], index=False, header=False)
df.to_csv('hqn_den.csv', columns = [0, 3], index=False, header=False)
df.to_csv('ddvp_den.csv', columns = [0, 4], index=False, header=False)


    














