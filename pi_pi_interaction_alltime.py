# -*- coding: utf-8 -*-
"""
Created on Thu Sep 29 16:48:03 2022
@author: SI H
"""
frame_number = 6001
t_inter = 50
num_mol = 200
atom = 14
row_number = (frame_number-1)/t_inter + 1

import numpy as np
import MDAnalysis as mda
import pandas as pd
from math import pi

def distance(p1, p2):
    res_dis = np.sqrt(np.sum(np.square(p1 - p2)))
    return res_dis 

def normal_vector(p1, p2, p3):
    l1 = p2-p1
    l2 = p3-p1
    ans = np.cross(l1,l2)
    return ans

def angle(v1, v2):
    top = np.dot(v1,v2)
    bot = np.sqrt(v1[0]**2 + v1[1]**2 + v1[2]**2) * np.sqrt(v2[0]**2 + v2[1]**2 + v2[2]**2)
    cos = top/bot
    ans = np.arccos(cos) * (180/pi)
    if ans > 90:
        ans = ans- 90
    return ans

def distribution(universe):  
    row = 0
    timestep = 0
    result = np.zeros((int(row_number), 7))
    
    for ts in universe.trajectory[0:frame_number:t_inter]:  
        cutoff_res = 5
        cutoff_angle = 5
        count = 0
        
        coor_o1 = []
        coor_o2 = []
        coor_c = []
        normal_v = []
        dihedral = []
        pos = []
        dist = []
        
        UNK = universe.select_atoms("resname UNK")
        for i in range(num_mol):
            stt_n = atom*i
            fin_n = atom*(i+1)-1
            for_UNK = UNK[stt_n:fin_n]
            co = for_UNK.center_of_mass()
            pos.append(list(co))
        
        for_UNK1 = universe.select_atoms("name O01")
        po1 = for_UNK1.positions
        for i in range(200):
            coor_o1.append(list(po1[i]))
                        
        for_UNK2 = universe.select_atoms("name O0C")
        po2 = for_UNK2.positions
        for i in range(200):
            coor_o2.append(list(po2[i]))
            
        for_UNK3 = universe.select_atoms("name C02")
        po3 = for_UNK3.positions
        for i in range(200):
            coor_c.append(list(po3[i]))

        for j in range(200):
            vector_1 = np.array(coor_o1[j])
            vector_2 = np.array(coor_o2[j])
            vector_3 = np.array(coor_c[j])
            nor = normal_vector(vector_1, vector_2, vector_3)
            normal_v.append(nor)
                
        di = []
        for i in range(200):
            for j in range(i+1, 200):
                dihedral = angle(normal_v[i], normal_v[j])
                di.append(dihedral)
                point_1 = np.array(pos[i])
                point_2 = np.array(pos[j])
                dis = distance(point_1, point_2)
                dist.append(dis)
                if dihedral <= cutoff_angle and dis < cutoff_res:
                    print("angle =", round(dihedral,2),"°", "/ distance =", round(dis,2), "Å")
                    count += 1
                else:
                    pass

        ##md_n 데이터 기록
        print("time =", timestep,"count =",count)
        result[row][0] = timestep
        result[row][1] = count
        print("\n")
        row+=1
        timestep += t_inter
        
    return result


##md_0 분석 시작
u0 = mda.Universe("md_0.tpr", "300hqn_md_0.xtc")
result0 = distribution(u0)

print("analysis finish")

last = result0[100:,1:].mean(axis=0)
lst_avg = np.zeros((1,1))

print("lst_avg =", last[0])
    
##최종 result 내보내기
df = pd.DataFrame(result0)
df.to_csv('pipi_alltime_cutoff.10.csv', index=False, header=False)

df = pd.DataFrame(last)
df.to_csv('last_lns_pipi_cutoff.10.csv', index=False, header=False)

    
    
    
    
    
    
    
    
    
    
    
    
    
    
    