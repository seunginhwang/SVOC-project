import numpy as np
import MDAnalysis as mda
import pandas as pd

atom = 18
num_mol = 200

frame_number = 3001
t_inter = 50
num_tpr = 1  #tpr의 마지막 number (갯수 x)


def distance(p1, p2):
    res_dis = np.sqrt(np.sum(np.square(p1 - p2)))
    return res_dis 

def ratio(universe):
    UNK = universe.select_atoms("resname UNK")
    
    z_arr = list(range(num_mol))
    row = 0
    timestep=0
    result = np.zeros((int(row_number), 6))

    for ts in universe.trajectory[0:frame_number:t_inter]:  
        cutoff_res1 = 5.2 + 0.01*(timestep/10)
        cutoff_res2 = 6.2 + 0.005*((timestep-1000)/20)
        z_arr_cutoff = 40 + 0.375*(timestep-1000)/40
        
        co = [0]
        pos = []
        count = []
        for i in range(num_mol):
            stt_n = atom*i
            fin_n = atom*(i+1)-1
            for_UNK = UNK[stt_n:fin_n]
            co = for_UNK.center_of_mass()
            pos.append(list(co))

            z = co[2]
            z_arr[i] = z
            z_sort = sorted(z_arr)

        # z_arr = UNK209~UNK408 z좌표만 받아놓은 리스트 0~199
        # z_sort = UNK z 높이순 리스트

        #min_z = min(z_arr)
        adsorp_n = 1
        count.append(z_arr.index(z_sort[0]))
        count_case = [0,0,0,0,0]
        
        for i in range(1, num_mol): #1~199
            sel_index = z_arr.index(z_sort[i])
            z_dis = z_arr[sel_index] - 24.5
            #z_dis2 = z_arr[sel_index] - min_z
            selected = sel_index
            
            for j in range(0,i):
                sell_index = z_arr.index(z_sort[j])
                point_1 = np.array(pos[sell_index])
                point_2 = np.array(pos[sel_index])
                res_dis = distance(point_1, point_2)
            
                if z_dis<3.5:
                    adsorp_n += 1
                    count.append(selected)
                    count_case[0] += 1
                    break
                 
                elif timestep <= 1000 and res_dis < cutoff_res1:
                    adsorp_n += 1
                    count.append(selected)
                    count_case[1] += 1
                    break
                elif timestep <= 1000 and z_arr[sel_index]<40:
                    adsorp_n += 1
                    count.append(selected)
                    count_case[2] += 1
                    break
                elif timestep > 1000 and res_dis < cutoff_res2:
                    adsorp_n += 1
                    count.append(selected)
                    count_case[3] += 1
                    break
                elif timestep > 1000 and z_arr[sel_index]<z_arr_cutoff:
                    adsorp_n += 1
                    count.append(selected)
                    count_case[4] += 1
                    break

                
        #print(count_case)  
        #adsorp_p = adsorp_n / num_mol * 100
       
        print(timestep)
        print("len=",len(count),"adsorp=",adsorp_n,"\n")
    
        ##md_n 데이터 기록
        result[row][0] = timestep
        result[row][1] = adsorp_n
        row+=1
        timestep += t_inter
        
    return result


######################################################################################################

row_number = (frame_number-1)/50 + 1
final = np.zeros((int(row_number), num_tpr+1))
ffinal = np.zeros((int(row_number), num_tpr+2))
time = np.zeros((int(row_number), 1))

##tpr 분석 시작

for i in range(0, num_tpr+1):
    result = ratio(mda.Universe(f"md_{i}" + '.tpr', f"300ddvp_md_{i}" + '.xtc'))
    for j in range (0, int(row_number)):
        final[j][i] = result[j][1] 
    
avg_y = final.mean(axis=1)
avg=np.zeros((int(row_number),2))

for j in range (0, int(row_number)):
    time[j][0] = result[j][0]
    avg[j][0] = result[j][0]
    avg[j][1] = avg_y[j]

ffinal = np.concatenate((time, final),axis=1)

##최종 result 내보내기
df = pd.DataFrame(avg)
df.to_csv('avg_result.csv', index=False, header=False)

df = pd.DataFrame(ffinal)
df.to_csv('all_result.csv', index=False, header=False)

for i in range(0, num_tpr+1):
    df.to_csv(f"result_{i}" + ".csv", columns = [0, i+1], index=False, header=False)

    
    





























