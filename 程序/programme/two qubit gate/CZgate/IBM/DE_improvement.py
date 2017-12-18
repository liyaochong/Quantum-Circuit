import numpy as np
import random
import time
from multiprocessing import Pool
import scipy.io as sio
import os

def evolution(evaluate_func,i,g,n,S,x_all,cr_list,f_list,x_u,x_l):
    
    if random.random() < S:#  Whole space 

        x_g_without_i = np.delete(x_all[g],i,0)
        np.random.shuffle(x_g_without_i)
        h_i = x_g_without_i[1]+f_list[g]*(x_g_without_i[2]-x_g_without_i[3])
        #变异操作后，h_i个体可能会过上下限区间，为了保证在区间以内对超过区间外的值赋值为相邻的边界值
        #先处理上边界，如果h_i[item]大于x_u则取x_u，如果小于则取h_i[item]
        h_i = [h_i[item] if h_i[item]<x_u[item] else x_u[item] for item in range(n)] 
        h_i = [h_i[item] if h_i[item]>x_l[item] else x_l[item] for item in range(n)] 

        #交叉操作，对变异后的个体，根据随机数与交叉阈值确定最后的个体
        # print(h_i)
        v_i = np.array([x_all[g][i][j] if (random.random() > cr_list[g] ) else h_i[j] for j in range(n) ])
        
    else:#subspace 
        index = int(n*random.random())  #choose the subspace
        v_i = x_all[g][i]
        x_g_without_i = np.delete(x_all[g],i,0)
        np.random.shuffle(x_g_without_i)
        v_i[index] = (x_g_without_i[1]+f_list[g]*(x_g_without_i[2]-x_g_without_i[3]))[index]
        v_i[index] = v_i[index] if v_i[index]<x_u[index] else x_u[index]
        v_i[index] = v_i[index] if v_i[index]>x_l[index] else x_l[index]
        v_i[index] = x_all[g][i][index] if (random.random() > cr_list[g]) else v_i[index]

    #根据评估函数确定是否更新新的个体
    value_vi = evaluate_func(v_i)
    return([v_i,value_vi])
    

def de(evaluate_func,n = 4, m_size = 20 , f = 0.5 , cr = 0.3 , S = 1 , iterate_time = 100 , x_l = np.array([0,1,0,2]),x_u = np.array([5,6,8,4]) ):
    num_CPU = 25
    #初始化
    x_all = np.zeros((iterate_time , m_size , n))#m_size为population ，n为dimension
    value = np.zeros((iterate_time , m_size ))
    f_list = np.zeros(iterate_time);f_list[0] = f
    cr_list = np.zeros(iterate_time);cr_list[0] = cr 
    
    if os.path.exists('result'):
        pass
    else:
        os.mkdir('result')

    initial = []
    p = Pool(num_CPU)
    for i in range(m_size):
        x_all[0][i] = x_l + random.random()*(x_u-x_l)
        initial.append(p.apply_async(evaluate_func,(x_all[0][i],)))
    value[0] = np.array([initial[i].get() for i in range(len(initial))])
    
    p.close()
    p.join()
    
    p = Pool(num_CPU)
    print('差分进化算法初始化完成')
    print('寻优参数维度为：',n)
    print('population为：',m_size)
    for g in range(iterate_time-1):
        print('第',g,'代')
        
        result = []
        for i in range(m_size):
            result.append(p.apply_async(evolution,(evaluate_func,i,g,n,S,x_all.copy(),cr_list.copy(),f_list.copy(),x_u.copy(),x_l.copy())))
        v_i = np.array([(result[i].get())[0] for i in range(len(result))])
        value_vi = np.array([(result[i].get())[1] for i in range(len(result))])
        for i in range(len(v_i)):
            x_all[g+1][i] = v_i[i] if value[g][i]>value_vi[i] else x_all[g][i]
            value[g+1][i] = value_vi[i] if value[g][i]>value_vi[i] else value[g][i]

        f_list[g+1] = 0.9+0.1*random.random() if random.random() < 0.1 else f_list[g]
        cr_list[g+1] = random.random() if random.random() < 0.1 else cr_list[g]
        print('Best parameters:',x_all[g+1][np.argmin(value[g+1])])    
        print('least cost function',np.min(value[g+1]))
        
        filename = './result/'+str(g)+'_'+time.strftime('%Y%m%d%X',time.localtime())+'.mat'
        sio.savemat(filename,{'x_all':x_all,'value':value,'f_list':f_list,'cr_list':cr_list,'best_para':x_all[g+1][np.argmin(value[g+1])],'min_fun':np.min(value[g+1])})
    evaluate_result = value[-1]
    best_x_i = x_all[-1][np.argmin(evaluate_result)]
    print('f_list:',f_list[-1])
    print('cr_list',cr_list[-1])
    print(evaluate_result)
    print('最小值：',np.min(evaluate_result))
    print('最佳参数：',best_x_i)
    p.close()
    p.join()

def evaluate_func(x):
    a = x[0]
    b = x[1]
    c = x[2]
    d = x[3]
    return(4*a**2 - 3*b + 5*c**3 - 6*d)
if __name__ == '__main__':
    de(evaluate_func) 