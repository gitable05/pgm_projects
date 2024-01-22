"""
Created on Mon Dec  3 08:31:21 2018

@author: Gene
"""

def data(train_fraction):
    
    data_file = open("mammographic_masses_data.txt","r")
    lines = data_file.read().split('\n')
    total = len(lines)
    total = total-1 #remove the last element that is empty
    
    data = []
    
    for i in lines:       
        j = i.split(',')
        jlen = len(j)
        val = []
        if jlen > 1:
            for k in range(jlen):
                if j[k] != '?':
                    elem = int(j[k])
                elif j[k] == '?':
                    elem = j[k]
            
                val.append(elem)
            
        data.append(val)
    
    data.remove(data[total])    
        
    train_set_no = int(train_fraction*total)
    
    train_data = []
    for i in range(train_set_no):
        train_data.append(data[i])
        
        
    test_data = []
    for i in range(train_set_no,total):
        test_data.append(data[i])

    data = []
    data_file.close()
    
    return train_data,test_data

###############################################################################
def simplify(data):
    
    len_TD = len(data)
    remove_B_eq_0 = []
    
    for i in range(len_TD):
        
        #BI-RADS
        if data[i][0] != '?':
            if data[i][0] == 0:
                remove_B_eq_0.append(i)
            #1,2,3 == benign assesment class
            elif data[i][0] <= 3:
                data[i][0] = 0
            #4,5,6 == malignant assesment class
            elif data[i][0] > 3:
                data[i][0] = 1
        
        #Age
        if data[i][1] != '?':
            # age <= 40
            if data[i][1] <= 40:                     
                data[i][1] = 0
            # age > 40
            else:                                      
                data[i][1] = 1
        
        #Shape
        if data[i][2] != '?':
            if data[i][2] <= 2:                     
                data[i][2] = 0
            else:                                      
                data[i][2] = 1
             
        #Mass
        if data[i][3] != '?':
            if data[i][3] <= 2:                     
                data[i][3] = 0
            else:                                      
                data[i][3] = 1
    
        #Density
        if data[i][4] != '?':
            if data[i][4] <= 2:                     
                data[i][4] = 0
            else:                                      
                data[i][4] = 1
                
    for i in range(len(remove_B_eq_0)):
        data.remove(data[remove_B_eq_0[i]])
        
    return data

###############################################################################

def divide(train_data):
    
    incomplete = []
    complete = []

    len_train = len(train_data)
    for i in range(len_train):
        if '?' in train_data[i]:
            incomplete.append(train_data[i])
        else:
            complete.append(train_data[i])
    
    return incomplete, complete

###############################################################################
def normalize(table,Node,Par):
    
    for i in range(Par):
        sumof = 0.
        for j in range(Node):
            sumof = sumof + table[j][i]
            
        for j in range(Node):
            table[j][i] = table[j][i]/sumof
            
    return table

################################################################################
def calculate_CPDs(complete_train_data):

    eps = 10e-10    #to avoid definitely zero probabilities
    len_comp_TD = len(complete_train_data)

    #PD for P(S)
    Prob_S = np.zeros(shape=(S_len,1)) + eps
    for i in range(len_comp_TD):
        S_val = complete_train_data[i][2]
        for k in range(S_len):
            if S_val == k:
                Prob_S[k] = Prob_S[k] + 1
            
    Prob_S = normalize(Prob_S,S_len,1)

    #PD for P(M)
    Prob_M = np.zeros(shape=(M_len,1)) + eps
    for i in range(len_comp_TD):
        M_val = complete_train_data[i][3]
        for k in range(M_len):
            if M_val == k:
                Prob_M[k] = Prob_M[k] + 1
            
    Prob_M = normalize(Prob_M,M_len,1)

    #PD for P(D)
    Prob_D = np.zeros(shape=(D_len,1)) + eps
    for i in range(len_comp_TD):
        D_val = complete_train_data[i][4]
        for k in range(D_len):
            if D_val == k:
                Prob_D[k] = Prob_D[k] + 1
            
    Prob_D = normalize(Prob_D,D_len,1)
    

    #Table-CPD for P(B|S,M,D)
    Cprob_B_SMD = np.zeros(shape=(B_len,S_len*M_len*D_len)) + eps
    
    for i in range(len_comp_TD):
        B_val = complete_train_data[i][0]
        for j in range(S_len):       
            for k in range(M_len):           
                for l in range(D_len):               
                    col_val = 4*j + 2*k + l #maps (j,k,l) to (0,1,...,7)
                    if (complete_train_data[i][2] == j and 
                        complete_train_data[i][3] == k and
                        complete_train_data[i][4] == l):
                        Cprob_B_SMD[B_val][col_val] = Cprob_B_SMD[B_val][col_val]+1

    Cprob_B_SMD = normalize(Cprob_B_SMD,B_len,S_len*M_len*D_len)

    #PD for P(A)
    Prob_A = np.zeros(shape=(A_len,1)) + eps
    for i in range(len_comp_TD):
        A_val = complete_train_data[i][1]
        for k in range(A_len):
            if A_val == k:
                Prob_A[k] = Prob_A[k] + 1
            
    Prob_A = normalize(Prob_A,A_len,1)


    #table-CPD for P(Se|A,B)
    Cprob_Se_BA = np.zeros(shape=(Se_len,B_len*A_len)) + eps
    
    for i in range(len_comp_TD):
        Se_val = complete_train_data[i][5]
        for j in range(B_len):
            for k in range(A_len):
                col_val = 2*j + k #maps (j,k) to (0,1,2,3)
                if (complete_train_data[i][0] == j and 
                    complete_train_data[i][1] == k):
                    Cprob_Se_BA[Se_val][col_val] = Cprob_Se_BA[Se_val][col_val]+1
            
    Cprob_Se_BA = normalize(Cprob_Se_BA,Se_len,B_len*A_len)
    
    return Prob_A, Prob_S, Prob_M, Prob_D, Cprob_B_SMD, Cprob_Se_BA

################################################################################
def combinations(max_miss_atts_val):
    
    if len(max_miss_atts_val) == 1:
        poss_val = []
        val = max_miss_atts_val[0]
        while val >= 0:
            poss_val.append([val])
            val = val-1          
        
    elif len(max_miss_atts_val) == 2:
        poss_val = []
        val0 = max_miss_atts_val[0]
        while val0 >= 0:
            val1 = max_miss_atts_val[1]
            while val1 >= 0:
                poss_val.append([val0,val1])
                val1=val1-1
            val0 = val0-1
                
    elif len(max_miss_atts_val) == 3:
        poss_val = []
        val0 = max_miss_atts_val[0]
        while val0 >= 0:
            val1 = max_miss_atts_val[1]
            while val1 >= 0:
                val2 = max_miss_atts_val[2]
                while val2 >=0:
                    poss_val.append([val0,val1,val2])
                    val2=val2-1
                val1 = val1-1
            val0 = val0-1
            
    elif len(max_miss_atts_val) == 4:
        poss_val = []
        val0 = max_miss_atts_val[0]
        while val0 >= 0:
            val1 = max_miss_atts_val[1]
            while val1 >= 0:
                val2 = max_miss_atts_val[2]
                while val2 >=0:
                    val3 = max_miss_atts_val[3]
                    while val3 >= 0:
                        poss_val.append([val0,val1,val2])
                        val3 = val3 - 1
                    val2=val2-1
                val1 = val1-1
            val0 = val0-1
                    
    return poss_val

###############################################################################
def Prob_vals_from(instances_set):
    
    P_vals = []
    for k in range(len(instances_set)):
            
        B_val = instances_set[k][0]
        A_val = instances_set[k][1]
        S_val = instances_set[k][2]
        M_val = instances_set[k][3]
        D_val = instances_set[k][4]
        Se_val = instances_set[k][5]
                        
        A_cont = Theta_A[A_val][0]
        S_cont = Theta_S[S_val][0]
        M_cont = Theta_M[M_val][0]
        D_cont = Theta_D[D_val][0]
            
        SMD_index = 4*S_val + 2*M_val + D_val
        B_SMD_cont = Theta_B_SMD[B_val][SMD_index]
        
        BA_index = 2*B_val + A_val
        Se_BA_cont = Theta_Se_BA[Se_val][BA_index]
        
        P = A_cont*S_cont*M_cont*D_cont*B_SMD_cont*Se_BA_cont
        
        P_vals.append([P])
            
    P_vals = normalize(P_vals,len(P_vals),1)

    return P_vals
################################################################################

import numpy as np
import matplotlib.pyplot as plt

# Mammographic mass data set from UC Irvine Machine Learning Repository
#               -- [CITE AS]

# Arrangement of data: |BI-RADS|Age|Shape|Margin|Density|Severity|

# [B]I-RADS assesment: 0 = incomplete, 1 = negative, 2 =  benign, 3 = probably benign
#                    4 = suspicious, 5 = highly suggestive of malignancy,
#                    6 = proven malignancy
# [A]ge: integer value
# [S]hape: 1 = Round, 2 = Oval, 3 = Lobular, 4 = Irregular
# [M]argin: 1 = Circumscribed, 2 = Microlobulated, 
#            3 = Obscured, 4 = Ill-defined, 5 = spiculated
# [D]ensity: 1 = High, 2 = Iso, 3 = low, 4 = Fat-containing
# [Se]verity: 0 = Benign, 1 = Malignant
#
# Shape, Margin and Density are BI RADS Attributes.
#
# The graph structure: {Shape,Margin,Density} -> BI RADS assesment, 
#                                {BI RADS, Age} -> Severity

B_len = A_len = S_len = M_len = D_len = Se_len = 2

#For simplify, {A,S,M,D} variables are changed into binaries
#BI-RADS = 0, instance is removed;
#BI-RADS <= 3 --> B = 0 ("Benign class"); BI-RADS  > 3 --> B = 1("Malignant class")
#Age < 40 --> A = 0; Age > 40 --> A = 1
#Shape = 1,2 ("Benign class") --> S = 0; Shape = 3,4 ("Malignant class") --> S = 1
#Margin = 1,2 ("Benign class") --> M = 0; Margin = 3,4 ("Malignant class") --> M = 1
#Density = 1,2 ("Benign class") --> D = 0; Density = 3,4 ("Malignant class") --> D = 1
#These changes are not fully justified medically but only for reducing the 
#possible values for these nodes, thereby reducing calculations.


fraction = 4/5 #of data as training set
train_data, test_data = data(fraction)

train_data = simplify(train_data)
test_data = simplify(test_data)
print('Total number of data=',len(train_data)+len(test_data))
print("Training data=",len(train_data),"Test data=",len(test_data))

#Separate instances that have complete data and that have missing data
incomplete_train_data, complete_train_data = divide(train_data)
train_data = [] #erase original list
print('      ')
print("Training data with complete attributes=",len(complete_train_data))
print("Training data with incomplete attribute/s=",len(incomplete_train_data))

# PART 1 -- COMPLETE DATA
# In calculate_CPDs we calculate the probabilities P(S), P(M), P(D), P(A)
# and the conditional probabilities P(B|S,M,D) and P(Se|B,A).
# They can be used as initial parameters for expectation maximization
(Prob_A, Prob_S, Prob_M, Prob_D, 
     Cprob_B_SMD, Cprob_Se_BA) = calculate_CPDs(complete_train_data)

print(' ')
print(Prob_A)
print(Prob_S)
print(Prob_M)
print(Prob_D)
print(Cprob_B_SMD)
print(Cprob_Se_BA)

#PART 2 -- INCOMPLETE DATA; EXPECTATION MAXIMIZATION ALGORITHM
len_incom_TD = len(incomplete_train_data)

atts_len = [B_len,A_len,S_len,M_len,D_len,Se_len] 

Theta_vals_ASMD = [] #for data output
Theta_vals_B_SMD = []
Theta_vals_Se_BA = []

#OPTION 1 ----- Initial theta parameters from complete data
Theta_S, Theta_M, Theta_D, Theta_A = (Prob_S.copy(), Prob_M.copy(),
                                     Prob_D.copy(), Prob_A.copy())
Theta_B_SMD, Theta_Se_BA, = Cprob_B_SMD.copy(), Cprob_Se_BA.copy()

#OPTION 2 ----- Arbitrary initial thete parameters
#Theta_S = np.zeros(shape=(S_len,1)) + 0.5
#Theta_M = np.zeros(shape=(M_len,1)) + 0.5
#Theta_D = np.zeros(shape=(D_len,1)) + 0.5
#Theta_A = np.zeros(shape=(A_len,1)) + 0.5
#Theta_B_SMD = normalize(np.zeros(shape=(B_len,S_len*M_len*D_len)) + 1.,
#                        B_len,S_len*M_len*D_len)
#Theta_Se_BA = normalize(np.zeros(shape=(Se_len,B_len*A_len)) + 1.,
#                        Se_len,B_len*A_len)

total_iterations = 100  
for i in range(total_iterations):
    
    Theta_vals_ASMD.append([Theta_A[0][0],Theta_S[0][0],Theta_M[0][0],Theta_D[0][0]])
    Theta_vals_B_SMD.append([Theta_B_SMD[0][0],Theta_B_SMD[0][1],Theta_B_SMD[0][2],
                        Theta_B_SMD[0][3],Theta_B_SMD[0][4],Theta_B_SMD[0][5],
                        Theta_B_SMD[0][6],Theta_B_SMD[0][7]])
    Theta_vals_Se_BA.append([Theta_Se_BA[0][0],Theta_Se_BA[0][1],
                        Theta_Se_BA[0][2],Theta_Se_BA[0][3]])
    
    
    #These store the M_theta[x,u]= sum_m P(x,u | o(m),theta)
    #We only deal with
    M_fS = np.zeros(shape=(S_len,1))
    M_fM = np.zeros(shape=(M_len,1))
    M_fD = np.zeros(shape=(D_len,1))
    M_fA = np.zeros(shape=(A_len,1))
    M_fB_SMD = np.zeros(shape=(B_len,S_len*M_len*D_len))
    M_fSe_BA = np.zeros(shape=(Se_len,B_len*A_len))
    
    for j in range(len_incom_TD):
    
        miss_atts = []
        obs_atts_val = []
        instance = incomplete_train_data[j]

        for k in range(6):
            if instance[k] == '?':
                miss_atts.append(k)
            else:
                obs_atts_val.append([k,instance[k]])
            
        max_miss_atts_val = []
        for k in range(len(miss_atts)):            
            max_miss_atts_val.append(atts_len[miss_atts[k]]-1)
            
        poss_val = combinations(max_miss_atts_val) 
        #all combinations of missing attributes
        
        instances_set = [] #will contain all possible instances given poss_val
        poss_instance = [0,0,0,0,0,0]
        
        len_pv = len(poss_val)
     
        for k in range(len_pv):
            count = 0
            for l in range(6):
                if instance[l] == '?':
                    poss_instance[l] = poss_val[k][count]
                    count = count + 1
                else:
                    poss_instance[l] = instance[l]
                    
            poss_instance_copy = poss_instance.copy()
            instances_set.append(poss_instance_copy)
        
        #E-step
        
        #Calculates joint distribution P(X,U|o(m),theta)
        P_vals = Prob_vals_from(instances_set)
        
        for k in range(len(instances_set)):
            
            B_val = instances_set[k][0]
            A_val = instances_set[k][1]
            S_val = instances_set[k][2]
            M_val = instances_set[k][3]
            D_val = instances_set[k][4]
            Se_val = instances_set[k][5]
            
            SMD_index = 4*S_val + 2*M_val + D_val
            BA_index = 2*B_val + A_val
            
            #for M_fA
            for l in range(A_len):
                if l == A_val:
                    M_fA[l][0] = M_fA[l][0] + P_vals[k]
                    
            #for M_fS
            for l in range(S_len):                
                if l == S_val:
                    M_fS[l][0] = M_fS[l][0] + P_vals[k]
                    
            #for M_fM
            for l in range(M_len):
                if l == M_val:
                    M_fM[l][0] = M_fM[l][0] + P_vals[k]
                    
            #for M_fD
            for l in range(D_len):
                if l == D_val:
                    M_fD[l][0] = M_fD[l][0] + P_vals[k]
                    
            #for M_fB_SMD
            for l in range(B_len):
                for m in range(S_len):
                    for n in range(M_len):
                        for o in range(D_len):
                            if (l == B_val and m == S_val 
                            and n == M_val and o == D_val):
                                SMD_index = 4*m + 2*n + o
                                M_fB_SMD[l][SMD_index] = (
                                    M_fB_SMD[l][SMD_index]
                                    + P_vals[k])
                                            
            #for M_fSe_BA
            for l in range(Se_len):
                for m  in range(B_len):
                    for n in range(A_len):
                        if (l == Se_val and m == B_val and n == A_val):
                            BA_index = 2*m + n
                            M_fSe_BA[l][BA_index] = (M_fSe_BA[l][BA_index] +
                                                    P_vals[k])
    
    #M-step
    Theta_A = normalize(M_fA,A_len,1)
    Theta_S = normalize(M_fS,S_len,1)
    Theta_M = normalize(M_fM,M_len,1)
    Theta_D = normalize(M_fD,D_len,1)
    Theta_B_SMD = normalize(M_fB_SMD,B_len,S_len*M_len*D_len)
    Theta_Se_BA = normalize(M_fSe_BA,Se_len,B_len*A_len)
    
#Last update values appending
Theta_vals_ASMD.append([Theta_A[0][0],Theta_S[0][0],Theta_M[0][0],Theta_D[0][0]])
Theta_vals_B_SMD.append([Theta_B_SMD[0][0],Theta_B_SMD[0][1],Theta_B_SMD[0][2],
                        Theta_B_SMD[0][3],Theta_B_SMD[0][4],Theta_B_SMD[0][5],
                        Theta_B_SMD[0][6],Theta_B_SMD[0][7]])
Theta_vals_Se_BA.append([Theta_Se_BA[0][0],Theta_Se_BA[0][1],
                        Theta_Se_BA[0][2],Theta_Se_BA[0][3]])


#print(Theta_values)
plt.plot(Theta_vals_ASMD)
plt.legend(['A=0','S=0','M=0','D=0'])
plt.xlim((0,10))
plt.ylim((0,1.1))
plt.xlabel('Iterations')
plt.ylabel('Theta')
plt.title('Age, Shape, Margin, Density' )
plt.show()

plt.plot(Theta_vals_B_SMD)
plt.legend(['B=0|0,0,0','B=0|0,0,1','B=0|0,1,0','B=0|0,1,1',
            'B=0|1,0,0','B=0|1,0,1','B=0|1,1,0','B=0|1,1,1'])
plt.xlabel('Iterations')
plt.ylabel('Theta')
plt.title('BI RADS | Shape, Mass, Density')
plt.show()

plt.plot(Theta_vals_Se_BA)
plt.legend(['Se=0|0,0','Se=0|0,1','Se=0|1,0','Se=0|1,1'])
plt.xlim((0,10))
plt.ylim((0,1.1))
plt.xlabel('Iterations')
plt.ylabel('Theta')
plt.title('Severity | BI RADS, Age')
plt.show()
        
print('A,S,M,D --->',Theta_vals_ASMD[total_iterations])
print('B|S,M,D --->',Theta_vals_B_SMD[total_iterations])
print('Se|B,A --->',Theta_vals_Se_BA[total_iterations])

#PART 3 -- CHECKING WITH TEST DATA
CP_B_SMD, CP_Se_BA = Theta_B_SMD.copy(), Theta_Se_BA.copy()
P_A, P_S, P_M, P_D = Theta_A.copy(), Theta_S.copy(), Theta_M.copy(), Theta_D.copy()

#With complete data
Agree_on_Se = 0.
Agree_on_B = 0.
sampling_size = 500

for h in range(sampling_size):

    for i in range(len(test_data)):
    
        instance = test_data[i]
        
        CP_B_SMD, CP_Se_BA = Theta_B_SMD.copy(), Theta_Se_BA.copy()
        P_A,P_S, P_M, P_D = Theta_A.copy(), Theta_S.copy(), Theta_M.copy(), Theta_D.copy()
    
        if instance[1] != '?':
            P_A[instance[1]][0] = 1.
            P_A[1-instance[1]][0] = 0.
        if instance[2] != '?':
            P_S[instance[2]][0] = 1.
            P_S[1-instance[2]][0] = 0.
        if instance[3] != '?':
            P_M[instance[3]][0] = 1.
            P_M[1-instance[3]][0] = 0.
        if instance[4] != '?':
            P_D[instance[4]][0] = 1.
            P_D[1-instance[4]][0] = 0.
    
        P_B_0_SMD = 0.
        for j in range(S_len):
            for k in range(M_len):
                for l in range(D_len):
                    SMD = 4*j+2*k+l
                    P_B_0_SMD = P_B_0_SMD + (CP_B_SMD[0][SMD]*
                                         P_S[j][0]*P_M[k][0]*P_D[l][0])
        
        r = np.random.uniform(0,1)
        if r < P_B_0_SMD:
            our_B = 0
        else:
            our_B = 1
                
        P_Se_0_BA = 0.
        for j in range(S_len):
            for k in range(M_len):
                for l in range(D_len):
                    for m in range(B_len):
                        for n in range(A_len):
                            SMD = 4*j+2*k+l
                            BA = 2*m + n
                            P_Se_0_BA = P_Se_0_BA + (CP_Se_BA[0][BA]*
                                                 CP_B_SMD[m][SMD]*
                                                 P_A[n][0]*P_S[j][0]*
                                                 P_M[k][0]*P_D[l][0])
                        
        r = np.random.uniform(0,1)
        if r < P_Se_0_BA:
            our_Se = 0
        else:
            our_Se = 1
        
        if our_Se == instance[5]:
            Agree_on_Se = Agree_on_Se + 1
        if our_B == instance[0]:
            Agree_on_B = Agree_on_B + 1
        
print(' ')
print('Comparison on predicting B and Se with',len(test_data),' test data')
Agree_on_B = round(Agree_on_B*100/(sampling_size*len(test_data)),2)
Agree_on_Se = round(Agree_on_Se*100/(sampling_size*len(test_data)),2)
print(Agree_on_B,'% agreement with test data on predicting BI RADS assesment')
print(Agree_on_Se,'% agreement with test data on predicting Severity')
    
                                                 
    
