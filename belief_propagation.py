  # -*- coding: utf-8 -*-
"""
Created on Sun Nov  4 09:36:19 2018

@author: Gene
"""
def Energy(j,i,spin,spin_n,J,h):
    
    Energy = -J[j,i]*spin*spin_n-h[j,i]*spin
    
    return Energy

###############################################################################

def Record_NN(Ny,Nx):
    
    NN_list = []
    
    #top edge
    NN_list.append([[0,1],[1,0]])                    #left-top corner
    i = 1
    while i < Nx-1:
        NN_list.append([[0,i-1],[0,i+1],[1,i]])
        i = i+1
    NN_list.append([[0,Nx-2],[1,Nx-1]])              #right-top corner
    
    #bulk
    j = 1
    while j < Ny-1:
        i = 0
        NN_list.append([[j,i+1],[j-1,i],[j+1,i]])    #left-most
        i = i+1
        while i < Nx-1:
            NN_list.append([[j,i-1],[j,i+1],[j-1,i],[j+1,i]])
            i = i+1
        NN_list.append([[j,i-1],[j-1,i],[j+1,i]])   #right-most 
        j = j + 1
    
    #bottom edge
    NN_list.append([[Ny-1,1],[Ny-2,0]])             #left-bottom corner
    i = 1
    while i < Nx-1:
        j = Ny-1
        NN_list.append([[j,i-1],[j,i+1],[j-1,i]])
        i = i + 1
    NN_list.append([[Ny-1,Nx-2],[Ny-2,Nx-1]])       #right-bottom corner
    
    return NN_list

##############################################################################
def all_transfer_matrices(Ny,Nx,J,h,T,NN):
    
    all_transfer_matrices = []
    
    for j in range(Ny):
        for i in range(Nx):
            q = i + j*Nx
            NN_len = len(NN[q])
            for k in range(NN_len):
                transfer_matrix = []
        
                for spin1 in range(1,-2,-2):
                    for spin2 in range(1,-2,-2):
                        val = m.exp(-Energy(j,i,spin1,spin2,J,h)/T)
                        transfer_matrix.append(round(val,6))
                            
                    
                all_transfer_matrices.append([[j,i],NN[q][k],transfer_matrix])                                  
    
    return all_transfer_matrices

##############################################################################

def messages_directions(j,i,Ny,Nx):
    
    Xdir = np.zeros((Ny,Nx-1))
    Ydir = np.zeros((Ny-1,Nx))
    
    #Directions of messages are translations of the first messages
    #we need for the target site
    
    #Direction along the x-axis in terms of Ny by (Nx-1) matrix: Xdir
    #Direction alonf the y-axis in terms of (Ny-1) by Nx matrix: Ydir
    
    #For example in a 3x3 lattice, and the target site is site (0,0)
    # (0,0) <- (0,1) <- (0,2)
    #   ^        ^        ^  
    # (1,0) <- (1,1) <- (1,2)
    #   ^        ^        ^
    # (2,0) <- (2,1) <- (2,2)
    #
    # Xdir = [[-1,-1],[-1,-1],[-1,1]]
    # Ydir = [[-1,-1,-1],[-1,-1,-1]]
    #
    # -1 : to the left or up
    # +1 : to the right or down
    
    for jp in range(Ny):
        for ip in range(Nx-1):
            if ip < i:
                Xdir[jp,ip] = 1   #to the right
            else:
                Xdir[jp,ip] = -1  #to the left
    for jp in range(Ny-1):
        for ip in range(Nx):
            if jp < j:
                Ydir[jp,ip] = 1  #to the top
            else:
                Ydir[jp,ip] = -1 #to the bottom
    
    
    return Xdir,Ydir

##############################################################################
def Tree_list(Ny,Nx,Xdir,Ydir,NN_list):
    
    tree_list = []
    target_list = []
    
    for j in range(Ny):
        for i in range(Nx):
            target_list.append([j,i])
            k = (i+j*Nx) 
                        
            NN_sites_len = len(NN_list[k])
            for k1 in range(NN_sites_len):
                wait_site = NN_list[k][k1]
                if wait_site not in target_list:
                    if wait_site[0]-j == 0 and wait_site[1]-i == 1: 
                        xdir = Xdir[j,i]
                        if xdir == -1:
                            tree_list.append([[j,i],wait_site])
                        else:
                            tree_list.append([wait_site,[j,i]])
                    elif wait_site[0]-j == 0 and wait_site[1]-i == -1:
                        xdir = Xdir[j,(i-1)]
                        if xdir == -1:
                            tree_list.append([[j,i],wait_site])
                        else:
                            tree_list.append([wait_site,[j,i]])
                    elif wait_site[0]-j == 1 and wait_site[1]-i == 0:
                        ydir = Ydir[j,i]
                        if ydir == -1:                            
                            tree_list.append([[j,i],wait_site])
                        else:
                            tree_list.append([wait_site,[j,i]])
                    elif wait_site[0]-j == -1 and wait_site[1]-i == 0:
                        ydir = Ydir[(j-1),i]
                        if ydir == -1:                            
                            tree_list.append([[j,i],wait_site])
                        else:
                            tree_list.append([wait_site,[j,i]])                   
    
    return tree_list

##############################################################################
def Tree_list_v2(T_list,Nx):
    
    # [target_site, wait_site];
    # but sites are referenced by number and not
    # by lattice coordinates
    #
    # Example: 3x3 lattice would have the labelling
    # 0 -- 1 -- 2
    # 3 -- 4 -- 5
    # 6 -- 7 -- 8 
    
    tree_list = []
    
    len_T_list = len(T_list)
    for k in range(len_T_list):
        pi = T_list[k][0][1]
        pj = T_list[k][0][0]
            
        qi = T_list[k][1][1]
        qj = T_list[k][1][0]
            
        s0 = (pi+pj*Nx)
        s1 = (qi+qj*Nx)
            
        tree_list.append([s0,s1])
    
    return tree_list

##############################################################################
def ordered_operations(j,i,tree_listv2,Nx):
    
    len_tree_list = len(tree_listv2)
    
    order = []
    order_rank = 1
    target_site = (i + j*Nx)    
    
    #Order of operations
    wait_lists = []
    for k in range(len_tree_list):
        if tree_listv2[k][0] == target_site:
            wait_lists.append(tree_listv2[k])
            
    order.append([order_rank,wait_lists])
    
    len_stop = 1
    len_order_old = 0
    
    while len_stop != 0:
    
        len_order = len(order)    
        order_rank = order_rank + 1
        
        
        for k in range(len_order_old,len_order):
            len_wait_lists = len(order[k][1])                        
                
            for l in range(len_wait_lists):
                target_site = order[k][1][l][1] 
                wait_lists = []
                for p in range(len_tree_list):
                    if tree_listv2[p][0] == target_site:
                        wait_lists.append(tree_listv2[p])
                                           
                order.append([order_rank,wait_lists])
         
        len_stop = len(wait_lists)    
        len_order_old = len_order
        
    #REMOVE EMPTY ELEMENTS
    
    len_order = len(order) 
    empty_elements = []
    for k in range(len_order):
        len_wl = len(order[k][1])
        if len_wl == 0:
            empty_elements.append(k)
            
    len_empty = len(empty_elements)        
    order = [order[k] for k in range(len_order-len_empty)]
    
    #REMOVE REDUNDANT ELEMENTS
    redundant_no = 1
    while redundant_no > 0:
        len_order = len(order)
        redundant_elements = []
        for k1 in range(len_order):
            for k2 in range(k1+1,len_order):
                if k2 != k1 and order[k2][0] == order[k1][0]:
                    if order[k2][1] == order[k1][1]:
                        redundant_elements.append([k1,k2])
                    
        redundant_no = len(redundant_elements)
        if redundant_no != 0:
            k_rem = redundant_elements[redundant_no-1][1]
            order.remove(order[k_rem])
        else:
            redundant_no = 0
        
    return order

################################################################################

def update_tree_list(k,tree_list,tr_matrix,target_site,wait_site):
    
    sub_tr_matrix = np.zeros((2,2))
    
    len_tree_list = len(tree_list)
                    
    count = 0
    for k_up in range(0,2):
        for k_down in range(0,2):
            val = tree_list[k][2][count]
            sub_tr_matrix[k_up,k_down] = val
            count = count+1 
 
    sub_tr_matrix = np.matmul(sub_tr_matrix,tr_matrix)
                                
    #finds the target_site_here 
    #and wait_site_here in the tree_listv2
    for s in range(len_tree_list):
        if target_site == tree_list[s][0] and wait_site == tree_list[s][1]:
            count = 0
            for k_up in range(0,2):
                for k_down in range(0,2):
                    #updates tree_listv2
                    tree_list[s][2][count] = sub_tr_matrix[k_up,k_down]                                        
                    count = count+1
    
    return tree_list

###########################################################################

def BP_algorithm(j,i,Ny,Nx,J,h,T,NN,ATM):
    
    #******STEP 1:
    
    #Directions of messages are "lattice translations" 
    #of the first messages we need for the target site
      
    Xdir, Ydir = messages_directions(j,i,Ny,Nx)
    #Xdir and Ydir make the graph directed in the expense of losing
    #some interactions with other nearest neighbors of a site
    
    #******STEP 2:
    
    tree_list = Tree_list(Ny,Nx,Xdir,Ydir,NN)
    #Make the tree: sites + messages direction
    #tree_list = |target_site|wait_site|
    #
    # target_site is the site queried
    # wait_site's are the NN sites directed to the target_sites
    
    
    tree_listv2 = Tree_list_v2(tree_list,Nx) 
    #easier to read tree_list
    #i.e., instead of lattice label [j,i] per site, 
    #we write [j,i] -> k = i + j*Nx ("particle" label)
    
    
    #*******STEP 3: THIS IS THE BP ALGORITHM
    
    
    order = ordered_operations(j,i,tree_listv2,Nx)
    #order of operations: list operations to conduct starting
    #from the largest rank (end nodes of the spanning tree) to
    #rank 1 which is the target_site
    #order = [rank,[[target_site,wait_site],...]
    
   
    #We retreive the transfer matrix from ATM 
    #for each element in the tree_list.
    #We rewrite tree_listv2 as [target_site,wait_site,transfer_matrix]
    
    len_ATM = len(ATM)
    len_tree_list = len(tree_list)
    for k in range(len_tree_list):
        target_site = tree_list[k][0]
        wait_site = tree_list[k][1]
        
        for s in range(len_ATM):
                
            ATM_ts = ATM[s][0]
            ATM_ws = ATM[s][1]
        
            if target_site == ATM_ts and wait_site == ATM_ws:
                transfer_matrix = ATM[s][2]
                
                tree_listv2[k] = [tree_listv2[k][0],tree_listv2[k][1],transfer_matrix]

    #We update tree_listv2 and by successive multiplications
    #of the transfer matrices dictated by the messages direction
    
    len_order = len(order)
    max_rank = order[len_order-1][0]
    
    rank = max_rank
    
    while rank >= 0:
        
        for k in range(len_order-1,0,-1):
                
            if order[k][0] == rank:
                
                len_wl = len(order[k][1])
                
                block_tr_matrix = np.zeros((2,2))
                    
                for p in range(len_wl):
                
                    target_site = order[k][1][p][0]
                    wait_site =  order[k][1][p][1]
                    
                    tr_matrix = np.zeros((2,2))
            
                    for s in range(len_tree_list):
                        if target_site == tree_listv2[s][0] and wait_site == tree_listv2[s][1]:
                       
                            count = 0
                            for k_up in range(0,2):
                                for k_down in range(0,2):
                                    tr_matrix[k_up,k_down] = tree_listv2[s][2][count]
                                    count = count+1 
                    
                    block_tr_matrix = block_tr_matrix + tr_matrix
                
                    
                    #In the next lower rank, we look for 
                    #wait_sites that are target_site of what we 
                    #previously calculated
                
                    for k_here in range(len_order-1,-1,-1):
                    
                        if order[k_here][0] == rank-1:   
                        
                            len_wl_here = len(order[k_here][1])
                        
                            for p_here in range(len_wl_here):
                                
                                wait_site_here = order[k_here][1][p_here][1]
                                target_site_here = order[k_here][1][p_here][0]
                                                        
                                if wait_site_here == target_site:
                                
                                    tree_listv2 = update_tree_list(k_here,tree_listv2,tr_matrix,
                                                                   target_site_here,wait_site_here)
                                    
                                                    
        rank = rank-1
        
    Belief_val = np.zeros(2)
        
    Belief_val[0] = m.exp(h[j,i])*(tree_listv2[0][2][0]+tree_listv2[0][2][1])
    Belief_val[1] = m.exp(-h[j,i])*(tree_listv2[0][2][2]+tree_listv2[0][2][3])
        
    Belief_val =  Belief_val/(Belief_val[0]+Belief_val[1]) 
        
    return Belief_val

    
###################################################################    

import time
start = time.time()

import numpy as np
import math as m
import matplotlib.pyplot as plt

#Parameters

Nx = 5    #no. of sites along the x-axis
Ny = 5    #no. of sites along the y-axis

#Temperature in units of Boltzmann constant kB
T = 10000.

#Interaction parameter
J0 = 1.
J = np.zeros((Ny,Nx))

#Option 1: Uniform interaction parameter
#--------------------------------------
J = J + J0        

#Option 2: Each site has randomly generated 
#          parameter but is the same for its nearest neighbors
#--------------------------------------
#for j in range(Ny):
#    for i in range(Nx):
#        J[j,i] = J0*np.random.uniform(-1,1)

#Field parameter
h0 = 1.
h = np.zeros((Ny,Nx))

#Option 1: Uniform field parameter
#--------------------------------------
h = h + h0

#Option 2: Each site has randomly generated field parameter 
#--------------------------------------
#h = h + h0*np.random.uniform(-1,1,size=(Ny,Nx))


NN = Record_NN(Ny,Nx) #records nearest neighbors of each site;
                      #in here so it does not have to be repeatedly calculated;
                      #needed for subparts of BP algorithm
                    
ATM = all_transfer_matrices(Ny,Nx,J,h,T,NN) #calculates transfer matrices
                                            #for all possible pairs

#BELIEF ALGORITHM

#Calculate belief per site per state

Belief_up = np.zeros((Ny,Nx))
Belief_down = np.zeros((Ny,Nx))


for j in range(Ny):
   for i in range(Nx):
        Belief_val = BP_algorithm(j,i,Ny,Nx,J,h,T,NN,ATM)
          
        Belief_up[j,i] = round(Belief_val[0],5)
        Belief_down[j,i] = round(Belief_val[1],5)
                      

#Magnetization       
Magnetization = Belief_up - Belief_down

end = time.time()
print('time elapsed =', round(end-start,5),' s')

plt.imshow(h,cmap='seismic')
plt.axis('off')
plt.title('Field')
plt.clim(vmin=-h0, vmax=h0)
plt.colorbar()
plt.show()

plt.imshow(Magnetization,cmap='seismic')
plt.axis('off')
plt.title('Magnetization')
plt.clim(vmin=-1, vmax=1)
plt.colorbar()
plt.show()