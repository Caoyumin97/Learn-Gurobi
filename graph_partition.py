#!/usr/bin/env python
# coding: utf-8

# # Libs

# In[1]:


import warnings
warnings.filterwarnings('ignore')

import gurobipy as grb
from gurobipy import *
import numpy as np
import geatpy as ea

import torch
from torch.autograd import Variable


# # Graph Partition
# number of groups (consider a 5x5 net)
group_num = 16

# number of nodes
node_num = 100

# variables
ajacency = np.round(np.random.rand(node_num, node_num))
weight = np.random.rand(node_num, node_num) * ajacency

# ======== DEFINE Model ========
m = Model("Graph_Partition")

# ======== DEFINE DECISION VARIABLE ========
X = m.addVars(['X_' + str(i) for i in range(group_num * node_num)],
                lb=0, obj=1, vtype=GRB.BINARY)
X = np.array(X.values()).reshape(group_num, node_num)
m.update()


# ======== DEFINE OBJECTIVE ========
obj = 0
for group_idx in range(group_num):
    
    # calculate group & cut connection state
    group_vec = X[group_idx, :].reshape(-1, 1)
    cut_vec = 1 - X[group_idx, :].reshape(-1, 1)
    group_connection_state = ajacency * np.matmul(group_vec, group_vec.T)
    cut_connection_state = ajacency * np.matmul(cut_vec, cut_vec.T)
    
    # Obj Term 1: positive gain from high connectivity coefficient
    obj += np.sum(group_connection_state * weight)
    
    # Obj Term 2: elimination loss from cuts
    obj -= np.sum(cut_connection_state * weight)
m.setObjective(obj, GRB.MAXIMIZE)

# ======== DEFINE CONSTRAINTS========
for group_idx in range(group_num):
    
    # calculate group & cut connection state
    group_vec = X[group_idx, :].reshape(-1, 1)
    group_connection_state = ajacency * np.matmul(group_vec, group_vec.T)
    
    # Constr 1: nodes in same group must connect
    m.addConstr(np.sum(group_connection_state) == np.sum(group_vec))
    
    # Constr 2: each node should be contained in a group
    m.addConstr(np.sum(X, axis=1) == 1)
    
    # Constr 2: 0 < node number in each group < max_num
    m.addConstr(np.sum(group_vec) >= 1)
    m.addConstr(np.sum(group_vec) <= 10)


# In[ ]:


m.params.NonConvex = 2
m.optimize()
# solution_x = [v.x for v in m.getVars()]

