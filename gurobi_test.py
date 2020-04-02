#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import gurobipy as grb
from gurobipy import *
import numpy as np
#import geatpy


# In[ ]:


class SPP():
    def __init__(self, num_intervals, num_pairs, num_links, impact_intervals):
        self.num_intervals = num_intervals
        self.num_pairs = num_pairs
        self.num_links = num_links
        self.imp_intervals = impact_intervals
    
    
    def build_model(self, prior_od, prior_link, assignment_fraction, mcr):
        
        assert prior_od.shape == (self.num_intervals, self.num_pairs)
        assert prior_link.shape == (self.num_intervals, self.num_links)
        assert assignment_fraction.shape == (self.num_pairs, self.num_links, self.imp_intervals)
        
        # create model
        spp = Model("SPP")
        
        # create varibles
        X = spp.addVars(['X_' + str(i) for i in range(self.num_intervals * self.num_pairs)],
                        lb=0, obj=1, vtype=GRB.CONTINUOUS)
        Y = spp.addVars(['Y_' + str(i) for i in range(self.num_intervals * self.num_links)],
                        lb=0, obj=1, vtype=GRB.CONTINUOUS)
        X = np.array(X.values()).reshape(self.num_intervals, self.num_pairs)
        Y = np.array(Y.values()).reshape(self.num_intervals, self.num_links)
        spp.update()
        
        # build objective
        X_hat = prior_od
        Y_hat = prior_link
        X_var = np.var(prior_od)
        Y_var = np.var(prior_link)
        
        obj = 0
        for k in range(self.num_intervals):
            for i in range(self.num_pairs):
                obj += (X[k,i] - X_hat[k,i]) * (X[k,i] - X_hat[k,i]) / X_var
        
        for k in range(self.num_intervals):
            for j in range(self.num_links):
                obj += (Y[k,j] - Y_hat[k,j]) * (Y[k,j] - Y_hat[k,j]) / Y_var
        
        spp.setObjective(obj, GRB.MINIMIZE)
        
        # add constraints
        # No.1: Linear Assignment
        T = self.imp_intervals
        A = assignment_fraction
        for k in range(self.num_intervals):
            t = 0
            rhs = 0
            while (k - t >= 0) & (t <= T):
                rhs += (np.dot(X[k - t,:], A[:,:,t - 1])).reshape(1,self.num_links)
                t += 1
            for j in range(self.num_links):
                spp.addConstr(Y[k,j] == rhs[0,j])
        
        # No.2: Maximum Change Ratio
        beta = mcr
        for k in range(1,self.num_intervals):
            for i in range(self.num_pairs):
                spp.addConstr(X[k,i] >= (1 - beta) * X[k - 1,i])
                spp.addConstr(X[k,i] <= (1 + beta) * X[k - 1,i])
        
        return spp

    
    def estimate(self, prior_od, prior_link, assignment_fraction, mcr):
        
        spp = self.build_model(prior_od, prior_link, assignment_fraction, mcr)
        spp.optimize()
        
        x_range = self.num_intervals * self.num_pairs
        solution_x = [v.x for v in spp.getVars()[:x_range]]
        solution_y = [v.x for v in spp.getVars()[x_range:]]
        
        return solution_x, solution_y


# In[ ]:


class PRA():
    def __init__(self, num_intervals, num_pairs, num_links, impact_intervals):
        self.num_intervals = num_intervals
        self.num_pairs = num_pairs
        self.num_links = num_links
        self.imp_intervals = impact_intervals
    
    
    def build_model(self, prior_od, prior_link, prior_pr, probe_od, assignment_fraction, assignment_fraction_pr, mcr):
        
        assert prior_od.shape == (self.num_intervals, self.num_pairs)
        assert prior_link.shape == (self.num_intervals, self.num_links)
        assert prior_pr.shape == (self.num_intervals, self.num_links)
        assert probe_od.shape == (self.num_intervals, self.num_pairs)
        assert assignment_fraction.shape == (self.num_pairs, self.num_links, self.imp_intervals)
        assert assignment_fraction_pr.shape == (self.num_pairs, self.num_links, self.imp_intervals)
        
        # create model
        pra = Model("PRA")
        
        # create varibles
        X = pra.addVars(['X_' + str(i) for i in range(self.num_intervals * self.num_pairs)],
                        lb=0, obj=1, vtype=GRB.CONTINUOUS)
        Y = pra.addVars(['Y_' + str(i) for i in range(self.num_intervals * self.num_links)],
                        lb=0, obj=1, vtype=GRB.CONTINUOUS)
        Theta = pra.addVars(['T_' + str(i) for i in range(self.num_intervals * self.num_links)],
                        lb=0, ub=1, obj=1, vtype=GRB.CONTINUOUS)
        H = pra.addVars(['H_' + str(i) for i in range(self.num_intervals * self.num_pairs)],
                        lb=0, ub=1, obj=1, vtype=GRB.CONTINUOUS)
        X = np.array(X.values()).reshape(self.num_intervals, self.num_pairs)
        Y = np.array(Y.values()).reshape(self.num_intervals, self.num_links)
        Theta = np.array(Theta.values()).reshape(self.num_intervals, self.num_links)
        H = np.array(H.values()).reshape(self.num_intervals, self.num_pairs)
        pra.update()
        
        # build objective
        X_hat = prior_od
        Y_hat = prior_link
        Theta_hat = prior_pr # OD penetration rates
        X_var = np.var(X_hat)
        Y_var = np.var(Y_hat)
        Theta_var = np.var(Theta_hat)
        
        obj = 0
        for k in range(self.num_intervals):
            for i in range(self.num_pairs):
                obj += (X[k,i] - X_hat[k,i]) * (X[k,i] - X_hat[k,i]) / X_var
        
        for k in range(self.num_intervals):
            for j in range(self.num_links):
                obj += (Y[k,j] - Y_hat[k,j]) * (Y[k,j] - Y_hat[k,j]) / Y_var
                obj += (Theta[k,j] - Theta_hat[k,j]) * (Theta[k,j] - Theta_hat[k,j]) / Theta_var
        
        pra.setObjective(obj, GRB.MINIMIZE)
        
        # add constraints
        # No.1: Linear Assignment
        T = self.imp_intervals
        A = assignment_fraction
        for k in range(self.num_intervals):
            t = 0
            rhs = 0
            while (k - t >= 0) & (t <= T):
                rhs += (np.dot(X[k - t,:], A[:,:,t - 1])).reshape(1,self.num_links)
                t += 1
            for j in range(self.num_links):
                pra.addConstr(0.9 * Y[k,j] <= rhs[0,j])
                pra.addConstr(1.9 * Y[k,j] >= rhs[0,j])
        
        # No.2: Penetration Rate Assignment
        Z = probe_od
        Rho = assignment_fraction_pr
        for k in range(self.num_intervals):
            for i in range(self.num_pairs):
                pra.addConstr(X[k,i] * H[k,i] == Z[k,i])
        
        for k in range(self.num_intervals):
            t = 0
            rhs = 0
            while (k - t >= 0) & (t <= T):
                rhs += (np.dot(H[k - t,:], Rho[:,:,t - 1])).reshape(1,self.num_links)
                t += 1
            for j in range(self.num_links):
                pra.addConstr(0.9 * Theta[k,j] <= rhs[0,j])
                pra.addConstr(1.1 * Theta[k,j] >= rhs[0,j])
        
        # No.3: Maximum Change Ratio
        beta = mcr
        for k in range(1,self.num_intervals):
            for i in range(self.num_pairs):
                pra.addConstr(X[k,i] >= (1 - beta) * X[k - 1,i])
                pra.addConstr(X[k,i] <= (1 + beta) * X[k - 1,i])
        
        return pra

    
    def estimate(self, prior_od, prior_link, prior_pr, probe_od, assignment_fraction, assignment_fraction_pr, mcr):
        
        pra = self.build_model(prior_od, prior_link, prior_pr, probe_od, assignment_fraction, assignment_fraction_pr, mcr)
        pra.params.NonConvex = 2
        pra.optimize()
        
        x_range = self.num_intervals * self.num_pairs
        y_range = self.num_intervals * (self.num_pairs + self.num_links)
        t_range = self.num_intervals * (self.num_pairs + 2 * self.num_links)
        solution_x = [v.x for v in pra.getVars()[:x_range]]
        solution_y = [v.x for v in pra.getVars()[x_range:y_range]]
        solution_t = [v.x for v in pra.getVars()[y_range:t_range]]
        solution_h = [v.x for v in pra.getVars()[t_range:]]
        
        return solution_x, solution_y, solution_t, solution_h


# In[ ]:


# params init
num_intervals = 18
num_pairs = 120
num_links = 4
impact_intervals = 3

# prior information
prior_od = 10 + 5 * np.random.rand(num_intervals, num_pairs)
prior_link = 100 + 10 * np.random.rand(num_intervals, num_links)
prior_pr = 0.1 + 0.02 * np.random.normal(0,1,(num_intervals, num_links))
probe_od = 1 + 1 * np.random.rand(num_intervals, num_pairs)
assignment_fraction = np.random.rand(num_pairs,num_links,impact_intervals)
assignment_fraction_pr = np.random.rand(num_pairs,num_links,impact_intervals)


# In[ ]:


spp = SPP(num_intervals, num_pairs, num_links, impact_intervals)
x, y = spp.estimate(prior_od, prior_link,assignment_fraction, 0.2)

