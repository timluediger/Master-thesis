
"""
Created on Wed Oct 11 15:03:26 2023

@author: Tim Lüdiger
tim@luediger-net.de 

----MASTER THESIS----

Endogenous Technological Change in a General-Equilibrium Integrated
Assessment Model.

Two-region simulation calibrated for OECD (Region 1) and non-OECD (Region 2) 
countries.

Note, that this model assumes that the following assumptions are fullfilled:
    1. alpha_1 + nu_1 = alpha_2
    2. 0 < rho < 1
    3. rho < min{1/(1-alpha_1) , 1/(1-alpha_2)}.
    

Turn on excel_output and plots for the program to save the output and create 
the plots. To save the simulation output in an Excel file change the file name
at the end of the code correspondingly.
"""
import numpy as np
import scipy.optimize
import matplotlib.pyplot as plt
import datetime
import pandas as pd
from openpyxl import load_workbook
import os

# Set working directory to save output
os.chdir("C:/Users/tim/Documents/Studium/WWU/VWL_Master/Master thesis/Output")

# Set to "yes" if output and plots should be saved
excel_output = 'yes'
plots = "yes"

### A. Parameters
 
FP = 5000                    # Maximum iteration for the fix point problem
T = 25                       # One model period = 10 years
tol = 1e-10                  # Tolerance level in solution accuracy

# Climate Parameters
phi = 0.0228                 # Share of geometrically decaying atmospheric CO2
phi_L = 0.2                  # Share of permanent atmospheric CO2
phi_0 = 0.393                # Share of immediately vanishing atmospheric CO2
gamma_1 = 0.000053           # Damage parameter region 1
gamma_2 = 0.000053           # Damage parameter region 2
S_bar = 581                  # Pre industrial level of atmospheric carbon [GtC]

# Production parameters 
# Final output
alpha_0 = 0.3                # Capital share in final output production 
nu_0 = 0.04                  # Energy cost share in final output production 
kappa = 0.5                  # Shift parameter in energy-CES-aggregator  
rho = 0.5                    # Elasticity of substitution in energy-aggregator
Q_0_1 = 4.02                 # Productivity parameter in region 1
Q_0_2 = 0.66                 # Productivity parameter in region 2

#Dirty energy
nu_1 = 0.2                   # Share of fossil fuel in dirty energy production
alpha_1 = 0.3                # Capital share in dirty energy production
Q_1_1_init_0 = 4.2           # Initial dirty machine quality in region 1
Q_1_2_init_0 = 12            # Initial dirty machine quality in region 2

#Clean Energy
Q_2_1_init_0 = 18            # Initial clean energy productivity in region 1
Q_2_2_init_0 = 26            # Initial clean energy productivity in region 2  
alpha_2 = 0.5                # Capital share in clean energy sector

# Resource sector
c = 0.000071                 # Fossil fuel extraction costs 
R_0 = 50000                  # Initial fossil fuel resource stock
zeta = 0.5835                # Emission factor of fossil fuel consumption

# Consumer parameters 
beta = 0.985**10             # Discount factor 
sigma = 1                    # Elasticity of intertemporal substitution 
N_bar_1_0 = 0.18             # Initial population region 1
N_bar_2_0 = 0.82             # Initial population region 2
g = 0.16                     # Labor supply growth rate 

# R&D parameters
iota = 1.076                 # Innovation step size
lam = 0.06                   # Technology parameter lambda     
omega = 1.6                  # R&D cost factor


# Auxiliary variables for scientist allocation and profit ratio
varphi =   (rho*(1-alpha_2))/(rho-1)
varphi_1 = (rho*(1-alpha_1))/(rho-1)
Gamma =    ((1-alpha_2)**(-(varphi-1)) / (1-alpha_1)) * ((alpha_2 \
           **((1-alpha_2*varphi)/(1-alpha_2)))/ alpha_1**(-((alpha_1*varphi \
           *(1-alpha_1)-(1-alpha_2))/ (1-alpha_1)*(1-alpha_2))))
                                                             
Gamma_opt = (((1/alpha_2)-1)*alpha_2**(1/(1-alpha_2))*(1-alpha_2)**(-1) \
            *alpha_2**(-alpha_2/(1-alpha_2))*(1-alpha_2)**(-(varphi-1)) \
            * alpha_2**(-alpha_2*varphi/(1-alpha_2))) / (((1/alpha_1)-1)\
            * alpha_1**(1/(1-alpha_1))*alpha_1** (alpha_1*nu_1 / \
              ((1-alpha_1-nu_1)*(1-alpha_1)))*alpha_1**(-alpha_1/ \
              (1-alpha_1-nu_1))*alpha_1**(-alpha_1*varphi / (1-alpha_1-nu_1)))                                                             

# Policy parameters
#Laissez-faire (no tax):
tau_bar_lf = 0
# Optimal carbon tax:
tau_bar_eff = phi_L/(1-beta) + (1-phi_L)*phi_0/(1-beta*(1-phi))
tau_bar_sub = phi_L/(1-beta) + (1-phi_L)*phi_0/(1-beta*(1-phi))

Policies = [tau_bar_lf, tau_bar_eff, tau_bar_sub]
P = len(Policies)
epsilon = 0.001              # Arbitrarily small number to ensure subsidy

# Initial values
K_bar_1_0 = 0.18 * (0.62)    # Initial capital stock region 1
K_bar_2_0 = 0.18 * (1-0.62)  # Initial capital stock region 2
v_0 = c                      # Initial fossil fuel resource price
S_1_0 = 744                  # Initial level of permanent CO2
S_2_0 = 112                  # Initial level of non-permanent CO2



## B. Functions employed in the simulation

# CES-Aggregator of clean and dirty energy
def G(E_f, E_c):
    E_value = (kappa*E_f**rho + (1-kappa)*E_c**rho)**(1/rho)
    return E_value        
   
# Final output production function
def F_0(K, N, E, Q):
    return Q*K**alpha_0*N**(1-alpha_0-nu_0)*E**nu_0

# Dirty energy output production function
def F_1(N, X, M, Q):
    return N**(1-alpha_1-nu_1) * X**(nu_1) * M**(alpha_1) * Q**(1-alpha_1)

# Clean energy production function
def F_2(N, M, Q):   
    return N**(1-alpha_2)*M**(alpha_2) * Q**(1-alpha_2)

# Climate damage function 
def Damage(S, S_bar, gamma):
    return 1 - np.exp(-gamma * (S - S_bar) )         
                
# Global mean surface temperature
def Temperature(S, S_bar):
    return 3*np.log(S/S_bar)/np.log(2)  

# Update variables for fix point iteration
def update_variables(E_1_1, E_1_2, E_2_1, E_2_2, S, D_1, D_2, Q_1_1, Q_1_2,
                     Q_2_1, Q_2_2, M_1_1, M_1_2, M_2_1, M_2_2, K_1, K_2, N_0_1,
                     N_0_2, N_1_1, N_1_2, N_2_1, N_2_2, X_1, X_2):      
    E_1_1_new = F_1(N_1_1, X_1, M_1_1, Q_1_1) #Dirty energy production region 1
    E_1_2_new = F_1(N_1_2, X_2, M_1_2, Q_1_2) #Dirty energy production region 2
    E_2_1_new = F_2(N_2_1, M_2_1, Q_2_1)      #Clean energy production region 1
    E_2_2_new = F_2(N_2_2, M_2_2, Q_2_2)      #Clean energy production region 2
    E_1_new = G(E_1_1, E_2_1)                 #CES-aggregator region 1
    E_2_new = G(E_1_2, E_2_2)                 #CES-aggregator region 2
    Y_1_new = (1-D_1)*F_0(K_1, N_0_1, E_1, Q_0_1)
    Y_2_new = (1-D_2)*F_0(K_2, N_0_2, E_2, Q_0_2)
    return  (E_1_1_new, E_1_2_new, E_2_1_new, E_2_2_new, E_1_new, E_2_new,
             Y_1_new, Y_2_new)


# Interior solution to scientist allocation
def s_interior(s, eta_1, eta_2, Q_1_init, Q_2_init):
    if p==0:
        s_int = (1 - (eta_2/eta_1) * Gamma * (nu_1/(v_0 + tau*zeta))** \
                (nu_1*varphi/(1-alpha_2)) * ((eta_1*(1-s)*(alpha_1** \
                (alpha_1/(1-alpha_1))-1)+1)* (1+ iota*eta_1*(1-s)))** \
                (varphi_1+1)/((eta_2*s*(alpha_2**(alpha_2/(1-alpha_2))-1)+1)* \
                (1+ iota*eta_2*s))**(varphi+1) * Q_2_init**(-varphi)/ \
                Q_1_init**(-varphi_1))
    if p==1 or p==2:
        s_int = (1-  (eta_2/eta_1) * Gamma_opt * (nu_1/(v_0 + tau*zeta))** \
                (nu_1*varphi/(1-alpha_2)) * ((1+ iota*eta_1*(1-s)))** \
                (varphi_1+1) / ((1+ iota*eta_2*s))**(varphi+1) * \
                Q_2_init**(-varphi)/ Q_1_init**(-varphi_1))
    return s_int

# Optimal R&D subsidy
def subsidy(eta_1, eta_2, Q_1_init, Q_2_init):
    if Q_2_init**varphi /Q_1_init**varphi_1 < (eta_2/eta_1) * Gamma_opt * \
       (nu_1/(v_0 + tau*zeta))**(nu_1*varphi/(1-alpha_2)) * \
       ((1+ iota*eta_2)) \
       **(-(varphi+1)):
        theta = 0
    else:
        theta = ((eta_2/eta_1) * Gamma_opt *(nu_1/(v_0 + tau*zeta)) \
                **(nu_1*varphi/(1-alpha_2)))**(-1) *((1+ iota*eta_2)) \
                **(varphi+1) * Q_2_init**varphi/Q_1_init**varphi_1 -1 + epsilon
    return theta

# Scientist allocation                                                     
def scientist_allocation(eta_1, eta_2, Q_1_init, Q_2_init, theta):
    if p==0 :
        if Q_2_init**varphi < (1+theta) * (eta_2/eta_1) * Gamma * \
           (nu_1/(v_0 + tau*zeta))**(nu_1*varphi/(1-alpha_2)) * \
           ((eta_2*(alpha_2**(alpha_2/(1-alpha_2))-1)+1)* (1+ iota*eta_2))** \
           (-(varphi+1)) * Q_1_init**(varphi_1):
            s = 1
        elif Q_1_init**(-varphi_1) > (1+theta) * (eta_2/eta_1) * Gamma * \
             (nu_1/(v_0 + tau*zeta))**(nu_1*varphi/(1-alpha_2)) *  \
             ((eta_1*(alpha_1**(alpha_1/(1-alpha_2))-1)+1)* (1+ iota*eta_1)) \
             **(varphi_1+1) * Q_2_init**(-varphi):        
            s = 0
        else: s = scipy.optimize.fsolve(s_interior, x0 = 0.5, args = (eta_1, 
                                                                      eta_2, 
                                                                 Q_1_init,
                                                                 Q_2_init) ) 
    elif p==2 or p==1:
        if Q_2_init**varphi < (1+theta) * (eta_2/eta_1) * Gamma_opt * \
           (nu_1/(v_0 + tau*zeta))**(nu_1*varphi/(1-alpha_2)) * \
           ((1+ iota*eta_2))**(-(varphi+1)) * Q_1_init**(varphi_1):
            s = 1
        elif Q_1_init**(-varphi_1) > (1+theta) * (eta_2/eta_1) * Gamma_opt * \
             (nu_1/(v_0 + tau*zeta))**(nu_1*varphi/(1-alpha_2)) *  \
             ((1+ iota*eta_1))**(varphi_1+1) * Q_2_init**(-varphi):        
            s = 0
        else: s = scipy.optimize.fsolve(s_interior, x0 = 0.5, args = (eta_1, 
                                                                      eta_2, 
                                                                Q_1_init,
                                                                Q_2_init))
    return s          

    
## D. Arrays to store variables in  
K_0_1_series =np.empty((T+1,P))
K_0_2_series =np.empty((T+1,P))
N_0_1_series =np.empty((T+1,P))
N_0_2_series =np.empty((T+1,P))
N_1_1_series =np.empty((T+1,P))
N_1_2_series =np.empty((T+1,P))
N_2_1_series =np.empty((T+1,P))
N_2_2_series =np.empty((T+1,P))
N_bar_1_series =np.empty((T+1,P))
N_bar_2_series =np.empty((T+1,P))
X_1_series =np.empty((T+1,P))
X_2_series =np.empty((T+1,P))
E_1_1_series =np.empty((T+1,P))
E_1_2_series =np.empty((T+1,P))
E_2_1_series =np.empty((T+1,P))
E_2_2_series =np.empty((T+1,P))
E_1_series =np.empty((T+1,P))
E_2_series =np.empty((T+1,P))
S_1_series =np.empty((T+1,P))
S_2_series =np.empty((T+1,P))
S_series =np.empty((T+1,P))
Z_series =np.empty((T+1,P))
D_1_series =np.empty((T+1,P))
D_2_series =np.empty((T+1,P))
TEMP_series =np.empty((T+1,P))
tau_series =np.empty((T+1,P))
Y_1_series =np.empty((T+1,P))
Y_2_series =np.empty((T+1,P))
C_series =np.empty((T+1,P))
K_series =np.empty((T+1,P))
r_series =np.empty((T+1,P))
w_1_series =np.empty((T+1,P))
w_2_series =np.empty((T+1,P))
q_series =np.empty((T+1,P))
W_1_series =np.empty((T+1,P))
W_2_series =np.empty((T+1,P))   
p_1_1_series =np.empty((T+1,P))
p_2_1_series =np.empty((T+1,P))
p_1_2_series =np.empty((T+1,P))
p_2_2_series =np.empty((T+1,P))
eta_1_1_series = np.empty((T+1,P))
eta_2_1_series = np.empty((T+1,P))
eta_1_2_series = np.empty((T+1,P))
eta_2_2_series = np.empty((T+1,P))
chi_1_1_series = np.empty((T+1,P))
chi_2_1_series = np.empty((T+1,P))
chi_1_2_series = np.empty((T+1,P))
chi_2_2_series = np.empty((T+1,P))
R_series =np.empty((T+1,P))
Q_1_1_series =np.empty((T+1,P))
Q_2_1_series =np.empty((T+1,P))
Q_1_2_series =np.empty((T+1,P))
Q_2_2_series =np.empty((T+1,P))
H_1_1_series =np.empty((T+1,P))
H_2_1_series =np.empty((T+1,P))
H_1_2_series =np.empty((T+1,P))
H_2_2_series =np.empty((T+1,P))
M_1_1_series =np.empty((T+1,P))
M_2_1_series =np.empty((T+1,P))
M_1_2_series =np.empty((T+1,P))
M_2_2_series =np.empty((T+1,P))
s_1_series =np.empty((T+1,P))
s_2_series =np.empty((T+1,P))
theta_1_series =np.empty((T+1,P))
theta_2_series =np.empty((T+1,P))

## E. Iteration through policy scenarios and time t
begin_time = datetime.datetime.now() # performance check
for p in range(len(Policies)):
    tau_bar = Policies[p] 
    print('p', p)
    # Laissez-faire:
    if(p == 0): #laissez-faire, no tax
        C_0 = 0.5 # Initial consumption
    # Optimal tax:
    if(p == 1):
        C_0= 0.5   
    # Optimal tax and subsidy:
    if(p == 2):
        C_0= 0.5
    C_max = C_0*1.12 # Upper bound of initial consumption
    C_min = C_0*0.88 # Lower bound of initial consumption
     
    #Initialization of model    
    finished = 0   
    count_crit = 100
    counter = 0
 
    while(finished == 0 and counter < count_crit):
        counter += 1
    
        C_0  = (C_max + C_min)/2
        
        # Initial guesses for the iteration
        N_bar_1 = N_bar_1_0
        N_bar_2 = N_bar_2_0
        K = K_bar_1_0 + K_bar_2_0
        S_1_0   = 744
        S_2_0   = 112
        D_1_0   = 0
        D_2_0   = 0
        S_1 = S_1_0
        S_2 = S_2_0
        Q_1_1_init = Q_1_1_init_0
        Q_1_2_init = Q_1_2_init_0
        Q_2_1_init = Q_2_1_init_0
        Q_2_2_init = Q_2_2_init_0
        v = v_0
        R = R_0
        C_0 = (C_max + C_min)/2
        C = C_0
    
        # Initial values for fix point search
        E_1_1 = 0.09703
        E_1_2 = 0.26575
        E_2_1 = 0.00626
        E_2_2 = 0.00674
        E_1 = 0.0381562
        E_2 = 0.0893
        Y_1 = 0.6023
        Y_2 = 0.1999

        for t in range(T):
            
            for iteration in range(FP):
                E_1_1_old, E_2_1_old, E_1_2_old, E_2_2_old, Y_1_old, Y_2_old, \
                E_1_old, E_2_old  = E_1_1, E_2_1, E_1_2, E_2_2, Y_1, Y_2, \
                E_1, E_2
                
                Q_1_1 = Q_1_1_init
                Q_1_2 = Q_1_2_init
                Q_2_1 = Q_2_1_init
                Q_2_2 = Q_2_2_init
                
                # Optimal carbon tax
                tau = tau_bar*(gamma_1*Y_1_old + gamma_2*Y_2_old)
                
                # Energy mix
                chi_1_1 = kappa * (E_1_1_old / E_1_old)**rho
                chi_1_2 = kappa * (E_1_2_old / E_2_old)**rho
                chi_2_1 = (1-kappa) * (E_2_1_old / E_1_old)**rho
                chi_2_2 = (1-kappa) * (E_2_2_old / E_2_old)**rho
                
                # Factor allocation for labor, capital and fossil resources
                n_0 = 1-alpha_0-nu_0
                n_1_1 = nu_0*(1-alpha_1-nu_1)*chi_1_1
                n_1_2 = nu_0*(1-alpha_1-nu_1)*chi_1_2
                n_2_1 = nu_0*(1-alpha_2)*chi_2_1
                n_2_2 = nu_0*(1-alpha_2)*chi_2_2
                N_0_1 = (n_0/(n_0 + n_1_1 + n_2_1))*N_bar_1
                N_1_1 = (n_1_1/(n_0 + n_1_1 + n_2_1))*N_bar_1
                N_2_1 = (n_2_1/(n_0 + n_1_1 + n_2_1))*N_bar_1
                N_0_2 = (n_0/(n_0 + n_1_2 + n_2_2))*N_bar_2
                N_1_2 = (n_1_2/(n_0 + n_1_2 + n_2_2))*N_bar_2
                N_2_2 = (n_2_2/(n_0 + n_1_2 + n_2_2))*N_bar_2
                
                K_1 = Y_1_old / (Y_1_old + Y_2_old) * K
                K_2 = Y_2_old / (Y_1_old + Y_2_old) * K
                
                X_1 = (nu_0 * nu_1 * Y_1_old * chi_1_1)/(c + zeta*tau)
                X_2 = (nu_0 * nu_1 * Y_2_old * chi_1_2)/(c + zeta*tau)
                
                # Climate variables and Emissions
                Z = zeta*(X_1 + X_2)            
                S = S_1 + phi_L*Z + (1-phi)*S_2 + (1-phi_L)*phi_0*Z
                D_1 = Damage(S, S_bar, gamma_1) 
                D_2 = Damage(S, S_bar, gamma_2) 
                
                # Interest rate
                r = (alpha_0*Y_1_old)/K_1
                
                # Regional R&D success probabilities
                eta_1_1 = (omega**(-1) * (lam/(1+lam))**lam *
                          (((1-alpha_1)*alpha_1**((1+alpha_1)/(1-alpha_1))* 
                          (nu_0*Y_1_old/E_1_1_old*kappa*
                          (E_1_1_old/E_1_old)**rho)**(1/(1-alpha_1))*
                          X_1**(nu_1/(1-alpha_1))*N_1_1**((1-alpha_1-nu_1)/
                          (1-alpha_1)))/r)**lam)
                
                eta_1_2 = (omega**(-1) * (lam/(1+lam))**lam *
                          (((1-alpha_1)*alpha_1**((1+alpha_1)/(1-alpha_1))* 
                          (nu_0*Y_2_old/E_1_2_old*kappa*
                          (E_1_2_old/E_2_old)**rho)**(1/(1-alpha_1))*
                          X_2**(nu_1/(1-alpha_1))*N_1_2**((1-alpha_1-nu_1)/
                          (1-alpha_1)))/r)**lam)
                
                eta_2_1 = (omega**(-1) * (lam/(1+lam))**lam * 
                          (((1-alpha_2)*alpha_2**((1+alpha_2)/(1-alpha_2))*
                          (nu_0*Y_1_old/E_2_1_old*(1-kappa)*(E_2_1_old/E_1_old)
                          **rho)**(1/(1-alpha_2))*N_2_1)/r)**lam)
                
                eta_2_2 = (omega**(-1) * (lam/(1+lam))**lam * 
                          (((1-alpha_2)*alpha_2**((1+alpha_2)/(1-alpha_2))*
                          (nu_0*Y_2_old/E_2_2_old*(1-kappa)*(E_2_2_old/E_2_old)
                          **rho)**(1/(1-alpha_2))*N_2_2)/r)**lam)
                
                # Subsidy to redirect innovation
                if p==2:
                    theta_1 = subsidy(eta_1_1, eta_2_1, Q_1_1_init, Q_2_1_init)
                    theta_2 = subsidy(eta_1_2, eta_2_2, Q_1_2_init, Q_2_2_init)
                else:
                    theta_1 =0
                    theta_2 =0
                
                # Scientist allocation

                s_1 = scientist_allocation(eta_1_1, eta_2_1,
                                       Q_1_1_init, Q_2_1_init, theta_1)
                s_2 = scientist_allocation(eta_1_2, eta_2_2,
                                       Q_1_2_init, Q_2_2_init, theta_2)
                
                # Update qualitivity level based on the law of motion
                Q_1_1 = Q_1_1 * (1 + iota*eta_1_1*(1-s_1))
                Q_2_1 = Q_2_1 * (1 + iota*eta_2_1*s_1)
                Q_1_2 = Q_1_2 * (1 + iota*eta_1_2*(1-s_2))
                Q_2_2 = Q_2_2 * (1 + iota*eta_2_2*s_2)
                
                # Determining sector specific machine demand 
                
                if p==0 :
                    # Region 1
                    M_1_1_m = (alpha_1**2 * nu_0*Y_1_old/E_1_1_old*kappa \
                              *(E_1_1_old/E_1_old)**rho)**(1/(1-alpha_1)) \
                              *X_1**(nu_1/(1-alpha_1))*N_1_1** \
                              ((1-alpha_1-nu_1) /(1-alpha_1))*Q_1_1 
                    M_1_1_c = (alpha_1 * nu_0*Y_1_old/E_1_1_old*kappa \
                              *(E_1_1_old/E_1_old)**rho)**(1/(1-alpha_1)) \
                              *X_1**(nu_1/(1-alpha_1))*N_1_1** \
                              ((1-alpha_1-nu_1) /(1-alpha_1))*Q_1_1 
                    M_2_1_m = (alpha_2**2 * nu_0*Y_1_old/E_2_1_old*(1-kappa) \
                              *(E_2_1_old/E_1_old)**rho)**(1/(1-alpha_2)) \
                              * N_2_1 * Q_2_1
                    M_2_1_c = (alpha_2 * nu_0*Y_1_old/E_2_1_old*(1-kappa) \
                              *(E_2_1_old/E_1_old)**rho)**(1/(1-alpha_2)) \
                              * N_2_1 * Q_2_1
                    M_1_1 = (eta_1_1*(1-s_1)*M_1_1_m + (1-eta_1_1*(1-s_1))\
                            *M_1_1_c)
                    M_2_1 = (eta_2_1*s_1*M_2_1_m + (1-eta_2_1*s_1)*M_2_1_c)
                    
                    # Region 2
                    M_1_2_m = (alpha_1**2 * nu_0*Y_2_old/E_1_2_old*kappa \
                              *(E_1_2_old/E_2_old)**rho)**(1/(1-alpha_1)) \
                              *X_2**(nu_1/(1-alpha_1))*N_1_2** \
                              ((1-alpha_1-nu_1) /(1-alpha_1))*Q_1_2 
                    M_1_2_c = (alpha_1 * nu_0*Y_2_old/E_1_2_old*kappa \
                              *(E_1_2_old/E_2_old)**rho)**(1/(1-alpha_1)) \
                              *X_2**(nu_1/(1-alpha_1))*N_1_2** \
                              ((1-alpha_1-nu_1) /(1-alpha_1))*Q_1_2 
                    M_2_2_m = (alpha_2**2 * nu_0*Y_2_old/E_2_2_old*(1-kappa) \
                              *(E_2_2_old/E_2_old)**rho)**(1/(1-alpha_2)) \
                              * N_2_2 * Q_2_2
                    M_2_2_c = (alpha_2 * nu_0*Y_2_old/E_2_2_old*(1-kappa) \
                              *(E_2_2_old/E_2_old)**rho)**(1/(1-alpha_2)) \
                              * N_2_2 * Q_2_2
                    M_1_2 = (eta_1_2*(1-s_2)*M_1_2_m + (1-eta_1_2*(1-s_2))\
                            *M_1_2_c)
                    M_2_2 = (eta_2_2*s_2*M_2_2_m + (1-eta_2_2*s_2)*M_2_2_c)
                elif p==2 or p==1:
                    M_1_1 = (alpha_1 * nu_0*Y_1_old/E_1_1_old*kappa \
                              *(E_1_1_old/E_1_old)**rho)**(1/(1-alpha_1)) \
                              *X_1**(nu_1/(1-alpha_1))*N_1_1** \
                              ((1-alpha_1-nu_1) /(1-alpha_1))*Q_1_1 
                    M_1_2 = (alpha_1 * nu_0*Y_2_old/E_1_2_old*kappa \
                              *(E_1_2_old/E_2_old)**rho)**(1/(1-alpha_1)) \
                              *X_2**(nu_1/(1-alpha_1))*N_1_2**\
                              ((1-alpha_1-nu_1) /(1-alpha_1))*Q_1_2 
                    M_2_1 = (alpha_2 * nu_0*Y_1_old/E_2_1_old*(1-kappa) \
                              *(E_2_1_old/E_1_old)**rho)**(1/(1-alpha_2)) \
                              * N_2_1 * Q_2_1
                    M_2_2 = (alpha_2 * nu_0*Y_2_old/E_2_2_old*(1-kappa) \
                              *(E_2_2_old/E_2_old)**rho)**(1/(1-alpha_2)) \
                              * N_2_2 * Q_2_2           
                        
                E_1_1, E_1_2, E_2_1, E_2_2, E_1, E_2, Y_1, Y_2 \
                = update_variables(E_1_1, E_1_2, E_2_1, E_2_2, S, D_1, D_2, 
                                   Q_1_1, Q_1_2, Q_2_1, Q_2_2, M_1_1, M_1_2,
                                   M_2_1, M_2_2, K_1, K_2, N_0_1, N_0_2, N_1_1,
                                   N_1_2, N_2_1, N_2_2, X_1, X_2)
                
                if(
                    abs(E_1_1 - E_1_1_old) < tol and  
                    abs(E_1_2 - E_1_2_old) < tol and  
                    abs(E_2_1 - E_2_1_old) < tol and 
                    abs(E_2_2 - E_2_2_old) < tol and
                    abs(E_1 - E_1_old) < tol and
                    abs(E_2 - E_2_old) < tol and
                    abs(Y_1 - Y_1_old) < tol and
                    abs(Y_2 - Y_2_old) < tol
                    ):
                    #print("Breaking the loop at i =", iteration)
                    break       
                
            # Resource stock
            R = R - X_1 - X_2
                        
            # Climate variables
            S_1 = S_1 + phi_L*Z
            S_2 = (1-phi)*S_2 + (1-phi_L)*phi_0*Z
            S = S_1 + S_2
            TEMP =  Temperature(S, S_bar) # 3*np.log(S/S_bar)/np.log(2)
    
            # Equilibrium factor prices and consumption 
            w_1 = (1-alpha_0 - nu_0)*Y_1/N_0_1
            w_2 = (1-alpha_0 - nu_0)*Y_2/N_0_2
            p_1_1 = nu_0*Y_1/E_1_1*kappa*(E_1_1/E_1)**rho
            p_2_1 = nu_0*Y_1/E_2_1*(1-kappa)*(E_2_1/E_1)**rho
            p_1_2 = nu_0*Y_2/E_1_2*kappa*(E_1_2/E_2)**rho
            p_2_2 = nu_0*Y_2/E_2_2*(1-kappa)*(E_2_2/E_2)**rho      
            C = C*(beta*r)**(1/sigma)
            
            # Energy mix
            chi_1_1 = p_1_1*E_1_1/(p_1_1*E_1_1 + p_2_1*E_2_1)
            chi_2_1 = p_2_1*E_2_1/(p_1_1*E_1_1 + p_2_1*E_2_1)
            chi_1_2 = p_1_2*E_1_2/(p_1_2*E_1_2 + p_2_2*E_2_2)
            chi_2_2 = p_2_2*E_2_2/(p_1_2*E_1_2 + p_2_2*E_2_2)
        
            # R&D expenditures in equilibrium
            H_1_1 = (lam/(1+lam)) * (((1-alpha_1)*alpha_1**((1+alpha_1) \
                    /(1-alpha_1))* (nu_0*Y_1/E_1_1*kappa*(E_1_1/E_1)**rho) \
                    **(1/(1-alpha_1))*X_1**(nu_1/(1-alpha_1))*N_1_1 \
                    **((1-alpha_1-nu_1)/(1-alpha_1)))/r) * Q_1_1 * (1-s_1)
            H_1_2 = (lam/(1+lam)) * (((1-alpha_1)*alpha_1**((1+alpha_1) \
                    /(1-alpha_1))* (nu_0*Y_2/E_1_2*kappa*(E_1_2/E_2)**rho) \
                    **(1/(1-alpha_1))*X_2**(nu_1/(1-alpha_1))*N_1_2 \
                    **((1-alpha_1-nu_1)/(1-alpha_1)))/r) * Q_1_2 * (1-s_2)
            H_2_1 = (lam/(1+lam))*(((1-alpha_2)*alpha_2**((1+alpha_2) \
                    /(1-alpha_2))* (nu_0*Y_1/E_2_1*(1-kappa)*(E_2_1/E_1) \
                    **rho)**(1/(1-alpha_2))*N_2_1)/r) * Q_2_1 * s_1
            H_2_2 = (lam/(1+lam))*(((1-alpha_2)*alpha_2**((1+alpha_2) \
                    /(1-alpha_2))* (nu_0*Y_2/E_2_2*(1-kappa)*(E_2_2/E_2) \
                    **rho)**(1/(1-alpha_2))*N_2_2)/r) * Q_2_2 * s_2
                             
            # Arrow-Debreu-prices
            if(t==0):
                q = 1
            else:
                q = q/r
            
            # Storing equilibrium values  
            K_series[t,p] = K  
            r_series[t,p] = r
            q_series[t,p] = q
            C_series[t,p] = C
            Z_series[t,p] = Z
            S_series[t,p] = S
            S_1_series[t,p] = S_1
            S_2_series[t,p] = S_2
            R_series[t,p] = R
            X_1_series[t,p] = X_1
            X_2_series[t,p] = X_2
            Y_1_series[t,p] = Y_1
            Y_2_series[t,p] = Y_2
            E_1_series[t,p] = E_1
            E_2_series[t,p] = E_2
            E_1_1_series[t,p] = E_1_1
            E_1_2_series[t,p] = E_1_2
            E_2_1_series[t,p] = E_2_1
            E_2_2_series[t,p] = E_2_2
            p_1_1_series[t,p] = p_1_1
            p_1_2_series[t,p] = p_1_2
            p_2_1_series[t,p] = p_2_1
            p_2_2_series[t,p] = p_2_2
            N_0_1_series[t,p] = N_0_1
            N_1_1_series[t,p] = N_1_1
            N_2_1_series[t,p] = N_2_1
            N_0_2_series[t,p] = N_0_2
            N_1_2_series[t,p] = N_1_2
            N_2_2_series[t,p] = N_2_2
            M_1_1_series[t,p] = M_1_1
            M_2_1_series[t,p] = M_2_1
            M_1_2_series[t,p] = M_1_2
            M_2_2_series[t,p] = M_2_2
            w_1_series[t,p] = w_1
            w_2_series[t,p] = w_2
            D_1_series[t,p] = D_1
            D_2_series[t,p] = D_2
            TEMP_series[t,p] = TEMP
            tau_series[t,p] = tau
            chi_1_1_series[t,p] = chi_1_1   
            chi_2_1_series[t,p] = chi_2_1   
            chi_1_2_series[t,p] = chi_1_2   
            chi_2_2_series[t,p] = chi_2_2   
            eta_1_1_series[t,p] = eta_1_1     
            eta_2_1_series[t,p] = eta_2_1 
            eta_1_2_series[t,p] = eta_1_2 
            eta_2_2_series[t,p] = eta_2_2 
            N_bar_1_series[t,p] = N_bar_1
            N_bar_2_series[t,p] = N_bar_2
            s_1_series[t,p] = s_1
            s_2_series[t,p] = s_2
            Q_1_1_series[t,p] = Q_1_1
            Q_2_1_series[t,p] = Q_2_1
            Q_1_2_series[t,p] = Q_1_2
            Q_2_2_series[t,p] = Q_2_2
            H_1_1_series[t,p] = H_1_1
            H_2_1_series[t,p] = H_2_1
            H_1_2_series[t,p] = H_1_2
            H_2_2_series[t,p] = H_2_2
            theta_1_series[t,p] = theta_1
            theta_2_series[t,p] = theta_2
            
            
            # Determining state variables for next period 
            K_old = K
            K = (Y_1 + Y_2 - C - c*(X_1 + X_2) - M_1_1 - M_1_2 - M_2_1 - M_2_2
                 - H_1_1 - H_2_1 - H_1_2 - H_2_2)
            N_bar_1 = N_bar_1*(1+g)
            N_bar_2 = N_bar_2*(1+g)
            Q_1_1_init = Q_1_1
            Q_1_2_init = Q_1_2
            Q_2_1_init = Q_2_1
            Q_2_2_init = Q_2_2
                    
            
            if(K<0.1*((Y_1 + Y_2 - c*(X_1 + X_2) - M_1_1 - M_1_2 - M_2_1 
                       - M_2_2 - H_1_1 - H_2_1 - H_1_2 - H_2_2))):
                C_max = C_0
                break
            if(C<0.1*(Y_1 + Y_2 - c*(X_1 + X_2) - M_1_1 - M_1_2 - M_2_1 
                       - M_2_2 - H_1_1 - H_2_1 - H_1_2 - H_2_2)):
                C_min = C_0
                break
            
            if(t==T):    
                finished = 1       


print("End of simulation! :-)")
print('Execution time is', datetime.datetime.now() - begin_time) 
        

## F. Plots


if plots == "yes":
    # Settings
    # Convert x-axis to years
    # Assuming each x value represents a decade starting from 2020
    x = np.arange(len(D_1_series[:20,0]*100))
    years = [2020 + i * 10 for i in x]  
    
    
    
    # 1.a. Share of green input in % (Region 1)
    plt.figure(1)
    plt.figure(dpi=500)
    plt.rc('figure', figsize=[6.22,3.84])
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
    plt.plot(years,chi_2_1_series[:20,0]*100 , label="Laissez-faire", linewidth=2.5)
    plt.plot(years,chi_2_1_series[:20,1]*100, label = "Efficient tax" , linewidth=2.5)
    plt.plot(years,chi_2_1_series[:20,2]*100, label = "Optimal", linewidth=2.5 )
    plt.xlabel(r't', loc="right", labelpad=1, fontsize = 18)
    plt.xticks(np.arange(2030, 2200, step=40), fontsize= 18)
    plt.gca().xaxis.set_label_coords(1.005, -0.02)  # Set the position of the x-axis label
    plt.ylabel("$\chi_2^1$ [in \%]", loc="top", rotation = 0, fontsize = 18)
    plt.yticks(np.arange(0,100,step = 20),fontsize=18)
    plt.gca().yaxis.set_label_coords(0.2, 0.9)
    plt.grid(axis ="y", which='both', color='gray', linestyle='--', linewidth=0.5)
    plt.xlim(2020, max(years))  
    plt.ylim(0,101)
    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['right'].set_visible(False)
    plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.15), shadow=False, ncol=2, fontsize = 16, frameon = False) 
    plt.savefig('Green_input_Re1.png', dpi=300, bbox_inches='tight')

    
    # 1.b. Share of green input in % (Region 2)
    plt.figure(2)
    plt.figure(dpi=500)
    plt.rc('figure', figsize=[6.22,3.84])
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
    plt.plot(years,chi_2_2_series[:20,0]*100 , label="Laissez-faire", linewidth=2.5)
    plt.plot(years,chi_2_2_series[:20,1]*100, label = "Efficient tax" ,linewidth=2.5)
    plt.plot(years,chi_2_2_series[:20,2]*100, label = "Optimal",linewidth=2.5 )
    plt.xlabel(r't', loc="right", labelpad=1, fontsize = 18)
    plt.xticks(np.arange(2030, 2200, step=40), fontsize= 18)
    plt.gca().xaxis.set_label_coords(1.005, -0.02)  # Set the position of the x-axis label
    plt.ylabel("$\chi_2^2$ [in \%]", loc="top", rotation = 0, fontsize = 18)
    plt.yticks(np.arange(0, 100, step=20),fontsize=18)
    plt.gca().yaxis.set_label_coords(0.2, 0.9)
    plt.grid(axis ="y", which='both', color='gray', linestyle='--', linewidth=0.5)
    plt.xlim(2020, max(years))  
    plt.ylim(0,101)
    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['right'].set_visible(False)
    plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.15), shadow=False, ncol=2, fontsize = 16, frameon = False) 
    plt.savefig('Green_input_Re2.png', dpi=300, bbox_inches='tight')
    
    # 2. Percentage change in net output
    
    delta_growth_1 = np.empty((21,P))
    delta_growth_2 = np.empty((21,P))
    for i in range(1,21):
        delta_growth_1[i] = ((Y_1_series[i]) - \
                            (Y_1_series[i-1])) / \
                            (Y_1_series[i-1])
        delta_growth_2[i] = ((Y_2_series[i]) - \
                            (Y_2_series[i-1])) / \
                            (Y_2_series[i-1])     
    adjust_cost_1 = (delta_growth_1[:,2]-delta_growth_1[1,0])*100
    adjust_cost_2 = (delta_growth_2[:,2]-delta_growth_2[1,0])*100
    
    plt.figure(3)
    plt.figure(dpi=500)
    plt.rc('figure', figsize=[6.22,3.84])
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
    # Calculate the width of each bar group
    bar_width = 0.4
    # Create an array of x positions for each group of bars
    l = np.arange(len(years))
    plt.bar(l-0.2, adjust_cost_1[1:], width=bar_width, color = "green", label='OECD countries')
    plt.bar(l+0.2, adjust_cost_2[1:], width=bar_width, color = "orange", label='NOECD countries')
    plt.xlabel(r't', loc="right", labelpad=1, fontsize = 12)
    custom_ticks = np.arange(0, len(years), 3)  # Set ticks at every alternate year
    custom_labels = np.arange(2020, 2201, 30)  # Generate labels from 2030 to 2220 in steps of 20
    plt.xticks(custom_ticks, custom_labels, fontsize=10)
    plt.gca().xaxis.set_label_coords(1.025, 0.02)  # Set the position of the x-axis label
    plt.ylabel("\%", loc="top", rotation = 0, fontsize = 12)
    plt.gca().yaxis.set_label_coords(0.05, 0.99)
    plt.grid(axis ="y", which='both', color='gray', linestyle='--', linewidth=0.5)
    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['right'].set_visible(False)
    plt.legend(loc='upper left',  fontsize = 10, frameon = False, bbox_to_anchor=(0.05, 0.88)) 
    plt.savefig('GDP_change_LF_Opt.png', dpi=300, bbox_inches='tight')
    
    
    # 3. Climate damages in % of GDP
    plt.figure(4)
    plt.figure(dpi=500)
    plt.rc('figure', figsize=[6.22,3.84])
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
    plt.plot(years,D_1_series[:20,0]*100 , label="Laissez-faire", linewidth=2.5)
    plt.plot(years,D_1_series[:20,1]*100, label = "Efficient tax" , linewidth=2.5)
    plt.plot(years,D_1_series[:20,2]*100, label = "Optimal" , linewidth=2.5)
    plt.xlabel(r't', loc="right", labelpad=1, fontsize = 12)
    plt.xticks(np.arange(2030, 2220, step=25), fontsize= 12)
    plt.gca().xaxis.set_label_coords(1.025, 0.02)  # Set the position of the x-axis label
    plt.ylabel("$D_t$ [\% of GDP]", loc="top", rotation = 0, fontsize = 12)
    plt.yticks(np.arange(0,12,step = 2),fontsize=12)
    plt.gca().yaxis.set_label_coords(0.25, 0.88)
    plt.grid(axis ="y", which='both', color='gray', linestyle='--', linewidth=0.5)
    plt.xlim(2020, max(years))  
    plt.ylim(0,12)
    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['right'].set_visible(False)
    plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.1), 
               shadow=False, frameon = False, ncol=3, fontsize = 12) 
    plt.savefig('Climate_damages.png', dpi=300, bbox_inches='tight')
    
    
    # 4. Global mean surface temperature
    plt.figure(5)
    plt.figure(dpi=500)
    plt.rc('figure', figsize=[6.22,3.84])
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
    plt.plot(years,TEMP_series[:20,0] - TEMP_series[0,0] , label="Laissez-faire", linewidth=2.5)
    plt.plot(years,TEMP_series[:20,1]- TEMP_series[0,1], label = "Efficient tax" , linewidth=2.5)
    plt.plot(years,TEMP_series[:20,2]- TEMP_series[0,2], label = "Optimal" , linewidth=2.5)
    plt.xlabel(r't', loc="right", labelpad=1, fontsize = 12)
    plt.xticks(np.arange(2030, 2220, step=25), fontsize= 12)
    plt.gca().xaxis.set_label_coords(1.025, 0.02)  # Set the position of the x-axis label
    plt.ylabel("°Celsius", loc="top", rotation = 0, fontsize = 12)
    plt.yticks(np.arange(0,7,step = 1),fontsize=12)
    plt.gca().yaxis.set_label_coords(0.15, 0.9)
    plt.grid(axis ="y", which='both', color='gray', linestyle='--', linewidth=0.5)
    plt.xlim(2020, max(years))  
    plt.ylim(0,7)
    plt.axhline(y=0.89, color='r',linestyle='dotted', label="2 °C-target")
    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['right'].set_visible(False)
    plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.1), 
               shadow=False, frameon = False, ncol=3, fontsize = 12) 
    plt.savefig('Temperature.png', dpi=300, bbox_inches='tight')
    
    # 5. Subsidy on green R&D in % of GDP 
    # Convert the subsidy rate theta to subsidy in % of regional GDP
    sub_percent_1 = np.empty((20))
    sub_percent_2 = np.empty((20))
    
    sub_percent_1[0] =\
    (theta_1_series[0,2] * (eta_2_1_series[0,2] * (1-alpha_2)*alpha_2 \
    **((1+alpha_2)/(1-alpha_2))*p_2_1_series[0,2]*N_2_1_series[0,2]\
    *Q_2_1_init_0))/Y_1_series[0,2]
        
    sub_percent_2[0] =\
    (theta_2_series[0,2] * (eta_2_2_series[0,2] * (1-alpha_2)*alpha_2 \
    **((1+alpha_2)/(1-alpha_2))*p_2_2_series[0,2]*N_2_2_series[0,2]\
    *Q_2_2_init_0))/Y_2_series[0,2]
    
    for i in range(1,20):
        sub_percent_1[i] = \
        (theta_1_series[i,2] * (eta_2_1_series[i,2] * (1-alpha_2)*alpha_2 \
        **((1+alpha_2)/(1-alpha_2))*p_2_1_series[i,2]*N_2_1_series[i,2]\
        *Q_2_1_series[i-1,2]))/Y_1_series[i,2]
            
        sub_percent_2[i] = \
        (theta_2_series[i,2] * (eta_2_2_series[i,2] * (1-alpha_2)*alpha_2 \
        **((1+alpha_2)/(1-alpha_2))*p_2_2_series[i,2]*N_2_2_series[i,2]\
        *Q_2_2_series[i-1,2]))/Y_2_series[i,2]
    
    
    plt.figure(6)
    plt.figure(dpi=500)
    plt.rc('figure', figsize=[6.22,3.84])
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
    plt.plot(years,sub_percent_1[:20]*100 , label="OECD-countries", color = "green",linewidth=2.5)
    plt.plot(years,sub_percent_2[:20]*100, label = "NOECD-countries" ,color = "orange", linewidth=2.5)
    plt.xlabel(r't', loc="right", labelpad=1, fontsize = 12)
    plt.xticks(np.arange(2030, 2220, step=25), fontsize= 12)
    plt.gca().xaxis.set_label_coords(1.025, 0.02)  # Set the position of the x-axis label
    plt.ylabel("\% ", loc="top", rotation = 0, fontsize = 12)
    plt.yticks(np.arange(0,2.5,step = 0.5),fontsize=12)
    plt.gca().yaxis.set_label_coords(0.05, 0.9)
    plt.grid(axis ="y", which='both', color='gray', linestyle='--', linewidth=0.5)
    plt.xlim(2020, max(years))  
    plt.ylim(0,2.5)
    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['right'].set_visible(False)
    plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.1), 
               shadow=False, frameon = False, ncol=3, fontsize = 12) 
    plt.savefig('Subsidy.png', dpi=300, bbox_inches='tight')
    
    # 6. Efficient carbon tax
    plt.figure(7)
    plt.figure(dpi=500)
    plt.rc('figure', figsize=[6.22,3.84])
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
    plt.plot(years,tau_series[:20,0]*10**15/(10**9)*12/44 , label="Laissez-faire", linewidth=2.5)
    plt.plot(years,tau_series[:20,1]*10**15/(10**9)*12/44 , label="Efficient tax", linewidth=2.5)
    plt.plot(years,tau_series[:20,2]*10**15/(10**9)*12/44 , label="Optimal", linewidth=2.5)
    plt.xlabel(r't', loc="right", labelpad=1, fontsize = 12)
    plt.xticks(np.arange(2030, 2220, step=25), fontsize= 12)
    plt.gca().xaxis.set_label_coords(1.025, 0.02)  # Set the position of the x-axis label
    plt.ylabel("[$\$$/tCO$_2$]", loc="top", rotation = 0, fontsize = 12)
    plt.yticks(np.arange(0,500,step = 50),fontsize=12)
    plt.gca().yaxis.set_label_coords(0.15, 0.9)
    plt.grid(axis ="y", which='both', color='gray', linestyle='--', linewidth=0.5)
    plt.xlim(2020, max(years))  
    plt.ylim(-4,500)
    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['right'].set_visible(False)
    plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.1), 
               shadow=False, frameon = False, ncol=3, fontsize = 12) 
    plt.savefig('Efficient_tax.png', dpi=300, bbox_inches='tight')
    
    # 7. Regional GDP
    plt.figure(8)
    plt.figure(dpi=500)
    plt.rc('figure', figsize=[6.22,3.84])
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
    plt.plot(years,Y_1_series[:20,0] , label="Laissez-faire OECD", linewidth=1.5)
    plt.plot(years,Y_1_series[:20,2], label = "Optimal OECD" , linewidth=1.5, linestyle = "--")
    plt.plot(years,Y_2_series[:20,0], label = "Laissez-faire NOECD" , linewidth=1.5)
    plt.plot(years,Y_2_series[:20,2], label = "Optimal NOECD" , linewidth=1.5,linestyle = "--")
    plt.xlabel(r't', loc="right", labelpad=1, fontsize = 12)
    plt.xticks(np.arange(2030, 2240, step=25), fontsize= 12)
    plt.gca().xaxis.set_label_coords(1.025, 0.02)  # Set the position of the x-axis label
    plt.ylabel("Quadr. U.S.-$\$$", loc="top", rotation = 0, fontsize = 12)
    plt.yticks(np.arange(0,8,step = 1),fontsize=12)
    plt.gca().yaxis.set_label_coords(0.23, 0.9)
    plt.grid(axis ="y", which='both', color='gray', linestyle='--', linewidth=0.5)
    plt.xlim(2020, max(years))  
    plt.ylim(0,8)
    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['right'].set_visible(False)
    plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.1), 
               shadow=False, frameon = False, ncol=3, fontsize = 12) 
    plt.savefig('GDP.png', dpi=300, bbox_inches='tight')
    
    
if excel_output == "yes":
    data = list(zip(K_series[:25,0], K_series[:25,1], K_series[:25,2],
                    r_series[:25,0], r_series[:25,1], r_series[:25,2],
                    q_series[:25,0], q_series[:25,1], q_series[:25,2],
                    C_series[:25,0], C_series[:25,1], C_series[:25,2],
                    Z_series[:25,0], Z_series[:25,1], Z_series[:25,2],
                    S_series[:25,0], S_series[:25,1], S_series[:25,2],
                    S_1_series[:25,0], S_1_series[:25,1], S_1_series[:25,2],
                    S_2_series[:25,0], S_2_series[:25,1], S_2_series[:25,2],
                    R_series[:25,0], R_series[:25,1], R_series[:25,2],
                    X_1_series[:25,0], X_1_series[:25,1], X_1_series[:25,2],
                    X_2_series[:25,0], X_2_series[:25,1], X_2_series[:25,2],
                    Y_1_series[:25,0], Y_1_series[:25,1], Y_1_series[:25,2],
                    Y_2_series[:25,0], Y_2_series[:25,1], Y_2_series[:25,2],
                    E_1_1_series[:25,0], E_1_1_series[:25,1], E_1_1_series[:25,2],
                    E_1_2_series[:25,0], E_1_2_series[:25,1], E_1_2_series[:25,2],
                    E_2_1_series[:25,0], E_2_1_series[:25,1], E_2_1_series[:25,2],
                    E_2_2_series[:25,0], E_2_2_series[:25,1], E_2_2_series[:25,2],
                    p_1_1_series[:25,0], p_1_1_series[:25,1], p_1_1_series[:25,2],
                    p_1_2_series[:25,0], p_1_2_series[:25,1], p_1_2_series[:25,2],
                    p_2_1_series[:25,0], p_2_1_series[:25,1], p_2_1_series[:25,2],
                    p_2_2_series[:25,0], p_2_2_series[:25,1], p_2_2_series[:25,2],
                    N_0_1_series[:25,0], N_0_1_series[:25,1], N_0_1_series[:25,2],
                    N_1_1_series[:25,0], N_1_1_series[:25,1], N_1_1_series[:25,2],
                    N_2_1_series[:25,0], N_2_1_series[:25,1], N_2_1_series[:25,2],
                    N_0_2_series[:25,0], N_0_2_series[:25,1], N_0_2_series[:25,2],
                    N_1_2_series[:25,0], N_1_2_series[:25,1], N_1_2_series[:25,2],
                    N_2_2_series[:25,0], N_2_2_series[:25,1], N_2_2_series[:25,2],
                    w_1_series[:25,0], w_1_series[:25,1], w_1_series[:25,2],
                    w_2_series[:25,0], w_2_series[:25,1], w_2_series[:25,2],
                    D_1_series[:25,0], D_1_series[:25,1], D_1_series[:25,2],
                    D_2_series[:25,0], D_2_series[:25,1], D_2_series[:25,2],
                    TEMP_series[:25,0], TEMP_series[:25,1], TEMP_series[:25,2],
                    tau_series[:25,0], tau_series[:25,1], tau_series[:25,2],
                    chi_1_1_series[:25,0], chi_1_1_series[:25,1], chi_1_1_series[:25,2],
                    chi_2_1_series[:25,0], chi_2_1_series[:25,1], chi_2_1_series[:25,2],
                    chi_1_2_series[:25,0], chi_1_2_series[:25,1], chi_1_2_series[:25,2],
                    chi_2_2_series[:25,0], chi_2_2_series[:25,1], chi_2_2_series[:25,2],
                    eta_1_1_series[:25,0], eta_1_1_series[:25,1], eta_1_1_series[:25,2],
                    eta_2_1_series[:25,0], eta_2_1_series[:25,1], eta_2_1_series[:25,2],
                    eta_1_2_series[:25,0], eta_1_2_series[:25,1], eta_1_2_series[:25,2],
                    eta_2_2_series[:25,0], eta_2_2_series[:25,1], eta_2_2_series[:25,2],
                    M_1_1_series[:25,0], M_1_1_series[:25,1], M_1_1_series[:25,2],
                    M_2_1_series[:25,0], M_2_1_series[:25,1], M_2_1_series[:25,2],
                    M_1_2_series[:25,0], M_1_2_series[:25,1], M_1_2_series[:25,2],
                    M_2_2_series[:25,0], M_2_2_series[:25,1], M_2_2_series[:25,2], 
                    N_bar_1_series[:25,0], N_bar_1_series[:25,1], N_bar_1_series[:25,2],
                    N_bar_2_series[:25,0], N_bar_2_series[:25,1], N_bar_2_series[:25,2],
                    H_1_1_series[:25,0], H_1_1_series[:25,1], H_1_1_series[:25,2],
                    H_2_1_series[:25,0], H_2_1_series[:25,1], H_2_1_series[:25,2],
                    H_1_2_series[:25,0], H_1_2_series[:25,1], H_1_2_series[:25,2],
                    H_2_2_series[:25,0], H_2_2_series[:25,1], H_2_2_series[:25,2],
                    s_1_series[:25,0], s_1_series[:25,1], s_1_series[:25,2],
                    s_2_series[:25,0], s_2_series[:25,1], s_2_series[:25,2],
                    Q_1_1_series[:25,0], Q_1_1_series[:25,1], Q_1_1_series[:25,2],
                    Q_2_1_series[:25,0], Q_2_1_series[:25,1], Q_2_1_series[:25,2],
                    Q_1_2_series[:25,0], Q_1_2_series[:25,1], Q_1_2_series[:25,2],
                    Q_2_2_series[:25,0], Q_2_2_series[:25,1], Q_2_2_series[:25,2],
                    theta_1_series[:25,0], theta_1_series[:25,1], theta_1_series[:25,2],
                    theta_2_series[:25,0], theta_2_series[:25,1], theta_2_series[:25,2]))
    df_a = pd.DataFrame(data)
    df_a.columns = ("Capital LF", "Capital TAX", "Capital OPT",
                    "Interest rate LF", "Interest rate TAX", "Interest rate OPT", 
                    "Discount factor LF", "Discount factor TAX", "Discount factor OPT",
                    "Consumption LF", "Consumption TAX", "Consumption OPT", 
                    "Emissions LF", "Emissions TAX", "Emissions OPT", 
                    "Carbon concentration LF", "Carbon concentration TAX", "Carbon concentration OPT",
                    "Permanent carbon LF", "Permanent carbon TAX", "Permanent carbon OPT",
                    "Non-permanent carbon LF", "Non-permanent carbon TAX", "Non-permanent carbon OPT",
                    "Resource stock LF", "Resource stock TAX", "Resource stock OPT", 
                    "Resource OECD LF", "Resource OECD TAX", "Resource OECD OPT", 
                    "Resource NOECD LF", "Resource NOECD TAX", "Resource NOECD OPT",
                    "GDP OECD LF", "GDP OECD TAX", "GDP OECD OPT",
                    "GDP NOECD LF", "GDP NOECD TAX", "GDP NOECD OPT",
                    "Dirty energy OECD LF", "Dirty energy OECD TAX", "Dirty energy OECD OPT",
                    "Dirty energy NOECD LF", "Dirty energy NOECD TAX", "Dirty energy NOECD OPT",
                    "Clean energy OECD LF", "Clean energy OECD TAX", "Clean energy OECD OPT",
                    "Clean energy NOECD LF", "Clean energy NOECD TAX", "Clean energy NOECD OPT",
                    "Dirty energy price OECD LF", "Dirty energy price OECD TAX", "Dirty energy price OECD OPT",
                    "Dirty energy price NOECD LF", "Dirty energy price NOECD TAX", "Dirty energy price NOECD OPT",
                    "Clean energy price OECD LF", "Clean energy price OECD TAX", "Clean energy price OECD OPT",
                    "Clean energy price NOECD LF", "Clean energy price NOECD TAX", "Clean energy price NOECD OPT",
                    "Labor final OECD LF", "Labor final OECD TAX", "Labor final OECD OPT",
                    "Labor clean OECD LF", "Labor clean OECD TAX", "Labor clean OECD OPT",
                    "Labor dirty OECD LF", "Labor dirty OECD TAX", "Labor dirty OECD OPT",
                    "Labor final NOECD LF", "Labor final NOECD TAX", "Labor final NOECD OPT",
                    "Labor clean NOECD LF", "Labor clean NOECD TAX", "Labor clean NOECD OPT",
                    "Labor dirty NOECD LF", "Labor dirty NOECD TAX", "Labor dirty NOECD OPT",
                    "Wages OECD LF", "Wages OECD TAX", "Wages OECD OPT",
                    "Wages NOECD LF", "Wages NOECD TAX", "Wages NOECD OPT",
                    "Climate damages OECD LF", "Climate damages OECD TAX", "Climate damages OECD OPT",
                    "Climate damages NOECD LF", "Climate damages NOECD TAX", "Climate damages NOECD OPT",
                    "Temperature LF", "Temperature TAX", "Temperature OPT", 
                    "Tax rate LF", "Tax rate TAX", "Tax rate OPT",
                    "Dirty energy share OECD LF", "Dirty energy share OECD TAX", "Dirty energy share OECD OPT",
                    "Clean energy share OECD LF", "Clean energy share OECD TAX", "Clean energy share OCED OPT",
                    "Dirty energy share NOECD LF", "Dirty energy share NOECD TAX", "Dirty energy share NOECD OPT",
                    "Clean energy share NOECD LF", "Clean energy share NOECD TAX", "Clean energy share NOCED OPT", 
                    "Success probability clean OECD LF", "Success probability clean OECD TAX", "Success probability clean OECD OPT",
                    "Success probability dirty OECD LF", "Success probability dirty OECD TAX", "Success probability dirty OECD OPT",
                    "Success probability clean NOECD LF", "Success probability clean NOECD TAX", "Success probability clean NOECD OPT",
                    "Success probability dirty NOECD LF", "Success probability dirty NOECD TAX", "Success probability dirty NOECD OPT",
                    "Machines clean OECD LF", "Machines clean OECD TAX", "Machines clean OECD OPT",
                    "Machines dirty OECD LF", "Machines dirty OECD TAX", "Machines dirty OECD OPT",
                    "Machines clean NOECD LF", "Machines clean NOECD TAX", "Machines clean NOECD OPT",
                    "Machines dirty NOECD LF", "Machines dirty NOECD TAX", "Machines dirty NOECD OPT",
                    "Population OECD", "Population OECD", "Population OECD",
                    "Population NOECD", "Population NOECD", "Population NOECD",
                    "R&D dirty OECD LF", "R&D dirty OECD TAX", "R&D dirty OECD OPT",
                    "R&D clean OECD LF", "R&D clean OECD TAX", "R&D clean OECD OPT",
                    "R&D dirty NOECD LF", "R&D dirty NOECD TAX", "R&D dirty NOECD OPT",
                    "R&D clean NOECD LF", "R&D clean NOECD TAX", "R&D clean NOECD OPT",
                    "Scientists OECD LF", "Scientists OECD TAX", "Scientists OECD OPT",
                    "Scientists NOECD LF", "Scientists NOECD TAX", "Scientists NOECD OPT",
                    "Quality dirty OECD LF", "Quality dirty OECD TAX", "Quality dirty OECD OPT",
                    "Quality clean OECD LF", "Quality clean OECD TAX", "Quality clean OECD OPT",
                    "Quality dirty NOECD LF", "Quality dirty NOECD TAX", "Quality dirty NOECD OPT",
                    "Quality clean NOECD LF", "Quality clean NOECD TAX", "Quality clean NOECD OPT",
                    "Subsidy OECD LF", "Subsidy OECD TAX", "Subsidy OECD OPT",
                    "Subsidy NOECD LF", "Subsidy NOECD TAX", "Subsidy NOECD OPT")
    
    parameters = [phi, phi_L, phi_0, gamma_1, gamma_2, S_bar, alpha_0, nu_0,
                  Q_0_1, Q_0_2, kappa, rho, nu_1, Q_1_1_init_0, Q_1_2_init_0,
                  Q_2_1_init_0, Q_2_2_init_0, alpha_2, zeta, c, R_0, beta, sigma, 
                  N_bar_1_0, N_bar_2_0, g, K_bar_1_0, K_bar_2_0, S_1_0, S_2_0,
                  lam, iota, omega
                  ]
    pnames = ('phi', 'phi_L', 'phi_0', 'gamma_1', 'gamma_2', 'S_bar',
              'alpha_0', 'nu_0', 'Q_0_1', 'Q_0_2', 'kappa', 'rho',
              'nu_1', 'Q_1_1_init', 'Q_1_2_init', 'Q_2_1_init', 'Q_2_2_init',
              'alpha_2', 'zeta', 'c', 'R_0', 'beta', 'sigma', 'N_bar_1_0', 
              'N_bar_2_0', 'g', 'K_bar_1_0', 'K_bar_2_0', 'S_1_0', 'S_2_0',
              "lambda", "iota", "omega")
    df_b = pd.DataFrame(list(zip(pnames, parameters)))    
    path = "Master_thesis_Tim_Luediger_Output.xlsx"
    book = load_workbook(path)
    writer = pd.ExcelWriter(path, engine = 'openpyxl')
    writer.book = book
    df_a.to_excel(writer, sheet_name = 'Results', startrow = 1, startcol = 0)
    df_b.to_excel(writer, sheet_name = 'Parameters', startrow = 1,
                  startcol = 0, header = False, index = False)
    writer.save()
    writer.close()
        



















