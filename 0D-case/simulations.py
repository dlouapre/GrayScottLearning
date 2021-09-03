#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt

def simulation(environment, network, N, dt):
    
    Binit = environment['B_init']
    B_tau = environment['B_tau']    
    
    bolus_period = environment['bolus_period']
    bolus_delta = environment['bolus_delta']
    S_bolus = environment['S_bolus']    
    S_tau = environment['S_tau']
    kT = environment['kT']
    T_tau = environment['T_tau']
    shift_time = environment['shift_time']

    kNT = environment['kNT']
    kBT = environment['kBT']
    N_tau = environment['N_tau']
    kBN = environment['kBN']
        
    M_tau = network['M_tau']
    L_tau = network['L_tau']    
    
    kM = network['kM']
    kL = network['kL']
    kA = network['kA']
    
    kP = network['kP']
    kD = network['kD']
    
    # Initialisation
    Ss = np.zeros(N)
    Ts = np.zeros(N)
    Ms = np.zeros(N)
    Ls = np.zeros(N)    
    Ns = np.zeros(N)
    Bs = np.zeros(N)
    
    kT_effs = np.zeros(N)
    
    Bs[0] = Binit
    
    for n in range(1,N):
    
        S, T, M, L, N, B = Ss[n - 1], Ts[n - 1], Ms[n - 1], Ls[n - 1], Ns[n - 1], Bs[n - 1]
        
        newS, newT, newM, newL, newN, newB = S, T, M, L, N, B
        
        # 1) Apply boluses

        bolus_period_n = bolus_period / dt
        bolus_delta_n = bolus_delta / dt
        
        # 1a) Stimulus bolus
        if(n%bolus_period_n==0):                
            newS = S_bolus
        
        # 1b) Toxin bolus
        # Evolution of kT_eff from time            
        if (n*dt < shift_time):            
            kT_eff = kT * (n*dt/shift_time)
        elif(n*dt < 2 * shift_time):
            kT_eff = kT * (2 - n*dt/shift_time)
        else:
            kT_eff = 0
            
        kT_effs[n] = kT_eff

        if(n%bolus_period_n == bolus_delta_n):
            newT += kT_eff * S_bolus
        
        
        # 2) B growth (with B0 max capacity)
        B0 = 1
        newB += B * (B0-B) * dt/B_tau
        
        # Apply decays
        newS -= S * dt / S_tau
        newT -= T * dt / T_tau
        newM -= M * dt / M_tau
        newL -= L * dt / L_tau
        newN -= N * dt / N_tau   
                            
        # Short term memory production
        # S => S + M
        newM += kM * S * dt

        # Longterm memory production
        # T + M => T + L
        newL += kL * T * M * dt
        newM -= kL * T * M * dt
        
        # Antidote production
        amount = 0
        # 1) DIRECT B + T => N + B + T        
        amount += kD * B * T * dt
        # 2) PRE-EMPTIVE B + T => O + B + T
        amount += kP * B * S * dt
        # 3) ASSOCIATIVE B + T => O + B + T
        amount += kA * B * L * S * dt
                
        
        MAX_RATE = 0.1
        if amount > MAX_RATE*dt:    
            #amount = MAX_RATE*dt
            print("MAX RATE !",end=None)

        newN += amount
        
        # Toxin degradation by antidote T + N => 
        newN -= kNT * T * N * dt
        newT -= kNT * T * N * dt
        
        # Toxin damages B
        newB -= kBT * T * B * dt
        newT -= kBT * T * B * dt
        
        # Antidote damages B
        newB -= kBN * N * B * dt
        
        # Register new values, remove negative
        Ss[n] = max(0, newS)
        Ts[n] = max(0, newT)
        Ms[n] = max(0, newM)        
        Ls[n] = max(0, newL)    
        Ns[n] = max(0, newN)
        Bs[n] = max(0, newB)

    return Ss, Ts, Ms, Ls, Ns, Bs, kT_effs




def makeplot(N, dt, Is, Ts, SIs, Ls, Os, Bs, YMAX = 1.25, plot_kT = False):
    
    Is = Is[:N]
    Ts = Ts[:N]
    SIs = SIs[:N]
    Ls = Ls[:N]
    Os = Os[:N]
    Bs = Bs[:N]
    
    ts = np.arange(0,dt*N,dt)
    
    
    
    plt.plot(ts,Is,color="navajowhite", label = "Stimulus S", ls='--')
    
    plt.plot(ts,Ts,color = "red", label="Toxin T")
    
    if plot_kT:
        
        kT = 1
        shift_time = 10000
        
        kT_effs = np.zeros(N)
        for n in range(N):
            if (n*dt < shift_time):            
                kT_effs[n] = kT * (n*dt/shift_time)
            elif(n*dt < 2 * shift_time):
                kT_effs[n] = kT * (2 - n*dt/shift_time)
            else:
                kT_effs[n] = 0
        
        plt.plot(ts,kT_effs, color = "black", ls = '-', label="Environment $\epsilon$")
   

    #plt.plot(ts, SIs,ls='--', color="dodgerblue", label = "$S_I$")
    
    if sum(Ls) > 0:
        plt.plot(ts,Ls,color="blue",label="Long term memory L")

    if sum(Os) > 0:
        plt.plot(ts,Os,color='forestgreen',label="Antidote N")
        
    if sum(Bs) > 0:
        plt.plot(ts,Bs,color="purple", label = "Target B", linewidth=3)
    
    plt.xlim(0,N*dt)
    plt.ylim(0,YMAX)
    plt.legend(loc="upper right")
    plt.xlabel("Time")
    plt.ylabel("Concentration")
    plt.tight_layout()
    
    
def linesearch(environment, network, N, dt, var, minvlog, maxvlog, stepvlog, doplot = False):
    
    v_list = 10**np.arange(minvlog , maxvlog + stepvlog, stepvlog)
    target_list = []
    best_v = 0
    max_target = 0
    
    for v in v_list:    
        network[var] = v
            
        Is, Ts, SIs, Ls, Os, Bs, _ = simulation(environment, network, N, dt)
        
        target = np.mean(Bs)
        target_list.append(target)
        
        if(target > max_target):
            max_target = target
            best_v = v
            print(v, target)
    
    return best_v, v_list, np.array(target_list)
    
        

