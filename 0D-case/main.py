#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt

from simulations import simulation, makeplot, linesearch

network0 = {'M_tau' : np.inf,
            'L_tau' : np.inf,            
            'kM' : 0,
            'kL' : 0,
            'kD' : 0,
            'kA' : 0,
            'kP' : 0}

networkD = {'M_tau' : np.inf,
            'L_tau' : np.inf,            
            'kM' : 0,
            'kL' : 0,
            'kD' : 1,
            'kA' : 0,
            'kP' : 0}

networkP = {'M_tau' : np.inf,
            'L_tau' : np.inf,            
            'kM' : 0,
            'kL' : 0,
            'kD' : 0,
            'kA' : 0,
            'kP' : 1}

networkA = {'M_tau' : 200,
            'L_tau' : 4000,            
            'kM' : 0.01,
            'kL' : 0.5,
            'kD' : 0,
            'kA' : 1,
            'kP' : 0}

environment = { 'B_init' : 1,
                'B_tau' : 500,
                'bolus_period' : 1000,
                'bolus_delta' : 400,
                'S_bolus' : 1,
                'S_tau' : 200,
                'T_tau' : 200,
                'N_tau' : 5000,
                'kT' : 1,
                'kBT' : 0.5,
                'kBN' : 0.001,
                'kNT' : 1,    
                'shift_time':10000}

environment0 = environment.copy()
environment0['B_init'] = 0

dt = 1    
N = 40000


plt.figure(figsize=(14,2.5))
Is, Ts, Ms, Ls, Ns, Bs, kT_effs = simulation(environment0, network0, N, dt)
makeplot(N, dt, Is, Ts, Ms, Ls, Ns, Bs, plot_kT = True)
plt.tight_layout()
plt.savefig("0D-environment.eps")

plt.figure(figsize=(14,2.5))
Is, Ts, Ms, Ls, Ns, Bs, kT_effs = simulation(environment, network0, N, dt)
makeplot(N, dt, Is, Ts, Ms, Ls, Ns, Bs)
plt.tight_layout()
plt.savefig("0D-no_network.eps")


print("Search best kD")
best_kD, kD_list, kD_targets = linesearch(environment, networkD, N, dt, 'kD', -3, -1, 0.1)
networkD['kD'] = best_kD

print("Search best kP")
best_kP, kP_list, kP_targets = linesearch(environment, networkP, N, dt, 'kP', -3, -1, 0.1)
networkP['kP'] = best_kP

print("Search best kA")
best_kA, kA_list, kA_targets = linesearch(environment, networkA, N, dt, 'kA', -3, -1, 0.1)
networkA['kA'] = best_kA

plt.figure()
plt.semilogx(kD_list, kD_targets, label = "Direct", ls = "--", marker="o", color = "black")
plt.semilogx(kP_list, kP_targets, label = "Pre-emptive", ls = "--", marker="P", color = "lightgrey")
plt.semilogx(kA_list, kA_targets, label = "Associative", ls = "--", marker="^", color = "grey")
plt.legend(loc = "upper left")
plt.xlabel("Reaction constant")
plt.ylabel("Average B concentration")
plt.ylim(0,1)
plt.savefig("0D-optimisation.eps")


FIGSIZE = (21,10)
plt.figure(figsize=FIGSIZE)

plt.subplot(3,1,1)
Is, Ts, Ms, Ls, Ns, Bs, kT_effs = simulation(environment, networkD, N, dt)
makeplot(N, dt, Is, Ts, Ms, Ls, Ns, Bs)
plt.title("Direct network")

plt.subplot(3,1,2)
Is, Ts, Ms, Ls, Ns, Bs, kT_effs = simulation(environment, networkP, N, dt)
makeplot(N, dt, Is, Ts, Ms, Ls, Ns, Bs)
plt.title("Pre-emptive network")

plt.subplot(3,1,3)
Is, Ts, Ms, Ls, Ns, Bs, kT_effs = simulation(environment, networkA, N, dt)
makeplot(N, dt, Is, Ts, Ms, Ls, Ns, Bs)
plt.title("Associative Learning network")

plt.savefig("0D-3_networks.eps")

plt.figure(figsize=(14,2.5))
plt.plot(np.arange(0, N * dt, dt), kT_effs, color = 'black', label = 'Environment $\epsilon$')
plt.plot(np.arange(0, N * dt, dt), Ls, color = 'blue', label = 'Long term memory L')
plt.xlabel("Time")
plt.ylabel("Concentration")
plt.legend(loc="upper left")
plt.savefig("0D-L_vs_epsilon.eps")




