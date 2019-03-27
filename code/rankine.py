# does the calculations for the seondary cooling loop
# adding the second turbine does not have the anticipated
# benefits, might be an error in this code


import math as math
import numpy as np
from pyXSteam.XSteam import XSteam
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import logging
import scipy.integrate as inte
from tabulate import tabulate as tabu

# steamTable = XSteam(XSteam.UNIT_SYSTEM_BARE)

steamTable = XSteam(XSteam.UNIT_SYSTEM_MKS)

# print(steamTable.hV_p(1.01325))

atm = 1.01325  # bar
# atm = 3
# nT = 0.9  # turbine efficiency
nT = 1
# nP = 0.9  # pump efficiency
nP = 1
rho = 997  # kg/m3
# rho = 1000
T = 20  # celsius
# max_press = 60  # 60 bar is close to the operating press of real nuke systems
max_press = 100 # maybe we can cheat a bit cause mdot is gonna be lower
h6 = steamTable.hV_p(max_press)
s6 = steamTable.sV_p(max_press)
h5 = steamTable.hL_p(max_press)

turbine_powers = np.linspace(35, 85, 6)  # mw
# pressures = np.linspace(2, 10, 9)  #  bar
pressures = np.linspace(20, 60, 100)  # bar, holds the low pressure turbine val
flow_tab_low = []
flow_tab_med = []
flow_tab_hig = []
flow_tab = []
pump_tab = []
qcore_tab = []
efficiency_tab = []
# mass_flows = []
# going to first calculate with isentropic pump and turbine

for power in turbine_powers:
    mass_flows_low = [power]
    mass_flows_med = [power]
    mass_flows_hig = [power]
    pump_powers_low = [power]
    pump_powers_hig = [power]
    pump_powers = [power]
    core_heats_low = [power]
    core_heats_hig = [power]
    core_heats = [power]
    eff = [power]
    power *= 10**3  # k watts
    for P in pressures:  # iterates over medium pressures
        h7 = steamTable.hV_p(P)
        s7 = steamTable.sV_p(P)
        s8s = steamTable.sV_p(P)

        h3 = steamTable.hL_p(P)
        h7s = steamTable.h_ps(P, s6)
        h7p = h6 - nT * (h6 - h7s)

        sg_atm = steamTable.sV_p(atm)
        sf_atm = steamTable.sL_p(atm)
        sfg_atm = sg_atm-sf_atm

        X8s = (s8s - sf_atm) / sfg_atm

        hf_atm = steamTable.hL_p(atm)
        hg_atm = steamTable.hV_p(atm)
        hfg_atm = hg_atm - hf_atm

        h8s = hf_atm + X8s * hfg_atm
        h8 = h7 - nT * (h7 - h8s)
        # h4 = h4s  # neglecting turbine inefficienct here

        a = np.array([[abs(h7p - h3), h3 - h7], [h7p - h6, abs(h8 - h7)]])
        b = np.array([0, power / nT])
        m_dots = np.linalg.solve(a, b)  # high, low
        m_dot_hig = m_dots[0]
        m_dot_low = m_dots[1]
        m_dot_med = m_dot_low - m_dot_hig

        mass_flows_low.append(m_dot_low)
        mass_flows_hig.append(m_dot_hig)
        mass_flows_med.append(m_dot_med)

        DeltaP_low = (P - atm) * 100000  # pa
        DeltaP_hig = (max_press - P) * 100000  # pa

        # pump_power = DeltaP * m_dot / rho
        pump_power_low = DeltaP_low * m_dot_low / rho / nP
        pump_power_low *= 10 ** -3  # kpa
        pump_powers_low.append(pump_power_low)

        pump_power_hig = DeltaP_hig * m_dot_hig / rho / nP
        pump_power_hig *= 10 ** -3  # kpa
        pump_powers_hig.append(pump_power_hig)

        pump_power = pump_power_low + pump_power_hig
        pump_powers.append(pump_power)

        s1 = steamTable.s_pt(atm, T)
        h1 = steamTable.h_pt(atm, T)
        h2 = h1 + pump_power_low / m_dot_low * nP
        # h2 = steamTable.h_ps(P, s1)
        h4 = h3 + pump_power_hig / m_dot_hig * nP

        Q_core_low = m_dot_low * (h3 - h2) * 10 ** -3 # mw
        core_heats_low.append(Q_core_low)
        Q_core_hig = m_dot_hig * (h6 - h4) * 10 ** -3 #mw
        core_heats_hig.append(Q_core_hig)
        Q_core = Q_core_hig + Q_core_low
        core_heats.append(Q_core)

        efficiency = power / Q_core * 10 ** -3
        eff.append(efficiency)

    efficiency_tab.append(eff)
    qcore_tab.append(core_heats)
    flow_tab.append(mass_flows_low)
    pump_tab.append(pump_powers)

for n in range(len(flow_tab)):
    la = str(flow_tab[n][0]) + 'MWe'
    plt.plot(pressures, flow_tab[n][1:len(flow_tab[n])+1], label=la)
plt.xlabel("Middle Pressure [Bar]")
plt.ylabel("Mass Flow Rate of Coolant [Kg/s]")
plt.legend()
plt.show()

for n in range(len(pump_tab)):
    la = str(pump_tab[n][0]) + 'MWe'
    plt.plot(pressures, pump_tab[n][1:len(pump_tab[n])+1], label=la)
plt.xlabel("Middle Pressure [Bar]")
plt.ylabel("Pump Power [KW]")
plt.legend()
plt.show()

for n in range(len(qcore_tab)):
    la = str(qcore_tab[n][0]) + 'MWe'
    plt.plot(pressures, qcore_tab[n][1:len(qcore_tab[n])+1], label=la)
plt.xlabel("Middle Pressure [Bar]")
plt.ylabel("Core Heat [MWth]")
plt.legend()
plt.show()

for n in range(len(efficiency_tab)):
    la = str(efficiency_tab[n][0]) + 'MWe'
    plt.plot(pressures, efficiency_tab[n][1:len(efficiency_tab[n])+1],
             label=la)
    # plt.plot(pressures, efficiency_tab[n][1:len(efficiency_tab[n])+1])
plt.xlabel("Middle Pressure [Bar]")
plt.ylabel("Wt/Qc")
# plt.legend()
plt.show()
