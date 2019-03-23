# does the calculations for the seondary cooling loop
# would be better if we add a second turbine

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
nT = 0.85  # turbine efficiency
nP = 0.85  # pump efficiency
rho = 997  # kg/m3
# rho = 1000
T = 20  # celsius

turbine_powers = np.linspace(35, 85, 6)  # mw
# pressures = np.linspace(2, 10, 9)  #  bar
pressures = np.linspace(10, 110, 100)  # bar
flow_tab = []
pump_tab = []
qcore_tab = []
efficiency_tab = []
# mass_flows = []
# going to first calculate with isentropic pump and turbine

for power in turbine_powers:
    mass_flows = [power]
    pump_powers = [power]
    core_heats = [power]
    eff = [power]
    power *= 10**3  # k watts
    for P in pressures:
        h3 = steamTable.hV_p(P)
        s3 = steamTable.sV_p(P)
        s4s = s3

        sg_atm = steamTable.sV_p(atm)
        sf_atm = steamTable.sL_p(atm)
        sfg_atm = sg_atm-sf_atm

        X4s = (s4s - sf_atm) / sfg_atm

        hf_atm = steamTable.hL_p(atm)
        hg_atm = steamTable.hV_p(atm)
        hfg_atm = hg_atm - hf_atm

        h4s = hf_atm + X4s * hfg_atm
        # h4 = h3 - nT*(h3 - h4s)
        h4 = h4s  # neglecting turbine inefficienct here

        m_dot = power / (h3 - h4)  # kg/s

        mass_flows.append(m_dot)

        DeltaP = (P - atm) * 100000  # pa

        pump_power = DeltaP * m_dot / rho
        #  pump_power = DeltaP * m_dot / rho / nP
        pump_power *= 10 ** -3  # kpa
        pump_powers.append(pump_power)

        s1 = steamTable.s_pt(atm, T)
        h1 = steamTable.h_pt(atm, T)
        h2 = steamTable.h_ps(P, s1)

        Q_core = m_dot*(h3-h2) * 10 ** -3  # mw
        core_heats.append(Q_core)

        efficiency = power / Q_core * 10 ** -3
        eff.append(efficiency)

    efficiency_tab.append(eff)
    qcore_tab.append(core_heats)
    flow_tab.append(mass_flows)
    pump_tab.append(pump_powers)

# print(tabu(flow_tab, tablefmt="github"))
# print(tabu(pump_tab, tablefmt="github"))
# print(tabu(qcore_tab, tablefmt="github"))

# myDict = {65: [1, 2, 3], 75: [4, 5, 6]}
# tuple = (85, [7, 8, 9])
# myDict[2] = dict(tuple)
# print(myDict[2])
# # print(myDict[85])
# print(myDict[65])
#
# print(qcore_tab[0])
# print(len(qcore_tab[0]))
# print(len(pressures))
# print(qcore_tab[0][1:len(qcore_tab[0])+1])
# plt.plot()
# print(nrg_axis)
# print(len(nrg_axis))
# print(rad_axis)

for n in range(len(flow_tab)):
    la = str(flow_tab[n][0]) + 'MWe'
    plt.plot(pressures, flow_tab[n][1:len(flow_tab[n])+1], label=la)
plt.xlabel("High Pressure [Bar]")
plt.ylabel("Mass Flow Rate of Coolant [Kg/s]")
plt.legend()
plt.show()

for n in range(len(pump_tab)):
    la = str(pump_tab[n][0]) + 'MWe'
    plt.plot(pressures, pump_tab[n][1:len(pump_tab[n])+1], label=la)
plt.xlabel("High Pressure [Bar]")
plt.ylabel("Pump Power [KW]")
plt.legend()
plt.show()

for n in range(len(qcore_tab)):
    la = str(qcore_tab[n][0]) + 'MWe'
    plt.plot(pressures, qcore_tab[n][1:len(qcore_tab[n])+1], label=la)
plt.xlabel("High Pressure [Bar]")
plt.ylabel("Core Heat [MWth]")
plt.legend()
plt.show()

for n in range(len(efficiency_tab)):
    la = str(efficiency_tab[n][0]) + 'MWe'
    plt.plot(pressures, efficiency_tab[n][1:len(efficiency_tab[n])+1], label=la)
    # plt.plot(pressures, efficiency_tab[n][1:len(efficiency_tab[n])+1])
plt.xlabel("High Pressure [Bar]")
plt.ylabel("Wt/Qc")
# plt.legend()
plt.show()
