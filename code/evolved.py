#tabulating many core parameters

import math as m
import numpy as np 
import matplotlib.pyplot as plt
from tabulate import tabulate as tabu

# constamts amd options

avogadro = 6.022 * 10 ** 23
sigma_f_235 = 2.156 # barns
sigma_a_238 = 0.404 # barns
sigma_a_235 = 2.844 # barns
nu = 5.297 / sigma_f_235
sigma_t_235 = 8.246
sigma_t_238 = 8.181
density_uo2 = 10970 #kg/m3
efficiency = 0.3 

power_options = np.linspace(35, 85, 11) #options in megawatts
power_electric = power_options[0] * 10**6 #watts
power_thermal = power_electric / efficiency
print(power_options)
# print(power_electric)
lifetime_options = np.linspace(1, 10, 10) #years
years = 2
lifetime = lifetime_options[years - 1] * 3.15 * 10 ** 7 #seconds 
print(lifetime_options)
print(lifetime)
# print(power_options)
# print(power_electric)
energy = power_thermal * lifetime # joules
print(energy)

# compute enrichment such that kinf = 1

beta = sigma_a_238 / 238
lamda = nu * sigma_f_235 / sigma_a_235
alpha = sigma_a_235 / 235 * (lamda + 1)
gamma = beta / alpha
enrichment = gamma / (1 + gamma)
enrichment_percent = enrichment*100

# compute u counts per m of fuel

m_u_per_m_uo2 = (238 * (1 - enrichment) + 235 * enrichment) / (238 * (1 - enrichment) + 235 * enrichment + 32)
n_u5_per_g_uo2 = enrichment * m_u_per_m_uo2 * avogadro / 235

# compute inner core volume

def vol_in(E): #wants e in joules
	fissions = E / (1.6 * 10 ** (-13)) / 200 
	fuel_mas = fissions / n_u5_per_g_uo2 * 10 ** (-3) #kg
	fuel_vol = fuel_mas / density_uo2 #m3
	fuel_rad = (fuel_vol / m.pi) ** (1/3) #m
	return(fuel_vol, fuel_rad)

def mdot_in(E, P): #wants e, p in joules, watts
	R = vol_in(E)[1] 
	H = R
	# H = 4/3*R
	r = 1/200 #radius of fuel rod in m
	cp = 130 #j/kg/k
	delta_t = 500 #k
	mdot = P / 200 * r / H / cp / delta_t
	return(mdot) 

# core geometry plotting

nrg_axis = []
vol_axis = []
rad_axis = []

for P in power_options:
	P *= 10 ** 6 #watts
	P /= efficiency #watts-thermal
	for T in lifetime_options: #years
		T *= 3.15 * 10 ** 7 #seconds
		E = P*T #joules
		e = E / 3.15 / (10**13) * efficiency #mwy
		nrg_axis.append(e) #Mw * years
		vol_axis.append(vol_in(E)[0])
		rad_axis.append(vol_in(E)[1])

# print(nrg_axis)
# print(len(nrg_axis))
# print(rad_axis)

# plt.plot(nrg_axis, rad_axis)
# plt.xlabel("Core energy, P*T [Mw*Y]")
# plt.ylabel("Core Volume [m3]")
# plt.show()

fig, ax1 = plt.subplots()
# t = np.arange(0.01, 10.0, 0.01)
# s1 = np.exp(t)
ax1.plot(nrg_axis, vol_axis, 'b-')
ax1.set_xlabel('Core energy [Mwe*Y]')
# Make the y-axis label, ticks and tick labels match the line color.
ax1.set_ylabel('Volume [m3]', color='b')
ax1.tick_params('y', colors='b')

ax2 = ax1.twinx()
# s2 = np.sin(2 * np.pi * t)
ax2.plot(nrg_axis, rad_axis, 'r.')
ax2.set_ylabel('Radius [m]', color='r')
ax2.tick_params('y', colors='r')

fig.tight_layout()
plt.show()

# print(vol_in(energy))

# row = []
tab = []

# cooling rate tabulation

for T in lifetime_options:
	row = [T]
	T *= 3.15 * 10 ** 7 #seconds
	for P in power_options: #years
		P *= 10 ** 6 #watts
		P /= efficiency #watts-thermal
		E = P*T
		row.append(mdot_in(E, P))
	tab.append(row)

print(tabu(tab, tablefmt="github"))