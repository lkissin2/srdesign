# tabulating many core parameters

import math as m
import numpy as np
import matplotlib.pyplot as plt
from tabulate import tabulate as tabu
import scipy.integrate as integrate
import scipy.special as special

# constamts amd options

avogadro = 6.022 * 10 ** 23
sigma_f_235 = 2.156 # barns
nu_sigma_f_235 = 5.297
nu_sigma_f_238 = 0.142
sigma_f_238 = 0.051
nu_238 = 0.142 / sigma_f_238
sigma_a_238 = 0.404 # barns
sigma_a_235 = 2.844 # barns
nu_235 = 5.297 / sigma_f_235
sigma_t_235 = 8.246
sigma_t_238 = 8.181
density_uo2 = 10970 #kg/m3
# density_uo2 = 20000 #kg/m3
efficiency = 0.3
rod_radius = 0.5*10**-2

power_options = np.linspace(35, 85, 11) #options in megawatts
power_electric = power_options[0] * 10**6 #watts
power_thermal = power_electric / efficiency
print(power_options)
# print(power_electric)

lifetime_options = np.linspace(1, 10, 5) #years
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
lamda = nu_235 * sigma_f_235 / sigma_a_235
alpha = sigma_a_235 / 235 * (lamda + 1)
gamma = beta / alpha
enrichment = gamma / (1 + gamma)
# print('en = ' + str(enrichment))
enrichment_percent = enrichment*100

# compute u counts per m of fuel

m_u_per_m_uo2 = (238 * (1 - enrichment) + 235 * enrichment) / (238 * (1 - enrichment) + 235 * enrichment + 32)
n_u5_per_g_uo2 = enrichment * m_u_per_m_uo2 * avogadro / 235
n_u_8_per_g_uo2 = (1- enrichment) * m_u_per_m_uo2 * avogadro / 238

# compute inner core volume

def vol_in(E): #wants e in joules
	fissions = E / (1.6 * 10 ** (-13)) / 200
	fuel_mas = fissions / n_u5_per_g_uo2 * 10 ** (-3) #kg
	fuel_vol = fuel_mas / density_uo2 #m3
	fuel_rad = (fuel_vol / m.pi) ** (1/3) #m
	return(fuel_vol, fuel_rad)

# print(vol_in(85*10**6*10*3.15*10**7))

def mdot_in(E, P): #wants e, p in joules, watts
	R = vol_in(E)[1]
	H = R
	# H = 4/3*R
	r = rod_radius #radius of fuel rod in m
	N = R / r
	cp = 130 #j/kg/k
	delta_t = 500 #k
	mdot = P / N / H / cp / delta_t
	return(mdot)

#outer core

def rad_out(wout):
	# density_uo2 = 10970
	wout /= 100
	m_u_per_m_uo2_out = (238 * (1 - wout) + 235 * wout) / (238 * (1 - wout) + 235 * wout + 32)
	n_u5_per_g_uo2_out = wout * m_u_per_m_uo2_out * avogadro / 235
	# Sigma_f_out = density_uo2 * n_u5_per_g_uo2_out * 1000 * sigma_f_235*10**-28
	n_u_8_per_g_uo2_out = (1 - wout) * m_u_per_m_uo2_out * avogadro / 238
	nu_macroscopic_fission_cross_section = density_uo2 * (n_u5_per_g_uo2_out * 1000 * nu_sigma_f_235*10**-28 + n_u_8_per_g_uo2_out * 1000 * nu_sigma_f_238*10**-28)
	macroscopic_absorption_cross_section = density_uo2 * (n_u5_per_g_uo2_out * 1000 * sigma_a_235*10**-28 + n_u_8_per_g_uo2_out * 1000 * sigma_a_238*10**-28)
	macroscopic_transport_cross_section =  density_uo2 * (n_u5_per_g_uo2_out * 1000 * sigma_t_235*10**-28 + n_u_8_per_g_uo2_out * 1000 * sigma_t_238*10**-28)
	diffusion_coefficient = 1 / 3 / macroscopic_transport_cross_section
	buckling = (nu_macroscopic_fission_cross_section - macroscopic_absorption_cross_section) / diffusion_coefficient
	# print('w = ' + str(wout))
	# print('b= ' + str(buckling))
	outer_radius = ( (2.405**2 + m.pi**2) / buckling ) ** (1/2)
	# outer_vol = m.pi*(outer_radius**3)
	# return outer_vol
	return outer_radius

for w in np.linspace(19.5, 99.5, 11):
	print('{0}, {1}'.format(w, rad_out(w)))

def pow_out(wout, E, P): #wants percent, mj, and mw
	wout /= 100
	P *= 10 ** 6
	E *= 10 ** 6
	R = rad_out(wout*100)
	rin = vol_in(E)[1]
	macroscopic_fission_cross_section_in = density_uo2 * (n_u5_per_g_uo2 * 1000 * sigma_f_235 * 10**-28 + n_u_8_per_g_uo2 * 1000 * sigma_f_238*10**-28)
	# print('macroscopic_fission_cross_section_in ' + str(macroscopic_fission_cross_section_in))
	energy_per_fission = 200 * 1.6 * 10 ** -13 #  joules
	phi_0 = P / vol_in(E)[0] / macroscopic_fission_cross_section_in / energy_per_fission
	# print(phi_0)
	m_u_per_m_uo2_out = (238 * (1 - wout) + 235 * wout) / (238 * (1 - wout) + 235 * wout + 32)
	# print('m_u_per_m_uo2_out ' + str(m_u_per_m_uo2_out))
	n_u5_per_g_uo2_out = wout * m_u_per_m_uo2_out * avogadro / 235
	# print('n_u5_per_g_uo2_out ' + str(n_u5_per_g_uo2_out))
	n_u_8_per_g_uo2_out = (1 - wout) * m_u_per_m_uo2_out * avogadro / 238
	# print('n_u_8_per_g_uo2_out ' + str(n_u_8_per_g_uo2_out))
	# Sigma_f_out = density_uo2 * n_u5_per_g_uo2_out * 1000 * sigma_f_235*10**-28
	macroscopic_fission_cross_section_out = density_uo2 * (n_u5_per_g_uo2_out * sigma_f_235 * 1000 * 10**(-28) + n_u_8_per_g_uo2_out * 1000 * sigma_f_238*10**-28)
	# print('macroscopic_fission_cross_section_out ' +str(macroscopic_fission_cross_section_out))
	int1 = integrate.quad(lambda r: r*special.jv(0, (r-rin) * 2.405 / R), rin, rin+R)[0]
	# print('int1 ' + str(int1))
	int2 = integrate.quad(lambda z: m.sin(z * m.pi / R), 0, R)[0]
	#  THIS IS PROBABLY WRONG IT IGNORES THE PART OF THE OUTER CORE THAT IS ON TOP OF THE INNER CORE
	# print('int2 ' + str(int2))
	power = 2 * m.pi * energy_per_fission * phi_0 * macroscopic_fission_cross_section_out * int1 * int2
	#something wrong with this function here????
	return power

# print('power outer ' + str(pow_out(19.5, 35*10*3.15*10**7, 35)))

for w in np.linspace(19.5, 99.5, 11):
	print('w = {0}, P = {1}'.format(w, pow_out(w, 35*15*3.15*10**7, 35)))

#why is this numbers imaginary or very large? code is almost identical to the envelope code so what the heck
print(rad_out(19.5))
print(rad_out(99.5))

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

fig, ax1 = plt.subplots()

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

# print('cooling rate tabulation')

for T in lifetime_options:
	row = [T]
	T *= 3.15 * 10 ** 7 #seconds
	for P in power_options: #years
		P *= 10 ** 6 #watts
		P /= efficiency #watts-thermal
		E = P*T
		row.append(mdot_in(E, P))
	tab.append(row)

for n in range(len(tab)):
    la = str(tab[n][0]) + 'y'
    plt.plot(power_options, tab[n][1:len(tab[n])+1], label=la)
plt.xlabel("Core Power [MWe]]")
plt.ylabel("Mass Flow Rate of Coolant [Kg/s]")
plt.legend()
plt.show()

# print(tabu(tab, tablefmt="github"))

tab = []

# print('inner outer power comparison')

for T in lifetime_options:
	row = [T]
	T *= 3.15 * 10 ** 7 #seconds
	for P in power_options: #years
		# P *= 10 ** 6 #watts
		P /= efficiency #watts-thermal
		E = P*T
		outer_power = pow_out(99.5, E = E, P = P)
		P *= 10**6
		ratio = outer_power / P
		outer_power /= (10**6)
		# row.append([round(outer_power, 6), round(ratio, 6)])
		row.append([outer_power, ratio])
	tab.append(row)

outer_power_axis = []
ratio_axis = []
inner_power_axis = []

power_options = np.linspace(85, 850, 11) #options in megawatts

for P in power_options: #years
	T = 10*3.15*10**7
	# P *= 10 ** 6 #watts
	P /= efficiency #watts-thermal
	E = P*T
	outer_power = pow_out(99.5, E = E, P = P)
	P *= 10**6
	ratio = outer_power / P
	outer_power /= (10**6)
	# row.append([round(outer_power, 6), round(ratio, 6)])
	inner_power_axis.append(P*10**-6)
	outer_power_axis.append(outer_power)
	ratio_axis.append(ratio)

fig, ax1 = plt.subplots()

ax1.plot(inner_power_axis, outer_power_axis, 'b-')
ax1.set_xlabel('Inner Core Power [MWth]')
# Make the y-axis label, ticks and tick labels match the line color.
ax1.set_ylabel('Outer Core Power [MWth]', color='b')
ax1.tick_params('y', colors='b')

ax2 = ax1.twinx()
# s2 = np.sin(2 * np.pi * t)
ax2.plot(inner_power_axis, ratio_axis, 'r.')
ax2.set_ylabel('Ratio', color='r')
ax2.tick_params('y', colors='r')

fig.tight_layout()
plt.show()

def k(wout, n, E): # n is number of pebbles, E is core energy in mj
	density_uo2 = 10970 #kg/m3
	wout /= 100
	r_in = vol_in(E*10**6)[1]
	r_out = rad_out(wout*100)
	# outer_vol = m.pi * (r_out**2 - r_in**2) * (r_out - r_in)
	outer_vol = m.pi * (r_out**2 - r_in**2) * r_in
	# bg = m.sqrt((2.405**2 / r_out**2) + (m.pi**2 / r_in ** 2))
	bg = m.sqrt((2.405**2 / r_out**2) + (m.pi**2 / r_out ** 2))
	r_pe = 0.5 * 10 ** -2 #  Pebble radius
	pe_vol = 4 / 3 * m.pi * r_pe ** 3
	# outer_vol = pe_vol * 100000
	max_n = m.pi / 3 / m.sqrt(2) * outer_vol / pe_vol #  max number of pebbles
	# mass_pe = density_uo2 * pe_vol
	density_uo2 *= n / max_n #  Updates the effective density of the fuel
	m_u_per_m_uo2_out = (238 * (1 - wout) + 235 * wout) / (238 * (1 - wout) + 235 * wout + 32)
	n_u5_per_g_uo2_out = wout * m_u_per_m_uo2_out * avogadro / 235
	n_u_8_per_g_uo2_out = (1 - wout) * m_u_per_m_uo2_out * avogadro / 238
	nu_sigf = density_uo2 * (n_u5_per_g_uo2_out * nu_sigma_f_235 * 1000 * 10**(-28) + n_u_8_per_g_uo2_out * 1000 * nu_sigma_f_238*10**-28)
	siga = density_uo2 * (n_u5_per_g_uo2_out * sigma_a_235 * 1000 * 10**(-28) + n_u_8_per_g_uo2_out * 1000 * sigma_a_238*10**-28)
	sigt = density_uo2 * (n_u5_per_g_uo2_out * sigma_t_235 * 1000 * 10**(-28) + n_u_8_per_g_uo2_out * 1000 * sigma_t_238*10**-28)
	if sigt < 10**-7:
		k = 0
	else:
		k = nu_sigf / (siga + 1/3/sigt * bg ** 2)
	return k

wout = 99.5
E = 35*10*3.15*10**7
r_in = vol_in(E*10**6)[1]
r_out = rad_out(wout)
# outer_vol = m.pi * (r_out**2 - r_in**2) * (r_out - r_in)
outer_vol = m.pi * (r_out**2 - r_in**2) * r_in
bg = m.sqrt((2.405**2 + m.pi**2) / (r_out ** 2))
r_pe = 0.5 * 10 ** -2 #  Pebble radius
pe_vol = 4 / 3 * m.pi * r_pe ** 3
max_n = m.pi / 3 / m.sqrt(2) * outer_vol / pe_vol #  max number of pebbles

k_axis = []
n_axis = np.linspace(0, max_n, 1000)
w_opt = np.linspace(19.5, 99.5, 11)
# n_axis = np.flip(n_axis, 0)

# power_options = np.linspace(35, 85, 11) #options in megawatts
# lifet = 10
# e_opts = power_options * lifet
# wout = 99.5
# # E = e_opts[-1]
# tab = []
#
# for E in e_opts:
# 	row = [E]
# 	E *= 3.15*10**7
# 	for n in n_axis: #years
# 		row.append(k(wout, n, E))
# 	tab.append(row)

# for n in range(len(tab)):
#     la = str(tab[n][0]) + 'MwY'
#     plt.plot(n_axis * -1, tab[n][1:len(tab[n])+1], label=la)
# plt.xlabel("Number of Pebbles")
# plt.ylabel("Keff")
# plt.legend()
# plt.show()

for n in n_axis:
	# print('k = ' + str(k(wout, n, E)))
	k_axis.append(k(wout, n, E))

n_axis *= -1
plt.plot(n_axis, k_axis)
plt.xlabel('number of pebbles')
plt.ylabel('keff')
plt.show()

def n(H, n0, t):
	n = n0 * (H - 9.81*t**2) / H
	if n != 0:
		n = n
	else:
		n = 0
	return n

t_axis = np.linspace(0, 0.3, 3601)
n_axis = []
k_axis = []

i = -1
for t in t_axis:
	i += 1
	n_axis.append(n(r_in, max_n, t))
	k_axis.append(k(wout, n_axis[i], E))

n_axis *= -1
plt.plot(t_axis, k_axis)
plt.xlabel('time [s]')
plt.ylabel('keff')
plt.show()
