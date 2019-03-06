#tabulating many core parameters

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
efficiency = 0.3 
rod_radius = 0.5*10**-2

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

print(vol_in(85*10**6*10*3.15*10**7))

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

def pow_out(wout, E, P): #wants percent, joules, and watts
	wout /= 100
	P *= 10 ** 6
	E *= 10 ** 6
	R = rad_out(wout*100)
	rin = vol_in(E)[1]
	macroscopic_fission_cross_section_in = density_uo2 * (n_u5_per_g_uo2 * 1000 * sigma_f_235 * 10**-28 + n_u_8_per_g_uo2 * 1000 * sigma_f_238*10**-28)
	print('macroscopic_fission_cross_section_in ' + str(macroscopic_fission_cross_section_in))
	energy_per_fission = 200 * 1.6 * 10 ** -13 #joules
	phi_0 = P / vol_in(E)[0] / macroscopic_fission_cross_section_in / energy_per_fission
	print(phi_0)
	m_u_per_m_uo2_out = (238 * (1 - wout) + 235 * wout) / (238 * (1 - wout) + 235 * wout + 32)
	print('m_u_per_m_uo2_out ' + str(m_u_per_m_uo2_out))
	n_u5_per_g_uo2_out = wout * m_u_per_m_uo2_out * avogadro / 235
	print('n_u5_per_g_uo2_out ' + str(n_u5_per_g_uo2_out))
	n_u_8_per_g_uo2_out = (1 - wout) * m_u_per_m_uo2_out * avogadro / 238
	print('n_u_8_per_g_uo2_out ' + str(n_u_8_per_g_uo2_out))
	# Sigma_f_out = density_uo2 * n_u5_per_g_uo2_out * 1000 * sigma_f_235*10**-28
	macroscopic_fission_cross_section_out = density_uo2 * (n_u5_per_g_uo2_out * sigma_f_235 * 1000 * 10**(-28) + n_u_8_per_g_uo2_out * 1000 * sigma_f_238*10**-28)
	print('macroscopic_fission_cross_section_out ' +str(macroscopic_fission_cross_section_out))
	int1 = integrate.quad(lambda r: r*special.jv(0, r * 2.405 / R), rin, R)[0]
	print('int1 ' + str(int1))
	int2 = integrate.quad(lambda z: m.sin(z * m.pi / R), 0, R)[0]
	print('int2 ' + str(int2))
	power = 2 * m.pi * energy_per_fission * phi_0 * macroscopic_fission_cross_section_out * integrate.quad(lambda r: r*special.jv(0, r * 2.405 / R), rin, R)[0] * integrate.quad(lambda z: m.sin(z * m.pi / R), 0, R)[0]
	#something wrong with this function here????
	return power

print('power outer ' + str(pow_out(99.5, 35*3.15*10**7, 35)))


#why is this numbers imaginary or very large? code is almost identical to the envelope code so what the heck
print(rad_out(19.5))
print(rad_out(99.5))



# for i in np.linspace(0, 100, 101):
# 	print(vol_out(i))

# print(vol_out(19.5))

# def power_outer(E, P, wout): #energy, power, outer core enrichment
# 	inner_volume = vol_in(E)[0]
# 	power_density = P / inner_volume
# 	Simga_f_in = density_uo2 * n_u5_per_g_uo2 * 1000 * sigma_f_235*10**-28
# 	phi_0 = power_options / Simga_f_in
# 	m_u_per_m_uo2_out = (238 * (1 - wout) + 235 * wout) / (238 * (1 - wout) + 235 * wout + 32)
# 	n_u5_per_g_uo2_out = wout * m_u_per_m_uo2 * avogadro / 235
# 	Sigma_f_out = density_uo2 * n_u5_per_g_uo2_out * 1000 * sigma_f_235*10**-28
# 	n_u_8_per_g_uo2 = (1 - enrichment) * m_u_per_m_uo2_out * avogadro / 238
# 	macroscopic_fission_cross_section = density_uo2 * n_u5_per_g_uo2 * 1000 * sigma_f_235*10**-28
# 	macroscopic_absorption_cross_section = density_uo2 * (n_u5_per_g_uo2 * 1000 * sigma_a_235*10**-28 + n_u_8_per_g_uo2 * 1000 * sigma_a_238*10**-28)
# 	macroscopic_transport_cross_section =  density_uo2 * (n_u5_per_g_uo2 * 1000 * sigma_t_235*10**-28 + n_u_8_per_g_uo2 * 1000 * sigma_t_238*10**-28)
# 	diffusion_coefficient = 1 / 3 / macroscopic_transport_cross_section
# 	buckling = (nu * macroscopic_fission_cross_section - macroscopic_absorption_cross_section) / diffusion_coefficient 
# 	outer_radius = ( (2.405**2 + m.pi**2) / buckling ) ** (1/2) 
# 	outer_power = Sigma_f_out * phi_0 * integrate.quad(lambda r: special.jv(0,r*2.405/ outer_radius), 0, outer_radius) * integrate.quad(lambda h: math.cos(h * math.pi / outer_radius), 0, outer_radius)
# 	return(outer_power)

# print(power_outer(35*10**6*3.15*10**7, 35*10**6, 19.5))

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

print('cooling rate tabulation')

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

tab = []

# print('inner outer power comparison')

# for T in lifetime_options:
# 	row = [T]
# 	T *= 3.15 * 10 ** 7 #seconds
# 	for P in power_options: #years
# 		# P *= 10 ** 6 #watts
# 		P /= efficiency #watts-thermal
# 		E = P*T
# 		outer_power = pow_out(19.5, E = E, P = P)
# 		P *= 10**6
# 		ratio = outer_power / P
# 		outer_power /= (10**6)
# 		# row.append([round(outer_power, 6), round(ratio, 6)])
# 		row.append([outer_power, ratio])
# 	tab.append(row)

# print(tabu(tab, tablefmt="github"))