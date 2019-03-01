# deriving the enrichment, mass of fuel, and coolant flow rate for earliest 
# some usefu constant from JAEA, fission-energy-spectrum-averaged cross section
import math as m
import numpy as np
avogadro = 6.022 * 10 ** 23
sigma_f_235 = 2.156 # barns
sigma_a_238 = 0.404 # barns
sigma_a_235 = 2.844 # barns
nu = 5.297 / sigma_f_235
sigma_t_235 = 8.246
sigma_t_238 = 8.181
density_uo2 = 10970 #kg/m3

# compute enrichment such that kinf = 1

beta = sigma_a_238 / 238
lamda = nu * sigma_f_235 / sigma_a_235
alpha = sigma_a_235 / 235 * (lamda + 1)
gamma = beta / alpha
enrichment = gamma / (1 + gamma)
enrichment_percent = enrichment*100
print('w = %s percent' %enrichment_percent)

# compute Uranium counts per mass of fuel

mass_uranium_per_mass_uo2 = (238 * (1 - enrichment) + 235 * enrichment) / (238 * (1 - enrichment) + 235 * enrichment + 32)
n_uranium_5_per_gram_uo2 = enrichment * mass_uranium_per_mass_uo2 * avogadro / 235
print('mu/mu02 = %s' %mass_uranium_per_mass_uo2)
print('nu5/mu02 = %s' %n_uranium_5_per_gram_uo2) #there will actually be more fissile atoms than this once Pu build up

# pick reactor power

power_options = [60*10**6, 60*10**7, 60*10**8, 60*10**5, 60*10**4, 60*10**3, 30*10**6, 30*10**7, 30*10**8, 30*10**5, 30*10**4, 30*10**3]
Power_electric = power_options[0]
thermal_efficiency = 0.3
power_thermal = Power_electric / thermal_efficiency
print('pth = %s' %power_thermal)

# uranium counting

fission_rate = power_thermal / (1.6*10**-13) / 200
fissions_per_year = fission_rate * 3.15*10**7
years_per_cycle = 2
fissions_per_cycle = years_per_cycle * fissions_per_year
print('fissions per cycle = %s' %fissions_per_cycle)

# core dimensions

mass_per_cycle = fissions_per_cycle / n_uranium_5_per_gram_uo2 / 1000
volume_per_cycle = mass_per_cycle / density_uo2
radius_fuel = (3 / 4 / m.pi * volume_per_cycle) ** (1/3)
height = radius_fuel * 4/3
check_vol = m.pi*radius_fuel**2*height
print('mass = %s kg' %mass_per_cycle)
print('volume = %s m3' %volume_per_cycle)
print('volume check = %s' %check_vol)
print('inner height = %s m' %height)
print('inner radius = %s' %radius_fuel)

# densities

power_density_mas = (mass_per_cycle / power_thermal)**-1
power_density_vol = (volume_per_cycle / power_thermal)**-1
power_density_hgt = power_thermal / height
print('p/m = %s' %power_density_mas)
energy_density_mass = power_thermal *6.3*10**7 / mass_per_cycle
print('e/m = %s' %energy_density_mass)

# thermal calculation

rod_radius = 10**-2
number_rods = m.floor(radius_fuel * 100)
heated_perimeter = 2 * m.pi * number_rods * rod_radius

cp = 130 # j/kg/k 
delta_t = 1000 # k
m_dot = power_thermal / cp / delta_t
v_dot = m_dot / 9000

print('m dot of coolant = %s' %m_dot)
print('v dot of coolant = %s' %v_dot)

# some outer core calculations
macroscopic_fission_cross_section = density_uo2 * n_uranium_5_per_gram_uo2 * 1000 * sigma_f_235*10**-28
phi = fission_rate / macroscopic_fission_cross_section
print('phi = %s' %phi)

total_neutrons = phi*m.pi*radius_fuel*height*2*10**4
t_half_pu = 24000*3.15*10**7
lamda_pu = m.log(2)/t_half_pu
number_pu_to_supply_phi = total_neutrons / lamda_pu
print(number_pu_to_supply_phi)

t_half_cf_250 = 13.08*3.15*10**7
lamda_cf_250 = m.log(2)/t_half_cf_250
number_cf_250_to_supply_phi = total_neutrons / lamda_cf_250
print(number_cf_250_to_supply_phi)

#radius of outer core as a function of w

# enrichment_options = np.linspace(19.5, 99.5, 161)
# enrichment = enrichment_options[0]
enrichment = 19.5/100
mass_uranium_per_mass_uo2 = (238 * (1 - enrichment) + 235 * enrichment) / (238 * (1 - enrichment) + 235 * enrichment + 32)
n_uranium_5_per_gram_uo2 = enrichment * mass_uranium_per_mass_uo2 * avogadro / 235
n_uranium_8_per_gram_uo2 = (1 - enrichment) * mass_uranium_per_mass_uo2 * avogadro / 238
macroscopic_fission_cross_section = density_uo2 * n_uranium_5_per_gram_uo2 * 1000 * sigma_f_235*10**-28
macroscopic_absorption_cross_section = density_uo2 * (n_uranium_5_per_gram_uo2 * 1000 * sigma_a_235*10**-28 + n_uranium_8_per_gram_uo2 * 1000 * sigma_a_238*10**-28)
macroscopic_transport_cross_section =  density_uo2 * (n_uranium_5_per_gram_uo2 * 1000 * sigma_t_235*10**-28 + n_uranium_8_per_gram_uo2 * 1000 * sigma_t_238*10**-28)
diffusion_coefficient = 1 / 3 / macroscopic_transport_cross_section
buckling = (nu * macroscopic_fission_cross_section - macroscopic_absorption_cross_section) / diffusion_coefficient 
print(buckling)
outer_radius = ( (2.405**2 + 9/16 * m.pi**2) / buckling ) ** (1/2) 
print(outer_radius)


def radius_finder(enrichment):
	mass_uranium_per_mass_uo2 = (238 * (1 - enrichment) + 235 * enrichment) / (238 * (1 - enrichment) + 235 * enrichment + 32)
	n_uranium_5_per_gram_uo2 = enrichment * mass_uranium_per_mass_uo2 * avogadro / 235
	n_uranium_8_per_gram_uo2 = (1 - enrichment) * mass_uranium_per_mass_uo2 * avogadro / 238
	macroscopic_fission_cross_section = density_uo2 * n_uranium_5_per_gram_uo2 * 1000 * sigma_f_235*10**-28
	macroscopic_absorption_cross_section = density_uo2 * (n_uranium_5_per_gram_uo2 * 1000 * sigma_a_235*10**-28 + n_uranium_8_per_gram_uo2 * 1000 * sigma_a_238*10**-28)
	macroscopic_transport_cross_section =  density_uo2 * (n_uranium_5_per_gram_uo2 * 1000 * sigma_t_235*10**-28 + n_uranium_8_per_gram_uo2 * 1000 * sigma_t_238*10**-28)
	diffusion_coefficient = 1 / 3 / macroscopic_transport_cross_section
	buckling = (nu * macroscopic_fission_cross_section - macroscopic_absorption_cross_section) / diffusion_coefficient ** (1/2)
	outer_radius = (2.405**2 + 9/16 * m.pi**2) / buckling ** (1/2)
	return(outer_radius)

# for val in enrichment_options:
# 	print("w = " + str(val) + " r = " + str(radius_finder(val)))