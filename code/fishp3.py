import math as math
import numpy as np
import matplotlib.pyplot as plt
# from tabulate import tabulate as tabu
import scipy.integrate as integrate
import scipy.special as special

# constamts amd options

av = 6.022 * 10 ** 23
sigma_f5 = 2.156  # barns
nu_sigma_f5 = 5.297
nu_sigma_f8 = 0.142
sigma_f8 = 0.051
nu8 = 0.142 / sigma_f8
sigma_a8 = 0.404  # barns
sigma_a5 = 2.844  # barns
nu5 = 5.297 / sigma_f5
sigma_t5 = 8.246
sigma_t8 = 8.181
# density_uo2 = 10970  # kg/m3
density = 19100  # metal u density kg/m3
efficiency = 0.3
rod_radius = 0.5*10**-2
yr = 3.154 * 10 ** 7
epf = 200 * 1.6 * 10 ** -13

# operating conditions

pow = 85  # mwe
powin_opts = np.linspace(1, pow, pow * 10 + 1)
# powin = 65.5  # mwe
powin_choice = 3.05  # set the inner core power here in mwe
powin = powin_opts[int(powin_choice * 10)]
powou = pow - powin  # mwe

pow_t = pow / efficiency  # mwth
powin_t = powin / efficiency
powou_t = powou / efficiency

life = 10 * yr
nrg = pow_t * life  # mjth
nrgin = powin_t * life
nrgout = powou_t * life

# compute enrichment such that kinf = 1

beta = sigma_a8 / 238
lamda = nu5 * sigma_f5 / sigma_a5
alpha = sigma_a5 / 235 * (lamda + 1)
gamma = beta / alpha
enrichment = gamma / (1 + gamma)
# print('en = ' + str(enrichment))
enrichment_percent = enrichment*100
print(enrichment_percent)

# compute inner core dimensions


def vol_in(E):  # wants e in mjth
    E *= 10 ** 6
    fissions = E / epf
    Awgt = enrichment * 235 + (1 - enrichment) * 238
    N5PerG = Awgt * av
    fuel_mas = fissions / N5PerG / enrichment * 10 ** (-3)  # kg
    fuel_vol = fuel_mas / density  # m3
    fuel_rad = (fuel_vol / 2 / math.pi) ** (1/3)  # m
    return(fuel_vol, fuel_rad)


print(vol_in(nrgin))

# compute outer core dimensions


def rad_out(wout):  # wants wout as a percent
    wout /= 100
    N5PerG = wout * av / 235
    # Sigma_f_out = density_uo2 * N5PerG * 1000 * sigma_f5*10**-28
    N8PerG = (1 - wout) * av / 238
    nu_Sigma_f = density * (N5PerG * nu_sigma_f5 +
                            N8PerG * nu_sigma_f8) * 10 ** -25
    Sigma_a = density * (N5PerG * sigma_a5 + N8PerG * sigma_a8) * 10 ** -25
    Sigma_tr = density * (N5PerG * sigma_t5 + N8PerG * sigma_t8) * 10 ** -25
    diffusion_coefficient = 1 / 3 / Sigma_tr
    buckling = (nu_Sigma_f - Sigma_a) / diffusion_coefficient
    # print('w = ' + str(wout))
    # print('b= ' + str(buckling))
    outer_radius = ((2.405**2 + math.pi**2 / 4) / buckling) ** (1/2)
    # outer_vol = m.pi*(outer_radius**3)
    # return outer_vol
    return outer_radius


print('rad = {0}'.format(rad_out(19.5)))

# Given an inner core power, how much heat does outer core produce?


def pow_out(wout, P):  # wants percent and mwe
    E = P * life / efficiency
    P *= 10 ** 6
    R = rad_out(wout)
    wout /= 100
    rin = vol_in(E)[1]
    N5PerGin = enrichment * av / 235
    N8PerGin = (1 - enrichment) * av / 238
    Sigma_f_in = density * (N5PerGin * sigma_f5 +
                            N8PerGin * sigma_f8) * 10 ** -25
    phi_0 = P / vol_in(E)[0] / Sigma_f_in / epf
    n5PerGout = wout * av / 235
    n8PerGout = (1 - wout) * av / 238
    ro = R + rin
    Sigma_f_out = density * (n5PerGout * sigma_f5 +
                             n8PerGout * sigma_f8) * 10 ** -25
    int1 = integrate.quad(lambda r: special.jv(0, (r-rin) * 2.405 / R),
                          rin, ro)[0]
    # print('int1 ' + str(int1))
    int2 = integrate.quad(lambda z: math.sin(z * math.pi / 2 * ro),
                          0, 2 / ro)[0]
#  THIS IS TOO BIG IT IGNORES THE PART OF THE
#  OUTER CORE THAT IS ON TOP OF THE INNER CORE
    # print('int2 ' + str(int2))
    power = 2 * math.pi * epf * phi_0 * Sigma_f_out * int1 * int2
    power /= 10 ** 6
    power *= efficiency
    # something wrong with this function here????
    return power


# there is no viable inner core power
# print('unlimited powaaaah')
# for p in powin_opts:
#     q = pow_out(19.5, p) + p - pow
#     print(q)
#     if q <= 1:
#         print('yay')
#     else:
#         print('no good')


def phi_0(pow, wout):
    pow *= 10 ** 6 / efficiency
    R = rad_out(wout)
    wout /= 100
    n5PerGout = wout * av / 235
    n8PerGout = (1 - wout) * av / 238
    Sigma_f_out = density * (n5PerGout * sigma_f5 +
                             n8PerGout * sigma_f8) * 10 ** -25
    int1 = integrate.quad(lambda r: special.jv(0, r * 2.405 / R),
                          0, R)[0]
    # print('int1 ' + str(int1))
    int2 = integrate.quad(lambda z: math.sin(z * math.pi / (2 * R)),
                          0, 2 * R)[0]
    # THIS IS PROBABLY WRONG IT IGNORES THE PART
    # OF THE OUTER CORE THAT IS ON TOP OF THE
    # INNER CORE
    # print('int2 ' + str(int2))
    phi_0 = pow / (2 * math.pi * epf * Sigma_f_out * int1 * int2)
    return phi_0


print(phi_0(85, 19.5))


def mdot_out(P, wout):  # wants p in mwth, wout in percent
    R = rad_out(wout)
    # P *= 10 ** 6 / efficiency
    P /= efficiency
    H = 2 * R
    r = rod_radius #radius of fuel rod in m
    n5PerGout = wout * av / 235
    n8PerGout = (1 - wout) * av / 238
    Sigma_f_out = density * (n5PerGout * sigma_f5 +
                             n8PerGout * sigma_f8) * 10 ** -25
    ppeb = phi_0(P, wout) * epf * 4 / 3 * math.pi * r ** 3 * Sigma_f_out
    N = P / ppeb * 10 ** 6
    # N = R / r
    cp = 130  # j/kg/k
    delta_t = 500  # k
    # mdot = ppeb / cp / delta_t / N
    mdot = P * 10 ** 6 / cp / delta_t
    vdot = mdot / 9000
    return(N, mdot, vdot)


print(mdot_out(85, 19.5))
# print(mdot_out(85, 9.5))
