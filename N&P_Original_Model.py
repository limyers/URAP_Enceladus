import math
import sys
from astropy import units as u
from astropy import constants as const

# Physical constants

mmw = 18.0e-3 * (u.kg / u.mol) # mean molecular weight of water (kg/mol)
R = 8.314 * (u.J / (u.mol * u.K)) # universal gas constant (J/mol * K)
P_0 = 1.4e11 * u.Pa # initial pressure (Pa)
L_LtV = 43655.0 * (u.J / u.kg) # latent heat from liquid to vapor (J/kg)
L_VtI = 2.80e6 * (u.J / u.kg) # latent heat from vapor to ice (J/kg)
g_E = 0.11 * (u.m / (u.s**2)) # Enceladus' gravity (m/s^2)
gamma = 1.333333 # heat capacity ratio of gas (C_p/C_v)
C_D = 0.002 # drag coefficient

sigma = const.sigma_sb
K = 2.4 * (u.W / (u.m * u.K))# thermal conductivity of surrounding ice (W/m * K)
C_p = 2000.0 * (u.J / u.kg / u.K) # specific heat of gas (J/kg/K)
C_v = 2000.0 * (u.J / u.kg / u.K) # specific heat of ice (J/kg/K)

A = 3.63e12 * u.Pa # constant for p_w (saturation vapor pressure at T_wall) equation (Pa)
B = 6147.0 * u.K # constant for p_w (saturation vapor pressure at T_wall) equation (K)
T_E = 68.0 * u.K # effective temperature at surface
T_bottom = 273 * u.K # temperature of subsurface ocean at bottom (K)

h = 1500 * u.m # height of the wall (m)
delta = 0.05 * u.m # crack width
s = 0.0 # solid ice mass fraction
ratio = 0.961231 # pressure at the bottom = ratio * saturation vapor pressure of liquid
dz = 0.001 # size of steps in z-direction
ddelta_dz = 0 #d/dz(delta)
output_number = 10000 # interval of outputs

def main():
    # Need to strip values of units so that it doesn't increase runtime exponentially
    val_P_0 = P_0.to_value(u.Pa)
    val_L_LtV = L_LtV.to_value(u.J / u.kg)
    val_R = R.to_value(u.J / (u.mol * u.K))
    val_mmw = mmw.to_value(u.kg / u.mol)
    val_h = h.to_value(u.m)
    val_dz = dz.to_value(u.m)
    val_A = A.to_value(u.Pa)
    val_B = B.to_value(u.K)
    val_K = K.to_value(u.W / (u.m * u.K))
    val_T_E = T_E.to_value(u.K)
    val_sigma = sigma_sb.value
    val_L_VtI = L_VtI.to_value(u.J / u.kg)

    # Initial conditions (at z = 0)
    T_flow = T_bottom.to_value(u.K) # temperature at surface (assigned to T_flow as to not rewrite T_E)
    P_l_bottom = val_P_0 * math.exp(-val_L_LtV / val_R / T_flow)
    P_g_bottom = P_l_bottom * ratio

    rho = P_g_bottom * val_mmw / val_R / T_flow/ (1.0 - s) # density
    u_energy = P_g_bottom / rho / (gamma - 1.0) # internal energy
    C = (P_l_bottom - P_g_bottom) / math.sqrt(2.0 * math.pi * val_R / val_mmw * T_flow) # rho * velocity
    velocity = C / rho # velocity

    # Changing variables
    E_mass_flux = 0.0 # mass flux
    w = 0 # step counter for while loop
    z = 0.0 # distance from bottom
    d = 0.0 # depth from surface
    dEsum = 0 # cumulative mass flux
    dFsum = 0 # cumulative heat flux
    P_current = P_g_bottom # tracks changing pressure as z increases

    #

    while P_current > 0.0 and z < (val_h + 200.0): # upper limit is wall height + 200 m
        w += 1

        if w % output_number == 0: # Logs output every 10,000 steps
            if w == output_number:
                print("# z, d, velocity, T_flow, rho, P_current")

            print(f"{z:13.5e} {d:13.5e} {velocity:13.5e} {T_flow:13.5e} {rho:13.5e} {P_current:13.5e}")

        P_here = P_current
        d = abs(val_h - z)
        P_new = 0.5 * P_current # arbitrary value so that the while loop is triggered
        flag = 0

        while abs((P_new - P_here) / P_new) > 0.001:
            if flag == 0:
                P_new = P_here # undoes fake guess of P_new = 1/2 * P_current
            P_here = P_new
            flag = 1

            E_mass_flux = 0.0

            FFF, T0, T_wall, E_mass_flux = find_F(P_here, T_surface, d, val_A, val_B, val_K, val_T_E, val_sigma, val_R, val_mmw, val_L_VtI)

            FFF *= 2 # factor of two because heat flux goes into both walls

            if FFF > 0.0:
                E_mass_flux = FFF / val_L_VtI
            else:
                print('F is negative')
                sys.exit(1)

        # Adding to cumulative fluxes
        dFsum += FFF * math.sqrt(1.0 + ddelta_dz**2.0) * val_dz
        dEsum += E_mass_flux * math.sqrt(1.0 + ddelta_dz**2.0) * val_dz

        # insert solving differential equations
        # P_new = P_current + dpdz * val_dz

        if P_new < 0.0:
            print('Pressure is negative. ...')
            sys.exit(1)


