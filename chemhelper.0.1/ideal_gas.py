# -*- coding: utf-8 -*-
"""
Created on Thu Jul 27 18:41:11 2017

@author: rosen
"""

from chemsolver import compounds
from chemsolver import conversions
from chemsolver import elements


def ideal_gas(solve_for, substance=None, P=None, V=None, n=None, T=None, m=None, decimals=2, Tunits='K', Vunits='L',
              Punits='atm'):
    R = 0.0821

    if (solve_for != 'P') and (Punits == 'torr'):
        P = round(P / 760, 2)
    if (solve_for != 'T') and (Tunits == 'C'):
        T = T + 273.5
    if (solve_for != 'V') & (Vunits != 'L'):
        V = round(V / 1000, 2)
    if m != None:
        n = round(m / compounds.mass(substance), 2)
    if solve_for == 'P':
        P = n * R * T / V
        return round(P, 2)
    if solve_for == 'n':
        n = P * V / (R * T)
        return round(n, 2)
    if solve_for == 'V':
        V = n * R * T / P
        return round(V, 2)
    if solve_for == 'T':
        T = P * V / (n * R)
        return (T, 2)


def ideal_gas_steps(solve_for, substance=None, P=None, V=None, n=None, T=None, m=None, decimals=2, Tunits='K',
                    Vunits='L', Punits='atm'):
    """
    :params
            solve_for: the value you want to solve for:
                        'P' = pressure, 'V' = volume in liters', 'n' = moles', 'T' = temperature in 'K'
            substance: the chemical
            P: pressure default None, units atm default
            V: volume, default None, units L
            n: moles, default None
            T: temperature, default None, units Kelvin
    """

    R = 0.0821
    ideal_gas_equation = 'PV = nRT'
    if (solve_for != 'P') and (Punits == 'torr'):
        print('P = {} torr * 1 atm/760 torr = {} atm'.format(P, round(P / 760), 2))
        P = round(P / 760, 2)
    if (solve_for != 'T') and (Tunits == 'C'):
        print('T = ({} + 273.5)K = {}'.format(T, T + 273.5))
        T = T + 273.5
    if (solve_for != 'V') & (Vunits != 'L'):
        print('V = {} mL * 1L/1000mL = {} L'.format(V, round(V / 1000, 2)))
        V = round(V / 1000, 2)
    if m != None:
        n = round(m / compounds.mass(substance), 2)
        print('n = {}g * 1 mole {}/{} g {} = {} moles {}'.format(m, substance, compounds.mass(substance), substance, n,
                                                                 substance))

    if solve_for == 'P':
        P = n * R * T / V
        print(ideal_gas_equation)
        print('P = nRT/V')
        print('P = {} moles * (0.0821 L*atm/mol*K) * {}K/{} L = {} atm'.format(n, T, round(V, 2), round(P, 2)))
    if solve_for == 'n':
        n = P * V / (R * T)
        print(ideal_gas_equation)
        print('n = PV/RT')
        print('n = ({} atm * {} L)/((0.0821 L*atm/mol*K) * {} K) = {} moles'.format(P, V, T, n))
    if solve_for == 'V':
        V = n * R * T / P
        print(ideal_gas_equation)
        print('V = nRT/P')
        print('V = {} moles * (0.0821 L*atm/mol*K) * {} K / {} atm = {} L'.format(n, T, P, round(V,2)))
    if solve_for == 'T':
        T = P * V / (n * R)
        print('T = PV/nR')
        print(ideal_gas_equation)
        print('T = {} atm * {} L/({} moles * 0.0821 L*atm/mole*K = {} K)'.format(P, V, n, T))
