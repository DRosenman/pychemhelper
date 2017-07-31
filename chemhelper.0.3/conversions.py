# -*- coding: utf-8 -*-
"""
Created on Thu Jul 27 17:34:21 2017

@author: rosen
"""

from chemhelper.elements import *
from chemhelper import compounds


def grams_to_moles(mass, chemical):
    '''
    :param mass: mass of sample of substance (in grams)
    :param chemical:  chemical formula of substance
    :returns: moles of the sample of the substance
    '''
    molar_mass = compounds.mass(chemical)
    moles = round(mass / molar_mass, 2)
    return moles


def grams_to_moles_steps(mass, chemical):
    '''
    :param mass: mass of sample of substance (in grams)
    :param chemical:  chemical formula of substance
    :prints: solution of conversion from grams to moles of the sample of the substance
    '''
    molar_mass = compounds.mass(chemical)
    moles = grams_to_moles(mass, chemical)
    print(str(mass) + 'g ' + chemical + ' x ' + ' 1 mole ' + chemical + '/' + str(
        molar_mass) + 'g ' + chemical + ' = ' + str(moles) + ' moles ' + chemical)


def moles_to_grams(moles, chemical):
    return moles * compounds.mass(chemical)


def moles_to_grams_steps(moles, chemical):
    grams = moles_to_grams(moles, chemical)
    compounds.mass_steps(chemical)
    mass = compounds.mass(chemical)
    print("{:.2f} moles x {:.2f} grams/mole = {:.2f} grams".format(moles, mass, grams))


def molarity(volume, mass=None, chemical=None, moles=None):
    '''
    :param volume: volume of solution (in L)
    :param mass: mass of the solute (in g), default = None
    :param chemical: the solute, default = None
    :param moles: # of moles of solute, default= None
    :return: molarity of solution (moles of substance/L of solution)
    '''
    if mass is None:
        return moles / volume
    else:
        moles = grams_to_moles(mass, chemical)
        return moles / volume


def molarity_steps(volume, mass=None, chemical=None, moles=None):
    '''
    :param volume: volume of solution (in L)
    :param mass: mass of the solute (in g), default = None
    :param chemical: the solute, default = None
    :param moles: # of moles of solute, default= None
    :prints: step by step determination of the molarity of the solution
    '''
    if moles is None:
        grams_to_moles_steps(mass, chemical)
        moles = grams_to_moles(mass, chemical)
        M = molarity(volume, mass, chemical)
    else:
        M = molarity(volume, moles=moles)
    print("{:.2f} moles/{:.2f} L = {:.2f} M".format(moles, volume, M))


def molarity_to_moles(M, volume):
    """
    :param M: the concentration of the solution in moles per L
    :param volume: the volume of the solution in Liters
    returns: moles of solute
    """
    return M * volume


def molarity_to_moles_steps(M, volume):
    print("{:.2f} moles/L x {:.2f} L = {:.2f} moles".format(M, volume, molarity_to_moles(M, volume)))


def molarity_to_grams(molarity, volume, chemical):
    '''
    :param molarity: molarity of solution (moles of solute/ L of solution)
    :param volume: volume of solution (in L)
    :param chemical: chemical formula of solute
    :return: mass of solute (in grams)
    '''
    molar_mass = compounds.mass(chemical)
    moles = molarity * volume
    return round(molar_mass * moles, 2)


def molarity_to_grams_steps(molarity, volume, chemical):
    moles = molarity_to_moles(molarity, volume)
    moles_to_grams_steps(moles, chemical)


