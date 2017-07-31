import pandas as pd
import re
from chemhelper import elements

def mass(chemical):
    '''

    :param chemical: compound chemical symbol (example: 'H2O')
    :return: molar mass of compound
    '''
    chemical = re.findall(r'([A-Z][a-z]*)(\d*)', chemical)
    m = 0
    for tup in chemical:
        if tup[1] == '':
            m+=elements.elements[tup[0]].mass
        else:
            m+= elements.elements[tup[0]].mass * int(tup[1])
    return round(m,2)

def mass_steps(chemical):
    """

    :param chemical: compound chemical symbol (example: 'H2O')
    :prints: Step by step calculation of compound's mass
    """
    chem = re.findall(r'([A-Z][a-z]*)(\d*)', chemical)
    m = 0
    numbers = []
    print(chemical)
    for tup in chem:
        if tup[1] == '':
            print(tup[0] + ': ' + str(elements.elements[tup[0]].mass))
            m+=elements. elements[tup[0]].mass
            numbers.append(elements.elements[tup[0]].mass)

        else:
            print(tup[0] + ': ' + tup[1] + 'x' + str(elements.elements[tup[0]].mass) + ' = ' + str(elements.
                elements[tup[0]].mass * int(tup[1])))
            m+= elements.elements[tup[0]].mass * int(tup[1])
            numbers.append(elements.elements[tup[0]].mass * int(tup[1]))
    n = len(numbers)
    i = 1

    string = str(numbers[0])
    while i < n:
        string = string + ' + ' +  str(numbers[i ])
        i+=1
    print( '\nMass ' + chemical + ' = ' + string + ' = ' + str(round(m,2)))


def composition(chemical):
    '''
    :param chemical: chemical compound symbol (example 'H2O')
    :return: pandas Series containing mass due to each element in compound and total mass of compound
    '''
    total_mass = mass(chemical)
    chem = re.findall(r'([A-Z][a-z]*)(\d*)', chemical)
    masses = []
    ind = []
    for tup in chem:
        if tup[1] == '':
            masses.append(elements.elements[tup[0]].mass)
            ind.append(tup[0])

        else:
            masses.append(elements.elements[tup[0]].mass * int(tup[1]))
            ind.append(tup[0])
    ind.append('Total:')
    masses.append(total_mass)
    return pd.Series(masses,index=ind,name= chemical)


def composition_steps(chemical):
    '''
    :param chemical: chemical compound symbol (example 'H2O')
    :prints: step by step solution of mass due to each element in compound and total mass
    '''
    mass_steps(chemical)
    print('')
    print(composition(chemical))


def percent_composition(chemical):
    '''
    :param chemical: chemical compound symbol (example 'H2O')
    :return: pandas Series containing percent composition of each element making up compound.
    '''
    return round((composition(chemical)/mass(chemical))* 100,3)


def percent_composition_steps(chemical):
    """
    :param chemical: chemical compound symbol (example 'H2O')
    :prints: step by step solution of solving the percent composition of each element making up the compound.
    """
    composition_steps(chemical)
    comp = composition(chemical)
    component_series = comp.drop('Total:')
    i=0
    for component in component_series:
        el = component_series.index[i]
        el_mass = component_series[i]
        chem_mass = mass(chemical)

        comp_str = "{}: 100 x {:.2f}g {}/{:.2f}g {} = {:.2f}%".format(el,el_mass,el,chem_mass,chemical,percent_composition(chemical)[i])
        print(comp_str)
        i +=1


