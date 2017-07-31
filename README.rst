
ChemSolver Package: Version 0.4
===============================

Purpose
-------

I created this package to aid me in tutoring chemistry students for
tutor.com. The package has many functions with the capability to solve
simple chemistry problems. Each of these functions of a corresponding
function with the same name with \*\*"\_steps"\*\* tacked on a at the
end of the function's name that provided a step by step solution to the
problem or step. The goal of this package is to enable me to quickly
check student's work and also display solutions, while being able to
focus on checking the student's understanding and catering to his or her
needs.

Installation
------------

.. code:: ruby

    pip install chemsolver

or available for download at
https://pypi.python.org/pypi?name=chemsolver&version=0.1&:action=display

Modules
-------

-  **chemsolver.elements**

   -  Importing the module

   -  Element Class

      -  Element Objects

         -  Element Object Attributes
         -  Examples

   -  element\_index
   -  element\_objects
   -  elements

-  **chemsolver.compounds**

   -  Importing the module
   -  Functions

      -  mass
      -  mass\_steps

         -  mass and mass\_steps examples

      -  composition
      -  composition\_steps

         -  composition and composition\_steps examples

      -  percent\_composition
      -  percent\_composition\_steps

         -  percent\_composition and percent\_composition\_steps
            examples

-  **chemsolver.conversions**

   -  Importing the module
   -  Functions

      -  grams\_to\_moles
      -  grams\_to\_moles\_steps

         -  grams\_to\_moles and grams\_to\_moles\_steps examples

      -  moles\_to\_grams
      -  moles\_to\_grams\_steps

         -  moles\_to grams and moles\_to\_grams\_steps examples

      -  molarity
      -  molarity\_steps

         -  molarity and molarity\_steps examples

      -  molarity\_to\_moles
      -  molarity\_to\_moles\_steps

         -  molarity\_to\_moles and molarity\_to\_moles\_steps examples

      -  molarity\_to\_grams
      -  molarity\_to\_grams\_steps

         -  molarity and molarity\_steps examples

-  **chemsolver.ideal\_gas**

   -  Importing the module
   -  Functions

      -  ideal\_gas
      -  ideal\_gas\_steps

         -  ideal\_gas and ideal\_gas\_steps examples

chemsolver.elements 
~~~~~~~~~~~~~~~~~~~~

-  **Importing the module**\ 

.. code:: ipython3

    from chemsolver import elements

-  **Element Class**
-  The elements module contains an **Element object** for each of the 93
   naturally occuring elements.

   .. code:: python

                class Element:
                    def __init__(self,symbol,name,atomic_number,mass):
                   '''

                   :param symbol: chemical symbol of element
                   :param name:  element name
                   :param atomic_number:  atomic number of element
                   :param mass: element mass
                   '''
                   self.symbol = symbol
                   self.name = name
                   self.mass = mass
                   self.atomic_number = atomic_number

                   def __str__(self):
                       return 'Element: ' + self.name + ', ' + self.symbol  +'\nAtomic Number: ' +str(self.atomic_number) + '\nMass: {}'.format(self.mass)

-  **Element Objects:** The variable names of the elements are the
   chemical symbols of the elements (example: hydrogen is H)

   -  **Attributes**:

      -  *self.symbol*: The chemical symbol of the element (example
      -  *self.name*: The name of the element
      -  *self.mass*: The molar mass of the element (also numerically
         equal to the mass of 1 atom of the element in atm)
      -  *self.atomic\_number*: The atomic number of the element

***Examples***:

.. code:: ipython3

    print(elements.H)


.. parsed-literal::

    Element: Hydrogen, H
    Atomic Number: 1
    Mass: 1.008
    

.. code:: ipython3

    elements.Ti.name




.. parsed-literal::

    'Titanium'



.. code:: ipython3

    elements.Ag.mass




.. parsed-literal::

    107.87



.. code:: ipython3

    elements.K.atomic_number




.. parsed-literal::

    19



.. code:: ipython3

    print(elements.Ca)


.. parsed-literal::

    Element: Calcium, Ca
    Atomic Number: 20
    Mass: 40.08
    

-  **element\_index:** List of element symbols.

``python       element_index = ["Ni","Cu","Zn","Ga","Ge","As","Se","Br","Kr","Rb","Sr","Y","Zr","Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd","In","Sn","Sb","Te","I","Xe","Cs","Ba","La","Ce","Pr","Nd","Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb","Lu","Hf","Ta","W","Re","Os","Ir","Pt","Au","Hg","Tl","Pb","Bi","Po","At","Rn","Fr","Ra","Ac","Th","Pa","U"]``

-  **element\_objects:** List of element objects

   .. code:: python

           element_objects = [H,He,Li,Be,B,C,N,O,F,Ne,Na,Mg,Al,Si,P,S,Cl,Ar,K,Ca,Sc,Ti,V,Cr,Mn,Fe,Co,Ni,Cu,Zn,Ga,Ge,As,Se,Br,Kr,Rb,Sr,Y,Zr,Nb,Mo,Tc,Ru,Rh,Pd,Ag,Cd,In,Sn,Sb,Te,I,Xe,Cs,Ba,La,Ce,Pr,Nd,Pm,Sm,Eu,Gd,Tb,Dy,Ho,Er,Tm,Yb,Lu,Hf,Ta,W,Re,Os,Ir,Pt,Au,Hg,Tl,Pb,Bi,Po,At,Rn,Fr,Ra,Ac,Th,Pa,U]

-  **elements**: a pandas series containing every element object.

chemsolver.compounds 
~~~~~~~~~~~~~~~~~~~~~

-  **Importing the module**

.. code:: ipython3

    from chemsolver import compounds

-  **Functions**

   ***mass***: determines the molar mass of a molecule and/or compound

\`\`\`python def mass(chemical): '''

::

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
    return round(m,2)```

***mass\_steps***: shows the calculation of the the molar mass of a
molecule and/or compound

\`\`\`python def mass\_steps(chemical): """

::

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
    ```

*Examples*:

.. code:: ipython3

    compounds.mass('HNO4')




.. parsed-literal::

    79.02



.. code:: ipython3

    compounds.mass_steps('HNO4')


.. parsed-literal::

    HNO4
    H: 1.008
    N: 14.01
    O: 4x16.0 = 64.0
    
    Mass HNO4 = 1.008 + 14.01 + 64.0 = 79.02
    

****composition****: returns the composition of compounds by mass

\`\`\`python def composition(chemical): ''' :param chemical: chemical
compound symbol (example 'H2O') :return: pandas Series containing mass
due to each element in compound and total mass of compound '''
total\_mass = mass(chemical) chem =
re.findall(r'([A-Z][a-z]\*)(:raw-latex:`\d*`)', chemical) masses = []
ind = [] for tup in chem: if tup[1] == '':
masses.append(elements.elements[tup[0]].mass) ind.append(tup[0])

::

            else:
                masses.append(elements.elements[tup[0]].mass * int(tup[1]))
                ind.append(tup[0])
        ind.append('Total:')
        masses.append(total_mass)
        return pd.Series(masses,index=ind,name= chemical) 
    ```

****composition\_steps****: step by step solution of the chemical
composition of a compound by mass

.. code:: python

    def composition_steps(chemical):
        '''
        :param chemical: chemical compound symbol (example 'H2O')
        :prints: step by step solution of mass due to each element in compound and total mass
        '''
        mass_steps(chemical)
        print('')
        print(composition(chemical))

*Examples:*

.. code:: ipython3

    compounds.composition('NH4')




.. parsed-literal::

    N         14.010
    H          4.032
    Total:    18.040
    Name: NH4, dtype: float64



.. code:: ipython3

    compounds.composition_steps('NH4')


.. parsed-literal::

    NH4
    N: 14.01
    H: 4x1.008 = 4.032
    
    Mass NH4 = 14.01 + 4.032 = 18.04
    
    N         14.010
    H          4.032
    Total:    18.040
    Name: NH4, dtype: float64
    

****percent\_composition****: calculates the percent composition of each
element in a compound.

.. code:: python

    def percent_composition(chemical):
        '''
        :param chemical: chemical compound symbol (example 'H2O')
        :return: pandas Series containing percent composition of each element making up compound.
        '''
        return round((composition(chemical)/mass(chemical))* 100,3)

****percent\_composition\_steps****: calculates and shoes, step by step,
the percent composition of each element in a compound \`\`\`python def
percent\_composition\_steps(chemical): """ :param chemical: chemical
compound symbol (example 'H2O') :prints: step by step solution of
solving the percent composition of each element making up the compound.
""" composition\_steps(chemical) comp = composition(chemical)
component\_series = comp.drop('Total:') i=0 for component in
component\_series: el = component\_series.index[i] el\_mass =
component\_series[i] chem\_mass = mass(chemical)

::

        comp_str = "{}: 100 x {:.2f}g {}/{:.2f}g {} = {:.2f}%".format(el,el_mass,el,chem_mass,chemical,percent_composition(chemical)[i])
        print(comp_str)
        i +=1```
        

*Examples:*

.. code:: ipython3

    compounds.percent_composition('KNO3')




.. parsed-literal::

    K          38.671
    N          13.856
    O          47.473
    Total:    100.000
    Name: KNO3, dtype: float64



.. code:: ipython3

    compounds.percent_composition_steps('KNO3')


.. parsed-literal::

    KNO3
    K: 39.1
    N: 14.01
    O: 3x16.0 = 48.0
    
    Mass KNO3 = 39.1 + 14.01 + 48.0 = 101.11
    
    K          39.10
    N          14.01
    O          48.00
    Total:    101.11
    Name: KNO3, dtype: float64
    K: 100 x 39.10g K/101.11g KNO3 = 38.67%
    N: 100 x 14.01g N/101.11g KNO3 = 13.86%
    O: 100 x 48.00g O/101.11g KNO3 = 47.47%
    

chemsolver.conversions 
~~~~~~~~~~~~~~~~~~~~~~~

-  **Importing the module**

.. code:: ipython3

    from chemsolver import conversions

-  **Functions**

***grams\_to\_moles***: Converts the mass of a sample of a substance to
the number of moles of that substance.

.. code:: python

    def grams_to_moles(mass, chemical):
        '''
        :param mass: mass of sample of substance (in grams)
        :param chemical:  chemical formula of substance
        :returns: moles of the sample of the substance
        '''
        molar_mass = compounds.mass(chemical)
        moles = round(mass / molar_mass, 2)
        return moles

***grams\_to\_moles\_steps***: Shows the step by step conversion of the
mass of a sample of a ubstance to the number of moles of that substance.

.. code:: python

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

*Examples*:

.. code:: ipython3

    conversions.grams_to_moles(36,'NO3')




.. parsed-literal::

    0.58



.. code:: ipython3

    conversions.grams_to_moles_steps(36,'NO3')


.. parsed-literal::

    36g NO3 x  1 mole NO3/62.01g NO3 = 0.58 moles NO3
    

***moles\_to\_grams***: From the number of moles of a sample of a
substance this function returns the mass of the sample in grams

.. code:: python

    def moles_to_grams(moles, chemical):
        """
        :param: moles: number of moles of the substance
        :param: chemical: the chemical symbol of the substance (as a string)
        """
        return moles * compounds.mass(chemical)

***moles\_to\_grams\_steps***: Displays the step by step conversion from
moles of a substance to grams of the substance

.. code:: python

    def moles_to_grams_steps(moles, chemical):
        grams = moles_to_grams(moles, chemical)
        compounds.mass_steps(chemical)
        mass = compounds.mass(chemical)
        print("{:.2f} moles x {:.2f} grams/mole = {:.2f} grams".format(moles, mass, grams))

*Examples*:

.. code:: ipython3

    conversions.moles_to_grams(moles = 1.5, chemical = 'H2O')




.. parsed-literal::

    27.03



.. code:: ipython3

    conversions.moles_to_grams_steps(moles = 1.5, chemical = 'H2O')


.. parsed-literal::

    H2O
    H: 2x1.008 = 2.016
    O: 16.0
    
    Mass H2O = 2.016 + 16.0 = 18.02
    1.50 moles x 18.02 grams/mole = 27.03 grams
    

***molarity***: Find the molarity (molar concentration) of a solute in a
solution.

.. code:: python

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

***molarity\_steps***: Demonstrates, step by step, the calculation of
the molarity of a solute in a solution.

.. code:: python

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

*Examples*:

.. code:: ipython3

    conversions.molarity(volume = 2, mass = 28.02, chemical = 'N2')




.. parsed-literal::

    0.5



.. code:: ipython3

    conversions.molarity_steps(volume = 2, mass = 28.02, chemical = 'N2')


.. parsed-literal::

    28.02g N2 x  1 mole N2/28.02g N2 = 1.0 moles N2
    1.00 moles/2.00 L = 0.50 M
    

.. code:: ipython3

    conversions.molarity(volume = 2, moles = 5)




.. parsed-literal::

    2.5



.. code:: ipython3

    conversions.molarity_steps(volume = 2, moles = 5)


.. parsed-literal::

    5.00 moles/2.00 L = 2.50 M
    

***molarity\_to\_moles***: Find the number of moles of a solute in a
solution from the molarity of the solution and the volume of the
solution (in liters).

.. code:: python

    def molarity_to_moles(M, volume):
        """
        :param M: the concentration of the solution in moles per L
        :param volume: the volume of the solution in Liters
        returns: moles of solute
        """
        return M * volume

***molarity\_to\_moles\_steps***: Shows the calculation of the number of
moles of a solute from the molarity of the solution and the volume of
the solution in liters.

.. code:: python

    def molarity_to_moles_steps(M, volume):
        print("{:.2f} moles/L x {:.2f} L = {:.2f} moles".format(M, volume, molarity_to_moles(M, volume)))

*Examples:*

.. code:: ipython3

    # M stands for molarity, volume is in Liters
    conversions.molarity_to_moles(M = 1.5, volume = 3.0 )




.. parsed-literal::

    4.5



.. code:: ipython3

    conversions.molarity_to_moles_steps(M = 1.5, volume = 3.0)


.. parsed-literal::

    1.50 moles/L x 3.00 L = 4.50 moles
    

***molarity\_to\_grams***: Find the mass of a solute, in grams, from the
molarity of the solution and the volume of the solution (in L).

.. code:: python

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

***molarity\_to\_grams\_steps***: Shows the step by step calculation of
the mass of a solute, in grams, from the molarity of the solution and
the volume of the solution (in L)

.. code:: python

    def molarity_to_grams_steps(molarity, volume, chemical):
        moles = molarity_to_moles(molarity, volume)
        moles_to_grams_steps(moles, chemical)

*Examples:*

.. code:: ipython3

    conversions.molarity_to_grams(molarity = 1.5, volume = 2, chemical = 'CH4')




.. parsed-literal::

    48.12



.. code:: ipython3

    conversions.molarity_to_grams_steps(molarity = 1.5, volume = 2, chemical = 'CH4')


.. parsed-literal::

    CH4
    C: 12.01
    H: 4x1.008 = 4.032
    
    Mass CH4 = 12.01 + 4.032 = 16.04
    3.00 moles x 16.04 grams/mole = 48.12 grams
    

chemsolver.ideal\_gas 
~~~~~~~~~~~~~~~~~~~~~~

-  **Importing the module**

.. code:: ipython3

    from chemsolver import ideal_gas

-  **Functions**

***ideal\_gas***: can be used to calculate the pressure (in atm), volume
(in liters), moles, mass (in grams), or temperature (in kelvin) of an
ideal gas.

.. code:: python

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

****ideal\_gas\_steps****: shows the step by step calculation of the
pressure, volume, moles, mass, or temperature of an ideal gas.

.. code:: python

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
            print('V = {} moles * (0.0821 L*atm/mol*K) * {} K / {} atm = {} L'.format(n, T, P, V))
        if solve_for == 'T':
            T = P * V / (n * R)
            print('T = PV/nR')
            print(ideal_gas_equation)
            print('T = {} atm * {} L/({} moles * 0.0821 L*atm/mole*K = {} K)'.format(P, V, n, T))

*Examples*:

.. code:: ipython3

    ideal_gas.ideal_gas(solve_for = 'P', V = 5, n = 4, T = 300)




.. parsed-literal::

    19.7



.. code:: ipython3

    ideal_gas.ideal_gas_steps(solve_for = "P", V = 5, n = 4, T = 300)


.. parsed-literal::

    PV = nRT
    P = nRT/V
    P = 4 moles * (0.0821 L*atm/mol*K) * 300K/5 L = 19.7 atm
    

.. code:: ipython3

    ideal_gas.ideal_gas(solve_for = 'V', substance = 'H2O', P = 1.5, m = 10.2, T = 300)




.. parsed-literal::

    9.36



.. code:: ipython3

    ideal_gas.ideal_gas_steps(solve_for = 'V',substance = 'H2O', P = 1.5, m = 10.2, T = 300)


.. parsed-literal::

    n = 10.2g * 1 mole H2O/18.02 g H2O = 0.57 moles H2O
    PV = nRT
    V = nRT/P
    V = 0.57 moles * (0.0821 L*atm/mol*K) * 300 K / 1.5 atm = 9.36 L
    
