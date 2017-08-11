
# ChemHelper Package: Version 0.4

## Purpose
I created this package to aid me in tutoring chemistry students for tutor.com. The package has many functions with the capability to solve simple chemistry problems. Each of these functions of a corresponding function with the same name with **"_steps"** tacked on a at the end of the function's name that provided a step by step solution to the problem or step. The goal of this package is to enable me to quickly check student's work and also display solutions, while being able to focus on checking the student's understanding and catering to his or her needs.

## Installation

```ruby
pip install chemhelper
``` 

or available for download at <a href = "https://pypi.python.org/pypi/chemhelper/0.4">https://pypi.python.org/pypi/chemhelper/0.4</a>

## Modules 
- **<a href="#chemhelper.elements" style="text-decoration:none">chemhelper.elements</a>**

    - <a href="#chemhelper.elements.import" style="text-decoration:none">Importing the module</a>

    - <a href="#chemhelper.elements.class" style="text-decoration:none">Element Class</a>
            
        - <a href="#chemhelper.elements.objects" style="text-decoration:none">Element Objects</a>
            - <a href="#chemhelper.elements.objects.attributes" style="text-decoration:none">Element Object Attributes</a>
            - <a href="#chemhelper.elements.objects.examples" style="text-decoration:none">Examples</a>
            <br><br>
     
    - <a href="#chemhelper.elements.element_index" style="text-decoration:none">element_index</a>
    - <a href="#chemhelper.elements.element_objects" style="text-decoration:none">element_objects</a>
    - <a href="#chemhelper.elements.elements" style="text-decoration:none">elements</a> <br> <br>
- **<a href = "#chemhelper.compounds" style="text-decoration:none">chemhelper.compounds</a>**
    - <a href = "#chemhelper.compounds.import" style="text-decoration:none">Importing the module</a>
    - <a href = "#chemhelper.compounds.functions" style="text-decoration:none">Functions</a>
        - <a href = "#chemhelper.compounds.functions.mass" style="text-decoration:none">mass</a>
        - <a href = "#chemhelper.compounds.functions.mass_steps" style="text-decoration:none">mass_steps</a>
            - <a href="#chemhelper.compounds.functions.examples1" style="text-decoration:none">mass and mass_steps examples</a>
            <br><br>
        - <a href = "#chemhelper.compounds.functions.composition" style="text-decoration:none">composition</a>
        - <a href = "#chemhelper.compounds.functions.composition_steps" style="text-decoration:none">composition_steps</a>
            - <a href="#chemhelper.compounds.functions.examples2" style="text-decoration:none">composition and composition_steps examples</a>
                <br><br>
        - <a href = "#chemhelper.compounds.functions.percent_composition" style="text-decoration:none">percent_composition</a>
        - <a href = "#chemhelper.compounds.functions.percent_composition_steps" style="text-decoration:none">percent_composition_steps</a>
            - <a href="#chemhelper.compounds.functions.examples3" style="text-decoration:none">percent_composition and percent_composition_steps examples</a> <br><br>
- **<a href = "#chemhelper.conversions" style="text-decoration:none">chemhelper.conversions</a>**
    - <a href = "#chemhelper.conversions.import" style="text-decoration:none">Importing the module</a>
    - <a href = "#chemhelper.conversions.functions" style="text-decoration:none">Functions</a>
        - <a href = "#chemhelper.conversions.functions.grams_to_moles" style="text-decoration:none">grams_to_moles</a>
        - <a href = "#chemhelper.conversions.functions.grams_to_moles_steps" style="text-decoration:none">grams_to_moles_steps</a>
            - <a href="#chemhelper.conversions.functions.examples1" style="text-decoration:none">grams_to_moles and grams_to_moles_steps examples</a>
            <br><br>
        - <a href = "#chemhelper.conversions.functions.moles_to_grams" style="text-decoration:none">moles_to_grams</a>
        - <a href = "#chemhelper.conversions.functions.moles_to_grams_steps" style="text-decoration:none">moles_to_grams_steps</a>
            - <a href = "#chemhelper.conversions.functions.examples2" style="text-decoration:none">moles_to grams and moles_to_grams_steps examples</a>
            <br><br>
        - <a href = "#chemhelper.conversions.functions.molarity" style="text-decoration:none">molarity</a>
        - <a href = "#chemhelper.conversions.functions.molarity_steps" style="text-decoration:none">molarity_steps</a>
            - <a href = "#chemhelper.conversions.functions.examples3" style="text-decoration:none">molarity and molarity_steps examples</a>
             <br><br>
        - <a href = "#chemhelper.conversions.functions.molarity_to_moles" style="text-decoration:none">molarity_to_moles</a>
        - <a href = "#chemhelper.conversions.functions.molarity_to_moles_steps" style="text-decoration:none">molarity_to_moles_steps</a>
            - <a href = "#chemhelper.conversions.functions.examples4" style="text-decoration:none">molarity_to_moles and molarity_to_moles_steps examples</a>
             <br><br>
        - <a href = "#chemhelper.conversions.functions.molarity_to_grams" style="text-decoration:none">molarity_to_grams</a>
        - <a href = "#chemhelper.conversions.functions.molarity_to_grams_steps" style="text-decoration:none">molarity_to_grams_steps</a>
            - <a href = "#chemhelper.conversions.functions.examples5" style="text-decoration:none">molarity and molarity_steps examples</a>        
            
            
- **<a href = "#chemhelper.ideal_gas" style="text-decoration:none">chemhelper.ideal_gas</a>**
    - <a href = "#chemhelper.ideal_gas.import" style="text-decoration:none">Importing the module</a>
    - <a href = "#chemhelper.ideal_gas.functions" style="text-decoration:none">Functions</a>
        - <a href = "#chemhelper.ideal_gas.functions.ideal_gas" style="text-decoration:none">ideal_gas</a>
        - <a href = "#chemhelper.ideal_gas.functions.ideal_gas_steps" style="text-decoration:none">ideal_gas_steps</a>
            - <a href="#chemhelper.ideal_gas.functions.examples" style="text-decoration:none">ideal_gas and ideal_gas_steps examples</a>

### chemhelper.elements <a name="chemhelper.elements"></a>

- **Importing the module**<a name="chemhelper.elements.import"></a>


```python
from chemhelper import elements
```

- **Element Class** <a name="chemhelper.elements.class"></a>
   - The elements module contains an **Element object** for each of the 93 naturally occuring elements.
```python
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
```
               
- **Element Objects:**  <a name="chemhelper.elements.objects"></a> The variable names of the elements are the chemical symbols of the elements (example: hydrogen is H)
    - **Attributes**: <a name = "chemhelper.elements.objects.attributes"></a> 
        - *self.symbol*: The chemical symbol of the element (example
        - *self.name*: The name of the element
        - *self.mass*: The molar mass of the element (also numerically equal to the mass of 1 atom of the element in atm)
        - *self.atomic_number*: The atomic number of the element  
          
__*Examples*__: <a name = "chemhelper.elements.objects.examples"></a>


```python
print(elements.H)
```

    Element: Hydrogen, H
    Atomic Number: 1
    Mass: 1.008
    


```python
elements.Ti.name
```




    'Titanium'




```python
elements.Ag.mass
```




    107.87




```python
elements.K.atomic_number
```




    19




```python
print(elements.Ca)
```

    Element: Calcium, Ca
    Atomic Number: 20
    Mass: 40.08
    

   - **element_index:** List of element symbols. <a name = "chemhelper.elements.element_index"></a>
    
```python
      element_index = ["Ni","Cu","Zn","Ga","Ge","As","Se","Br","Kr","Rb","Sr","Y","Zr","Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd","In","Sn","Sb","Te","I","Xe","Cs","Ba","La","Ce","Pr","Nd","Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb","Lu","Hf","Ta","W","Re","Os","Ir","Pt","Au","Hg","Tl","Pb","Bi","Po","At","Rn","Fr","Ra","Ac","Th","Pa","U"]```
      
      
- **element_objects:** List of element objects <a name ="chemhelper.elements.element_objects"></a>
```python
        element_objects = [H,He,Li,Be,B,C,N,O,F,Ne,Na,Mg,Al,Si,P,S,Cl,Ar,K,Ca,Sc,Ti,V,Cr,Mn,Fe,Co,Ni,Cu,Zn,Ga,Ge,As,Se,Br,Kr,Rb,Sr,Y,Zr,Nb,Mo,Tc,Ru,Rh,Pd,Ag,Cd,In,Sn,Sb,Te,I,Xe,Cs,Ba,La,Ce,Pr,Nd,Pm,Sm,Eu,Gd,Tb,Dy,Ho,Er,Tm,Yb,Lu,Hf,Ta,W,Re,Os,Ir,Pt,Au,Hg,Tl,Pb,Bi,Po,At,Rn,Fr,Ra,Ac,Th,Pa,U]
```        
- **elements**: <a name ="chemhelper.elements.elements"></a> a pandas series containing every element object.

### chemhelper.compounds <a name = "chemhelper.compounds"></a>

- **Importing the module** <a name="chemhelper.compounds.import"></a>


```python
from chemhelper import compounds
```

- **Functions** <a name = "chemhelper.compounds.functions"></a>

    __*mass*__: <a name = "chemhelper.compounds.functions.mass"></a>determines the molar mass of a molecule and/or compound
     
```python 
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
    return round(m,2)```
 
    
   
   
   
   __*mass_steps*__: <a name = "chemhelper.compounds.functions.mass_steps"></a> shows the calculation of the the molar mass of a molecule and/or compound
     
```python 
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
    ```
    
*Examples*: <a name = "chemhelper.compounds.functions.examples1"></a>


```python
compounds.mass('HNO4')
```




    79.02




```python
compounds.mass_steps('HNO4')
```

    HNO4
    H: 1.008
    N: 14.01
    O: 4x16.0 = 64.0
    
    Mass HNO4 = 1.008 + 14.01 + 64.0 = 79.02
    

   __**composition**__: <a name = "chemhelper.compounds.functions.composition"></a> returns the composition of compounds by mass
     
```python
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
    ```
   
 __**composition_steps**__: <a name = "chemhelper.compounds.functions.composition_steps"></a> step by step solution of the chemical composition of a compound by mass
    
```python
def composition_steps(chemical):
    '''
    :param chemical: chemical compound symbol (example 'H2O')
    :prints: step by step solution of mass due to each element in compound and total mass
    '''
    mass_steps(chemical)
    print('')
    print(composition(chemical))
```

*Examples:* <a name = "chemhelper.compounds.functions.examples2"></a>


```python
compounds.composition('NH4')
```




    N         14.010
    H          4.032
    Total:    18.040
    Name: NH4, dtype: float64




```python
compounds.composition_steps('NH4')
```

    NH4
    N: 14.01
    H: 4x1.008 = 4.032
    
    Mass NH4 = 14.01 + 4.032 = 18.04
    
    N         14.010
    H          4.032
    Total:    18.040
    Name: NH4, dtype: float64
    

__**percent_composition**__: <a name = "chemhelper.compounds.functions.percent_composition"></a> calculates the percent composition of each element in a compound.
```python
def percent_composition(chemical):
    '''
    :param chemical: chemical compound symbol (example 'H2O')
    :return: pandas Series containing percent composition of each element making up compound.
    '''
    return round((composition(chemical)/mass(chemical))* 100,3)
```

__**percent_composition_steps**__: <a name = "chemhelper.compounds.functions.percent_composition_steps"></a> calculates and shoes, step by step, the percent composition of each element in a compound
```python
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
        i +=1```
        
*Examples:* <a name = "chemhelper.compounds.functions.examples3"></a>


```python
compounds.percent_composition('KNO3')
```




    K          38.671
    N          13.856
    O          47.473
    Total:    100.000
    Name: KNO3, dtype: float64




```python
compounds.percent_composition_steps('KNO3')
```

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
    

### chemhelper.conversions <a name = "chemhelper.conversions"></a>

- **Importing the module** <a name = "chemhelper.conversions.import"></a>


```python
from chemhelper import conversions
```

- **Functions** <a name = "chemhelper.conversions.functions"></a>

__*grams_to_moles*__: <a name = "chemhelper.conversions.functions.grams_to_moles"></a> Converts the mass of a sample of a substance to the number of moles of that substance.

    
```python
def grams_to_moles(mass, chemical):
    '''
    :param mass: mass of sample of substance (in grams)
    :param chemical:  chemical formula of substance
    :returns: moles of the sample of the substance
    '''
    molar_mass = compounds.mass(chemical)
    moles = round(mass / molar_mass, 2)
    return moles
```

__*grams_to_moles_steps*__:  <a name = "chemhelper.conversions.functions.grams_to_moles_steps"></a> Shows the step by step conversion of the mass of a sample of a ubstance to the number of moles of that substance.

```python
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
```


*Examples*: <a name = "chemhelper.conversions.functions.examples1"></a>


```python
conversions.grams_to_moles(36,'NO3')
```




    0.58




```python
conversions.grams_to_moles_steps(36,'NO3')
```

    36g NO3 x  1 mole NO3/62.01g NO3 = 0.58 moles NO3
    

__*moles_to_grams*__: <a name = "chemhelper.conversions.functions.moles_to_grams"></a> From the number of moles of a sample of a substance this function returns the mass of the sample in grams
```python
def moles_to_grams(moles, chemical):
    """
    :param: moles: number of moles of the substance
    :param: chemical: the chemical symbol of the substance (as a string)
    """
    return moles * compounds.mass(chemical)
```

__*moles_to_grams_steps*__: <a name = "chemhelper.conversions.functions.moles_to_grams_steps"></a> Displays the step by step conversion from moles of a substance to grams of the substance 

```python
def moles_to_grams_steps(moles, chemical):
    grams = moles_to_grams(moles, chemical)
    compounds.mass_steps(chemical)
    mass = compounds.mass(chemical)
    print("{:.2f} moles x {:.2f} grams/mole = {:.2f} grams".format(moles, mass, grams))
```

*Examples*: <a name = "chemhelper.conversions.functions.examples2"></a>


```python
conversions.moles_to_grams(moles = 1.5, chemical = 'H2O')
```




    27.03




```python
conversions.moles_to_grams_steps(moles = 1.5, chemical = 'H2O')
```

    H2O
    H: 2x1.008 = 2.016
    O: 16.0
    
    Mass H2O = 2.016 + 16.0 = 18.02
    1.50 moles x 18.02 grams/mole = 27.03 grams
    

__*molarity*__: <a name = "chemhelper.conversions.functions.molarity"></a> Find the molarity (molar concentration) of a solute in a solution. 

```python
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
```

__*molarity_steps*__: <a name = "chemhelper.conversions.functions.molarity_steps"></a> Demonstrates, step by step, the calculation of the molarity of a solute in a solution.

```python
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
```

*Examples*: <a name = "chemhelper.conversions.functions.examples3"></a>


```python
conversions.molarity(volume = 2, mass = 28.02, chemical = 'N2')
```




    0.5




```python
conversions.molarity_steps(volume = 2, mass = 28.02, chemical = 'N2')
```

    28.02g N2 x  1 mole N2/28.02g N2 = 1.0 moles N2
    1.00 moles/2.00 L = 0.50 M
    


```python
conversions.molarity(volume = 2, moles = 5)
```




    2.5




```python
conversions.molarity_steps(volume = 2, moles = 5)
```

    5.00 moles/2.00 L = 2.50 M
    

__*molarity_to_moles*__: <a name = "chemhelper.conversions.functions.molarity_to_moles"></a> Find the number of moles of a solute in a solution from the molarity of the solution and the volume of the solution (in liters).

```python
def molarity_to_moles(M, volume):
    """
    :param M: the concentration of the solution in moles per L
    :param volume: the volume of the solution in Liters
    returns: moles of solute
    """
    return M * volume
```

__*molarity_to_moles_steps*__: <a name = "chemhelper.conversions.functions.molarity_to_moles_steps"></a> Shows the calculation of the number of moles of a solute from the molarity of the solution and the volume of the solution in liters.

```python
def molarity_to_moles_steps(M, volume):
    print("{:.2f} moles/L x {:.2f} L = {:.2f} moles".format(M, volume, molarity_to_moles(M, volume)))
```

*Examples:* <a name = "chemhelper.conversions.functions.examples4"></a>


```python
# M stands for molarity, volume is in Liters
conversions.molarity_to_moles(M = 1.5, volume = 3.0 )
```




    4.5




```python
conversions.molarity_to_moles_steps(M = 1.5, volume = 3.0)
```

    1.50 moles/L x 3.00 L = 4.50 moles
    

__*molarity_to_grams*__: <a name = "chemhelper.conversions.functions.molarity_to_grams"></a> Find the mass of a solute, in grams, from the molarity of the solution and the volume of the solution (in L).
```python
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
```

__*molarity_to_grams_steps*__: <a name = "chemhelper.conversions.functions.molarity_to_grams_steps"></a> Shows the step by step calculation of the mass of a solute, in grams, from the molarity of the solution and the volume of the solution (in L)
```python
def molarity_to_grams_steps(molarity, volume, chemical):
    moles = molarity_to_moles(molarity, volume)
    moles_to_grams_steps(moles, chemical)
```

*Examples:* <a name = "chemhelper.conversions.functions.examples5"></a>


```python
conversions.molarity_to_grams(molarity = 1.5, volume = 2, chemical = 'CH4')
```




    48.12




```python
conversions.molarity_to_grams_steps(molarity = 1.5, volume = 2, chemical = 'CH4')
```

    CH4
    C: 12.01
    H: 4x1.008 = 4.032
    
    Mass CH4 = 12.01 + 4.032 = 16.04
    3.00 moles x 16.04 grams/mole = 48.12 grams
    

### chemhelper.ideal_gas <a name = "chemhelper.ideal_gas"></a>

- **Importing the module** <a name = "chemhelper.ideal_gas.import"></a>


```python
from chemhelper import ideal_gas
```

- **Functions** <a name = "chemhelper.ideal_gas.functions"></a>

__*ideal_gas*__: <a name = "chemhelper.ideal_gas.functions.ideal_gas"></a> can be used to calculate the pressure (in atm), volume (in liters), moles, mass (in grams), or temperature (in kelvin) of an ideal gas.
```python
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
```

__**ideal_gas_steps**__: <a name = "chemhelper.ideal_gas.functions.ideal_gas_steps"></a> shows the step by step calculation of the pressure, volume, moles, mass, or temperature of an ideal gas.

```python
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
```

*Examples*: <a name = 'chemhelper.ideal_gas.functions.examples'></a>


```python
ideal_gas.ideal_gas(solve_for = 'P', V = 5, n = 4, T = 300)
```




    19.7




```python
ideal_gas.ideal_gas_steps(solve_for = "P", V = 5, n = 4, T = 300)
```

    PV = nRT
    P = nRT/V
    P = 4 moles * (0.0821 L*atm/mol*K) * 300K/5 L = 19.7 atm
    


```python
ideal_gas.ideal_gas(solve_for = 'V', substance = 'H2O', P = 1.5, m = 10.2, T = 300)
```




    9.36




```python
ideal_gas.ideal_gas_steps(solve_for = 'V',substance = 'H2O', P = 1.5, m = 10.2, T = 300)
```

    n = 10.2g * 1 mole H2O/18.02 g H2O = 0.57 moles H2O
    PV = nRT
    V = nRT/P
    V = 0.57 moles * (0.0821 L*atm/mol*K) * 300 K / 1.5 atm = 9.36 L
    
