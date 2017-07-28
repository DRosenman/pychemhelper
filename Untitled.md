
# ChemSolver Package

## Purpose
I created this package to aid me in tutoring chemistry students for tutor.com. The package has many functions with the capability to solve simple chemistry problems. Each of these functions of a corresponding function with the same name with **"\_steps"** tacked on a at the end of the function's name that provided a step by step solution to the problem or step. The goal of this package is to enable me to quickly check student's work and also display solutions, while being able to focus on checking the student's understanding and catering to his or her needs.

## Modules
### chemsolver.elements

- **Element Class**
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
                    return ('Element: {}, {}\nAtomic Number: {}\nMass:
                   {}'.format(self.name,self.symbol,self.atomic_number,self.mass)||
                   ```
               
   - **Element Objects:** The variable names of the elements are the chemical symbols of the elements (example: hydrogen is H)
       - **Attributes**: 
          - *self.symbol*: The chemical symbol of the element (example
          - *self.name*: The name of the element
          - *self.mass*: The molar mass of the element (also numerically equal to the mass of 1 atom of the element in atm)
          - *self.atomic_number*: The atomic number of the element  
      
- **element_index:** List of element symbols.
    - List of element symbols
```python
      element_index = ["Ni","Cu","Zn","Ga","Ge","As","Se","Br","Kr","Rb","Sr","Y","Zr","Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd","In","Sn","Sb","Te","I","Xe","Cs","Ba","La","Ce","Pr","Nd","Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb","Lu","Hf","Ta","W","Re","Os","Ir","Pt","Au","Hg","Tl","Pb","Bi","Po","At","Rn","Fr","Ra","Ac","Th","Pa","U"]```
      
      
- **element_objects:** List of element objects
```python
        element_objects = [H,He,Li,Be,B,C,N,O,F,Ne,Na,Mg,Al,Si,P,S,Cl,Ar,K,Ca,Sc,Ti,V,Cr,Mn,Fe,Co,Ni,Cu,Zn,Ga,Ge,As,Se,Br,Kr,Rb,Sr,Y,Zr,Nb,Mo,Tc,Ru,Rh,Pd,Ag,Cd,In,Sn,Sb,Te,I,Xe,Cs,Ba,La,Ce,Pr,Nd,Pm,Sm,Eu,Gd,Tb,Dy,Ho,Er,Tm,Yb,Lu,Hf,Ta,W,Re,Os,Ir,Pt,Au,Hg,Tl,Pb,Bi,Po,At,Rn,Fr,Ra,Ac,Th,Pa,U]
        
- **elements**: a pandas series containing every element object.

### chemsolver.compounds

- **Functions**
    - *mass*: returns the molar mass (or mass in atm) of chemical compounds
    ```python
    def mass(chemical):
    '''
    :param chemical: compound chemical symbol (example: 'H2O')
    :return: molar mass of compound
    '''
    chemical = re.findall(r'([A-Z][a-z]*)(\d*)', chemical)
    m = 0
    for tup in chemical:
        if tup[1] == '': m+=elements. element_series[tup[0]].mass
        else:
            m+= elements. element_series[tup[0]].mass * int(tup[1])
    return round(m,2)
    ```
    
    - *mass_steps*: shows a step by step calculation of the molar mass of a compound
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
            m+= elements.element[tup[0]].mass * int(tup[1])
            numbers.append(elements.elements[tup[0]].mass * int(tup[1]))
    n = len(numbers)
    i = 1

    string = str(numbers[0])
    while i < n:
        string = string + ' + ' +  str(numbers[i ])
        i+=1 
    print( '\nMass ' + chemical + ' = ' + string + ' = ' + str(round(m,2)))
    ```
     
    - *composition*: returns the mass composition of chemical elements.
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
     
 


```python
import elements
```


```python
elements.elements
```




    H     Element: Hydrogen, H\nAtomic Number: 1\nMass: ...
    He    Element: Helium, He\nAtomic Number: 2\nMass: 4...
    Li    Element: Lithium, Li\nAtomic Number: 3\nMass: ...
    Be    Element: Berylium, Be\nAtomic Number: 4\nMass:...
    B      Element: Boron, B\nAtomic Number: 5\nMass: 10.81
    C     Element: Carbon, C\nAtomic Number: 6\nMass: 12.01
    N     Element: Nitrogen, N\nAtomic Number: 7\nMass: ...
    O      Element: Oxygen, O\nAtomic Number: 8\nMass: 16.0
    F     Element: Fluorine, F\nAtomic Number: 9\nMass: ...
    Ne    Element: Neon, Ne\nAtomic Number: 10\nMass: 20.18
    Na    Element: Sodium, Na\nAtomic Number: 11\nMass: ...
    Mg    Element: Magnesium, Mg\nAtomic Number: 12\nMas...
    Al    Element: Aluminum, Al\nAtomic Number: 13\nMass...
    Si    Element: Silicon, Si\nAtomic Number: 14\nMass:...
    P     Element: Phosphorus, P\nAtomic Number: 15\nMas...
    S     Element: Sulfur, S\nAtomic Number: 16\nMass: 3...
    Cl    Element: Chlorine, Cl\nAtomic Number: 17\nMass...
    Ar    Element: Argon, Ar\nAtomic Number: 18\nMass: 3...
    K     Element: Potassium, K\nAtomic Number: 19\nMass...
    Ca    Element: Calcium, Ca\nAtomic Number: 20\nMass:...
    Sc    Element: Scandium, Sc\nAtomic Number: 21\nMass...
    Ti    Element: Titanium, Ti\nAtomic Number: 22\nMass...
    V     Element: Vanadium, V\nAtomic Number: 23\nMass:...
    Cr    Element: Chromium, Cr\nAtomic Number: 24\nMass...
    Mn    Element: Manganese, Mn\nAtomic Number: 25\nMas...
    Fe    Element: Iron, Fe\nAtomic Number: 26\nMass: 55.85
    Co    Element: Cobalt, Co\nAtomic Number: 27\nMass: ...
    Ni    Element: Nickel, Ni\nAtomic Number: 28\nMass: ...
    Cu    Element: Copper, Cu\nAtomic Number: 29\nMass: ...
    Zn    Element: Zinc, Zn\nAtomic Number: 30\nMass: 65.39
                                ...                        
    Eu    Element: Europium, Eu\nAtomic Number: 63\nMass...
    Gd    Element: Gadolinium, Gd\nAtomic Number: 64\nMa...
    Tb    Element: Terbium, Tb\nAtomic Number: 65\nMass:...
    Dy    Element: Dysprosium, Dy\nAtomic Number: 66\nMa...
    Ho    Element: Holmim, Ho\nAtomic Number: 67\nMass: ...
    Er    Element: Erbium, Er\nAtomic Number: 68\nMass: ...
    Tm    Element: Thulium, Tm\nAtomic Number: 69\nMass:...
    Yb    Element: Ytterbium, Yb\nAtomic Number: 70\nMas...
    Lu    Element: Lutetium, Lu\nAtomic Number: 71\nMass...
    Hf    Element: Hafnium, Hf\nAtomic Number: 72\nMass:...
    Ta    Element: Tantalum, Ta\nAtomic Number: 73\nMass...
    W     Element: Tungsten, W\nAtomic Number: 74\nMass:...
    Re    Element: Rhenium, Re\nAtomic Number: 75\nMass:...
    Os    Element: Osmium, Os\nAtomic Number: 76\nMass: ...
    Ir    Element: Iridium, Ir\nAtomic Number: 77\nMass:...
    Pt    Element: Platinum, Pt\nAtomic Number: 78\nMass...
    Au    Element: Gold, Au\nAtomic Number: 79\nMass: 19...
    Hg    Element: Mercury, Hg\nAtomic Number: 80\nMass:...
    Tl    Element: Thallium, Tl\nAtomic Number: 81\nMass...
    Pb    Element: Lead, Pb\nAtomic Number: 82\nMass: 207.2
    Bi    Element: Bismuth, Bi\nAtomic Number: 83\nMass:...
    Po    Element: Polonium, Po\nAtomic Number: 84\nMass...
    At    Element: Astatine, At\nAtomic Number: 85\nMass...
    Rn     Element: Radon, Rn\nAtomic Number: 86\nMass: 222
    Fr    Element: Francium, Fr\nAtomic Number: 87\nMass...
    Ra    Element: Radium, Ra\nAtomic Number: 88\nMass: 226
    Ac    Element: Actinium, Ac\nAtomic Number: 89\nMass...
    Th    Element: Thorium, Th\nAtomic Number: 90\nMass:...
    Pa    Element: Proctactinium, Pa\nAtomic Number: 91\...
    U     Element: Uranium, U\nAtomic Number: 92\nMass: ...
    Length: 92, dtype: object




```python

```
