# ChemSolver Package

## Purpose
I created this package to aid me in tutoring chemistry students for tutor.com. The package has many functions with the capability to solve simple chemistry problems. Each of these functions of a corresponding function with the same name with **"_steps"** tacked on a at the end of the function's name that provided a step by step solution to the problem or step. The goal of this package is to enable me to quickly check student's work and also display solutions, while being able to focus on checking the student's understanding and catering to his or her needs.

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
                    return 'Element: ' + self.name + ', ' + self.symbol  +'\nAtomic Number: ' +str(self.atomic_number) + '\nMass: {}'.format(self.mass)
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
    - __*mass*__: returns the molar mass (or mass in atm) of chemical compounds
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
    
    - __*mass_steps*__: shows a step by step calculation of the molar mass of a compound
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
     
    - __*composition*__: returns the mass composition of chemical elements.
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
     
 
