# -*- coding: utf-8 -*-
"""
Created on Wed Jul 26 20:27:31 2017

@author: rosen
"""

import pandas as pd
import numpy as np

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



H = Element('H','Hydrogen',1,1.008)
He = Element('He','Helium',2,4.0026)
Li = Element('Li','Lithium',3,6.94)
Be = Element('Be','Berylium',4,9.01)
B = Element('B','Boron',5,10.81)
C = Element('C','Carbon',6,12.01)
N = Element('N','Nitrogen',7,14.01)
O = Element('O','Oxygen',8, 16.00)
F = Element('F','Fluorine',9,19.00)
Ne = Element('Ne','Neon',10,20.18)
Na = Element('Na','Sodium',11,22.99)
Mg = Element('Mg','Magnesium',12,24.31)
Al = Element('Al','Aluminum',13,26.98)
Si = Element('Si','Silicon',14,28.09)
P = Element('P','Phosphorus',15,30.97)
S = Element('S','Sulfur',16,32.07)
Cl = Element('Cl','Chlorine',17,35.45)
Ar = Element('Ar','Argon',18,39.95)
K = Element('K','Potassium',19,39.10)
Ca = Element('Ca','Calcium',20,40.08)
Sc = Element('Sc','Scandium',21,44.96)
Ti = Element('Ti','Titanium',22,47.87)
V = Element('V','Vanadium',23,50.94)
Cr = Element('Cr','Chromium',24,52.00)
Mn = Element('Mn','Manganese',25,54.94)
Fe = Element('Fe','Iron',26,55.85)
Co = Element('Co','Cobalt',27,58.93)
Ni = Element('Ni','Nickel',28,58.69)
Cu = Element('Cu','Copper',29,63.55)
Zn = Element('Zn','Zinc',30,65.39)
Ga = Element('Ga','Gallium',31,69.72)
Ge = Element('Ge','Germanium',32,72.64)
As = Element('As','Arsenic',33,74.92)
Se = Element('Se','Selenium',34,78.96)
Br = Element('Br','Bromine',35,79.90)
Kr = Element('Kr','Krypton',36,83.8)
Rb = Element('Rb','Rubidium',37,85.47)
Sr = Element('Sr','Strontium',38,87.62)
Y = Element('Y','Yttrium',39,88.91)
Zr = Element('Zr','Zirconium',40,91.22)
Nb = Element('Nb','Niobium',41,92.91)
Mo = Element('Mo','Molybdenum',42,95.94)
Tc = Element('TC','Technetium',43,98)
Ru = Element('Ru','Ruthenium',44,101.07)
Rh = Element('Rh','Rhodium',45,102.91)
Pd = Element('Pd','Palladium',46,106.42)
Ag = Element('Ag','Silver',47,107.87)
Cd = Element('Cd','Cadmium',48,112.411)
In = Element('In','Indium',49,114.82)
Sn = Element('Sn','Tin',50,118.71)
Sb = Element('Sb','Antimony',51,121.76)
Te = Element('Te','Tellurium',52,127.6)
I = Element('I','Iodine',53,126.90)
Xe = Element('Xe','Xenon',54,131.29)
Cs = Element('Cs','Cesium',55,132.91)
Ba = Element('Ba','Barium',56,137.33)
La = Element('La','Lanthanum',57,138.91)
Ce = Element('Ce','Cerium',58,140.12)
Pr = Element('Pr','Praseodymium',59,140.91)
Nd = Element('Nd','Neodymium',60,144.24)
Pm = Element('Pm','Promethium',61,145)
Sm = Element('Sm','Samarium',62,150.36)
Eu = Element('Eu','Europium',63,151.96)
Gd = Element('Gd','Gadolinium',64,157.25)
Tb = Element('Tb','Terbium',65,158.93)
Dy = Element('Dy','Dysprosium',66,162.5)
Ho = Element('Ho','Holmim',67,164.93)
Er = Element('Er','Erbium',68,167.26)
Tm = Element('Tm','Thulium',69,168.93)
Yb = Element('Yb','Ytterbium',70,173.04)
Lu = Element('Lu','Lutetium',71,174.97)
Hf = Element('Hf','Hafnium',72,178.49)
Ta = Element('Ta','Tantalum',73,180.95)
W = Element('W','Tungsten',74,183.84)
Re = Element('Re','Rhenium',75,186.21)
Os = Element('Os','Osmium',76,190.23)
Ir = Element('Ir','Iridium',77,192.22)
Pt = Element('Pt','Platinum',78,195.08)
Au = Element('Au','Gold',79,196.97)
Hg = Element('Hg','Mercury',80,200.59)
Tl = Element('Tl','Thallium',81,204.38)
Pb = Element('Pb','Lead',82,207.2)
Bi = Element('Bi','Bismuth',83,208.98)
Po = Element('Po','Polonium',84,209)
At = Element('At','Astatine',85,210)
Rn = Element('Rn','Radon',86,222)
Fr = Element('Fr','Francium',87,223)
Ra = Element('Ra','Radium',88,226)
Ac = Element('Ac','Actinium',89,227)
Th = Element('Th','Thorium',90,232.04)
Pa = Element('Pa','Proctactinium',91,231.04)
U = Element('U','Uranium',92,238.03)


element_index = ['H','He','Li','Be','B','C','N','O','F','Ne','Na','Mg','Al','Si','P','S','Cl','Ar','K','Ca','Sc','Ti','V','Cr','Mn','Fe','Co','Ni',"Cu","Zn","Ga","Ge","As","Se","Br","Kr","Rb","Sr","Y","Zr","Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd","In","Sn","Sb","Te","I","Xe","Cs","Ba","La","Ce","Pr","Nd","Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb","Lu","Hf","Ta","W","Re","Os","Ir","Pt","Au","Hg","Tl","Pb","Bi","Po","At","Rn","Fr","Ra","Ac","Th","Pa","U"]

element_objects = [H,He,Li,Be,B,C,N,O,F,Ne,Na,Mg,Al,Si,P,S,Cl,Ar,K,Ca,Sc,Ti,V,Cr,Mn,Fe,Co,Ni,Cu,Zn,Ga,Ge,As,Se,Br,Kr,Rb,Sr,Y,Zr,Nb,Mo,Tc,Ru,Rh,Pd,Ag,Cd,In,Sn,Sb,Te,I,Xe,Cs,Ba,La,Ce,Pr,Nd,Pm,Sm,Eu,Gd,Tb,Dy,Ho,Er,Tm,Yb,Lu,Hf,Ta,W,Re,Os,Ir,Pt,Au,Hg,Tl,Pb,Bi,Po,At,Rn,Fr,Ra,Ac,Th,Pa,U]

elements = pd.Series(element_objects,index=element_index)