import pandas as pd
import os
import matplotlib.pyplot as plt
import math
import numpy as np
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 unused import

import math
import numpy as np

dir_path = os.path.dirname(os.path.realpath(__file__))

E = 72 # GPa
G = 28 # GPa
v = E / (2 * G) - 1

inner_d = 0.01250 # m
outer_d = 0.02750

area = math.pi / 4 * (outer_d ** 2 - inner_d ** 2)
I = np.pi / 4 * ((outer_d/2) ** 4 - (inner_d/2) ** 4)
J = np.pi / 2 * ((outer_d/2) ** 4 - (inner_d/2) ** 4)

class Data: # Defines the Data class, which contains all of the data from a given .csv file
    def __init__(self, length):

        df = pd.read_csv(dir_path + '\\' + length + r'cm_test.csv')
        self.length = length + ' cm'
        self.dec_length = float(length) / 100 # length in m

        self.load = df['P (N)']

        e0 = df['microstrain3'] # Pulls data only once for computational speed
        e45 = df['microstrain4']
        e90 = df['microstrain5']
        self.exx = e0
        self.eyy = e90
        self.exy = e45 - (e0 + e90)/2

        self.sxx = 1000 * E*(self.exx * 10 ** -6 + v*self.eyy * 10 ** -6)/(1 - v**2) # Calculates stress values, in MPa
        self.txy = 1000 * 2*G*self.exy * 10 ** -6
        self.syy = 1000 * E * (self.eyy * 10 ** -6 + v * self.exx * 10 ** -6) / (1 - v ** 2)

        self.theor_sxx = self.load * 0.11250 * (outer_d/2) / I * 10 ** -6
        self.theor_txy = self.load * self.dec_length * outer_d/2 / J * 10 ** -6
        self.theor_syy = -self.load / area * 10 ** -6

        self.p1 = (self.sxx + self.syy)/2 + ((self.sxx-self.syy) ** 2 /4 + self.txy ** 2) ** 0.5

        self.p1_theor = (self.theor_sxx + self.theor_syy) / 2 + \
                    ((self.theor_sxx - self.theor_syy) ** 2 / 4 + self.theor_txy ** 2) ** 0.5

        self.theta_1 = np.arctan2(2 * self.txy, self.sxx - self.syy) / 2 * 180/math.pi
        self.theta_1_theor = np.arctan2(2 * self.theor_txy, self.theor_sxx - self.theor_syy) / 2 * 180/math.pi
        # Value (in degrees) checked to correspond to p1

L7cm = Data('7') # Imports .csv file data
L10cm = Data('10')
L15cm = Data('15')

### c ----------------------------------

plt_array = [L7cm, L10cm, L15cm]

for i in plt_array: # For each
    plt.figure()
    exx = plt.scatter(i.load, i.exx) # Plots each strain with a tag to create legend
    exy = plt.scatter(i.load, i.exy)
    eyy = plt.scatter(i.load, i.eyy)

    plt.grid(b=True, which='major', color='#666666', linestyle='-')

    plt.legend((exx, exy, eyy), ('ε_xx', 'ε_xy', 'ε_yy')) # Creates legend (matches respective terms)

    plt.title('Strain vs. Load '+ '(Torsion Arm Length ' + i.length + ')')
    plt.ylabel('Engineering Strain (microstrain)')
    plt.xlabel('Applied Load (N)')

### d ----------------------------------

plt.figure()
_7 = plt.scatter(i.load, L7cm.exy) # Plots each strain with a tag to create legend
_10 = plt.scatter(i.load, L10cm.exy)
_15 = plt.scatter(i.load, L15cm.exy)

plt.grid(b=True, which='major', color='#666666', linestyle='-')

plt.legend((_7, _10, _15), ('7 cm', '10 cm', '15 cm'))  # Creates legend (matches respective terms)

plt.title('ε_xy vs. Load for 3 Torsion Arm Lengths')
plt.ylabel('Engineering Strain (microstrain)')
plt.xlabel('Applied Load (N)')

### e ----------------------------------

for i in plt_array: # For each
    plt.figure()
    sxx = plt.scatter(i.load, i.sxx, marker='x', color = 'r') # Plots each stress with a tag to create legend
    txy = plt.scatter(i.load, i.txy, marker='o', color = 'r')
    syy = plt.scatter(i.load, i.syy, marker='+', color = 'r')

    theor_sxx = plt.scatter(i.load, i.theor_sxx, marker='x', color = 'b')
    theor_txy = plt.scatter(i.load, i.theor_txy, marker='o', color = 'b')
    theor_syy = plt.scatter(i.load, i.theor_syy, marker='+', color = 'b')

    plt.grid(b=True, which='major', color='#666666', linestyle='-')

    plt.legend((sxx, txy, syy, theor_sxx, theor_txy, theor_syy), ('σ_xx', 'σ_xy', 'σ_yy', 'Theoretical σ_xx', \
                    'Theoretical σ_xy', 'Theoretical σ_yy')) # Creates legend (matches respective terms)

    plt.title('Stress vs. Load '+ '(Torsion Arm Length ' + i.length + ')')
    plt.ylabel('Engineering Stress (MPa)')
    plt.xlabel('Applied Load (N)')

### f ----------------------------------

plt.figure()
_7 = plt.scatter(L7cm.load, L7cm.p1, marker='x', color = 'r')  # Plots each principal stress with a tag to create legend
_10 = plt.scatter(L10cm.load, L10cm.p1, marker='o', color = 'r')
_15 = plt.scatter(L15cm.load, L15cm.p1, marker='+', color = 'r')

_7_t = plt.scatter(L7cm.load, L7cm.p1_theor, marker='x', color = 'b')  # Plots each principal stress with a tag to create legend
_10_t = plt.scatter(L10cm.load, L10cm.p1_theor, marker='o', color = 'b')
_15_t = plt.scatter(L15cm.load, L15cm.p1_theor, marker='+', color = 'b')

plt.grid(b=True, which='major', color='#666666', linestyle='-')

plt.legend((_7, _10, _15, _7_t, _10_t, _15_t), ('7 cm', '10 cm', '15 cm', '7 cm theoretical', '10 cm theoretical',\
                                                '15 cm theoretical'))  # Creates legend (matches respective terms)

plt.title('Maximum Principal Stress vs. Load for 3 Torsion Arm Lengths')
plt.ylabel('Engineering Stress (MPa)')
plt.xlabel('Applied Load (N)')

### g ----------------------------------

plt.figure()
_7 = plt.scatter(L7cm.load, L7cm.theta_1, marker='x', color = 'r')  # Plots each principal stress with a tag to create legend
_10 = plt.scatter(L10cm.load, L10cm.theta_1, marker='o', color = 'r')
_15 = plt.scatter(L15cm.load, L15cm.theta_1, marker='+', color = 'r')

_7_t = plt.scatter(L7cm.load, L7cm.theta_1_theor, marker='x', color = 'b')  # Plots each principal stress with a tag to create legend
_10_t = plt.scatter(L10cm.load, L10cm.theta_1_theor, marker='o', color = 'b')
_15_t = plt.scatter(L15cm.load, L15cm.theta_1_theor, marker='+', color = 'b')

plt.grid(b=True, which='major', color='#666666', linestyle='-')

plt.legend((_7, _10, _15, _7_t, _10_t, _15_t), ('7 cm', '10 cm', '15 cm', '7 cm theoretical', '10 cm theoretical',\
                                                '15 cm theoretical'))  # Creates legend (matches respective terms)

plt.title('Principal Stress 1 Direction vs. Load for 3 Torsion Arm Lengths')
plt.ylabel('Angle (degrees)')
plt.xlabel('Applied Load (N)')

plt.show()
