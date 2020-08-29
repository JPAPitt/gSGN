"""
Python code to generate csv files that contain integrals of the basis functions
needed to generate the matrices of the Finite Element Method.

by Jordan Pitt 11/10/2018
"""

"""
############################# IMPORTS #########################################
"""


from sympy import * #Related website: https://www.sympy.org/en/index.html
from IPython.display import display


# symbols
a1, a2, b1, b2, c1, c2 = symbols('a_1 a_2 b_1 b_2 c_1 c_2',real=True)

A = Matrix([[a1,a2]])
B = Matrix([[b1,b2]])
C = Matrix([[c1,c2]])

Int = (A.dot(B)) * (A.dot(C))
Int = Int.expand()






