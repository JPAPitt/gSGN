"""
Python code to generate csv files that contain integrals of the basis functions
needed to generate the matrices of the Finite Element Method.

by Jordan Pitt 11/10/2018
"""

"""
############################# IMPORTS #########################################
"""


from sympy import * #Related website: https://www.sympy.org/en/index.html
import csv


"""
########################## Function Definitions ###############################
"""

#
# ---------------- Basis Function Generators ----------------------------------
#

def lineardiscP(x,x0,x1):
    # Generates the basis functions \psi_{j+1/2}^- from the thesis.
    # x is just a variable and x0 and x1 are the left and right cell edge 
    # respectively. Returns a tuple: (basis function, left edge, right edge)
    return ((x - S(x1)) / (S(x0) - S(x1)) ,S(x0), S(x1))

def lineardiscM(x,x0,x1):
    # Generates the basis functions \psi_{j+1/2}^+ from the thesis.
    # x is just a variable and x0 and x1 are the left and right cell edge 
    # respectively. Returns a tuple: (basis function, left edge, right edge)
    return ((x - S(x0)) / (S(x1) - S(x0)) ,S(x0), S(x1))

        
def quadcontEM(x,x0,x1,x2):    
    # Generates the basis functions \phi_{j-1/2} from the thesis.
    # x is just a variable and x0 and x2 are the left and right cell edge 
    # with x1 the cell midpoint x_{j} respectively.
    # Returns a tuple: (basis function, left edge, right edge)
    sexp = (x - S(x1))*(x - S(x2)) / ((S(x0) - S(x1))*(S(x0) - S(x2)))
    return (sexp,S(x0),S(x2))

def quadcontEP(x,x0,x1,x2):    
    # Generates the basis functions \phi_{j+1/2} from the thesis.
    # x is just a variable and x0 and x2 are the left and right cell edge 
    # with x1 the cell midpoint x_{j} respectively.
    # Returns a tuple: (basis function, left edge, right edge)
    sexp = (x - S(x1))*(x - S(x2)) / ((S(x0) - S(x1))*(S(x0) - S(x2)))
    return (sexp,S(x2),S(x0))


def quadcontM(x,xm1,x0,x1):
    # Generates the basis functions \phi_{j+1/2} from the thesis.
    # x is just a variable and x0 and x2 are the left and right cell edge 
    # with x1 the cell midpoint x_{j} respectively.
    # Returns a tuple: (basis function, left edge, right edge)
    fexp = (x - S(xm1))*(x - S(x1)) / ((S(x0) - S(xm1))*(S(x0) - S(x1)))
    return (fexp,S(xm1),S(x1))

def Cubic1(x,x0,x1,x2,x3,loc1): 
    
    # Generates the basis functions \gamma_{j-1/2},\gamma_{j-1/6},
    # \gamma_{j+1/6},\gamma_{j+1/2} from the thesis.
    # x is just a variable and x0 and x3 are the left and right cell edge 
    # while x1 and x2 are x_{j-1/2} and x_{j+1/2} respectively.
    # loc1 determines which of the \gamma basis functions are reconstructed
    # by determining where the basis function is 1 so that
    # loc1 = x0 => reconstruct \gamma_{j-1/2},
    # loc2 = x1 => reconstruct \gamma_{j-1/6},etc.
    # Returns a tuple: (basis function, left edge, right edge)
    if(loc1 == x0):
        exp = (x - S(x1))*(x - S(x2))*(x - S(x3))  / ((S(x0) - S(x1))*(S(x0) - S(x2))*(S(x0) - S(x3)))
    elif(loc1 == x1):
        exp = (x - S(x0))*(x - S(x2))*(x - S(x3))  / ((S(x1) - S(x0))*(S(x1) - S(x2))*(S(x1) - S(x3)))
    elif(loc1 == x2):
        exp = (x - S(x1))*(x - S(x0))*(x - S(x3))  / ((S(x2) - S(x1))*(S(x2) - S(x0))*(S(x2) - S(x3)))
    else:
        exp = (x - S(x1))*(x - S(x2))*(x - S(x0))  / ((S(x3) - S(x1))*(S(x3) - S(x2))*(S(x3) - S(x0)))
    
    return (exp,S(x0),S(x3))


def phideriv(x,phi):
    # Given a basis function tuple: (basis function, left edge, right edge)
    # returns its derivative
        return (diff(phi[0],x),phi[1],phi[2])        




#
# ---------------- CSV File Generator ----------------------------------
#
def ExpressionConvert(s,hyperlist,x,Cterm,Ctermname,Clb,Cub,filen):
    #This is a recursive function that integrates all the basis function
    #combinations and then writes them to a csv file.
    
    # s is a string representing the basis function combination to integrate
    
    # hyperlist is a list asociating the 'key' charachters of the string s
    #   to their basis function class
    
    # x is the symbolic variable 'x'
    
    # Cterm is the current term; consisting of multiplications of the various
    #   basis functions
    
    # Clb and Cub are the lower and upper bounds on the integration domain
    
    # filen is the name of the csv file to write to
    
    
    # hyperlist has is a list of tuples of (char,functionlist,functionname)
    
    n = len(hyperlist)
    
    # If length of s is zero, then function is complete and all basis functions
    # have been multiplied together, this is the stopping condition of the recursive function
    if (len(s) == 0):
        # Integrate then write to the csv file
        Intv = integrate(Cterm,(x,Clb,Cub))
        filen.writerow([str(Cterm),str(Clb),str(Cub), str(Intv)] + Ctermname)        

    # If Cterm is none then this basis function is the first in our experession
    # this is the initialisation of the recursive function
    elif(Cterm is None):   
        
        # Find what basis functions the last charachter in 's' corresponds to
        for i in range(n):
            if(s[-1] == hyperlist[i][0]):
                m = len(hyperlist[i][1])
                
                #For all these basis functions Cterm becomes
                # the asscociated basis fuction and we recurse
                for j in range(m):
    
                    ExpressionConvert(s[:-1],hyperlist,x,hyperlist[i][1][j][0],[hyperlist[i][2][j]],Clb,Cub,filen)
    
    
    # Since it is neither the initial step or the stopping conditions, we recruse down the string
    else:
        
        # Find what basis functions the last charachter in 's' corresponds to
        for i in range(n):
            if(s[-1] == hyperlist[i][0]):
                m = len(hyperlist[i][1])
                
                #For all these basis functions Cterm becomes
                # the asscociated basis fuction and we recurse
                for j in range(m):
                    ExpressionConvert(s[:-1],hyperlist,x,Cterm*hyperlist[i][1][j][0],[hyperlist[i][2][j]] + Ctermname,Clb,Cub,filen)
                    

"""
########################## Code Body ###############################
"""

#
# ---------------- Defintions ----------------------------------
#
x = Symbol('x')

# special locations in the cell in the \xi-space
# xjm1o2 = x_{j - 1/2}
xjm1o2 = "-1"
xjm1o6 = "-1/3"
xj = "0"
xjp1o6 = "1/3"
xjp1o2 = "1"


# Define all the basis function tuples we need
# phijm1o2 = \phi_{j-1/2} in the thesis
# dphijm1o2 =  derivative of \phi_{j-1/2} in the thesis
# wijm1o2
wjm1o2p = Cubic1(x,xjm1o2,xjm1o6,xjp1o6,xjp1o2,xjm1o2)
wjm1o6 = Cubic1(x,xjm1o2,xjm1o6,xjp1o6,xjp1o2,xjm1o6)
wjp1o6 = Cubic1(x,xjm1o2,xjm1o6,xjp1o6,xjp1o2,xjp1o6)
wjp1o2m = Cubic1(x,xjm1o2,xjm1o6,xjp1o6,xjp1o2,xjp1o2)

dwjm1o2p =  phideriv(x,wjm1o2p)
dwjm1o6 =  phideriv(x,wjm1o6)
dwjp1o6 =  phideriv(x,wjp1o6)
dwjp1o2m =  phideriv(x,wjp1o2m)

# Names for the basis functions ot be used in CSV file
wsNAME = ["wjm1o2","wjm1o6","wjp1o6","wjp1o2"]
dwsNAME = ["dwjm1o2","dwjm1o6","dwjp1o6","dwjp1o2"]

# #List of the basis functions
ws = [wjm1o2p,wjm1o6,wjp1o6,wjp1o2m]
dws = [dwjm1o2p,dwjm1o6,dwjp1o6,dwjp1o2m]

#hyperlist assosciates the charachter to the basis function type
hyperlist = []

hyperlist.append( ('w',ws,wsNAME))
hyperlist.append( ('d',dws,dwsNAME))

# #
# # ---------------- CSV File Generation Runs ----------------------------------
# #

# # For integral of Gt need to integrate all combinations of psis and phis that
# # are non-zero over the cell, hence s = 'ph', (t is the test function)
s = "Ge.csv"
with open(s,'w') as file2:
      writefile2 = csv.writer(file2, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)

      writefile2.writerow(['Integrand','lower bound of integral','upper bound of integral', 'Value of Integral', 'W Term 1', 'P Term 1'])        
      ExpressionConvert('ww',hyperlist,x,None,None,-1,1,writefile2)      


s = "uhe.csv"
with open(s,'w') as file2:
      writefile2 = csv.writer(file2, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)

      writefile2.writerow(['Integrand','lower bound of integral','upper bound of integral', 'Value of Integral', 'W Term 1' , 'P Term 1', 'P Term 2'])        
      ExpressionConvert('www',hyperlist,x,None,None,-1,1,writefile2)  

s = "h3uxe.csv"
with open(s,'w') as file2:
      writefile2 = csv.writer(file2, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)

      writefile2.writerow(['Integrand','lower bound of integral','upper bound of integral', 'Value of Integral', 'W Term 1', 'W Term 2' , 'W Term 3' , 'D Term 1', 'D Term 2'])        
      ExpressionConvert('wwwdd',hyperlist,x,None,None,-1,1,writefile2)  

