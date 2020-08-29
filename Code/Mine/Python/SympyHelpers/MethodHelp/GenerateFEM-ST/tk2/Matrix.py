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
    
    return exp


def PrintElements(Matrix,MatName,CoeffTerm):
    for i in range(Matrix.shape[0]):
        for j in range(Matrix.shape[1]):
            term = MatName + '['+str(i) + ',' + str(j)+']' + " = " + CoeffTerm + '('+ str(Matrix[i,j]) + ')'
            print(term)



#
# ---------------- Defintions ----------------------------------
#
x,dx = symbols('x dx')

# special locations in the cell in the \xi-space
# xjm1o2 = x_{j - 1/2}
xjm1o2 = "-1"
xjm1o6 = "-1/3"
xj = "0"
xjp1o6 = "1/3"
xjp1o2 = "1"

hjmh,hjms,hjps,hjph = symbols('hjmh,hjms,hjps,hjph')
ujmh,ujms,ujps,ujph = symbols('ujmh,ujms,ujps,ujph')
Gjmh,Gjms,Gjps,Gjph = symbols('Gjmh,Gjms,Gjps,Gjph')

Fjmh = Function('wjmh')
Fjms = Function('wjms')
Fjps = Function('wjps')
Fjph = Function('wjph')


wjmh = Fjmh(x)
wjms = Fjms(x)
wjps = Fjps(x)
wjph = Fjph(x)

#wjmh,wjms,wjps,wjph = symbols('wjmh,wjms,wjps,wjph')

wjmh_v = Cubic1(x,xjm1o2,xjm1o6,xjp1o6,xjp1o2,xjm1o2)
wjms_v = Cubic1(x,xjm1o2,xjm1o6,xjp1o6,xjp1o2,xjm1o6)
wjps_v = Cubic1(x,xjm1o2,xjm1o6,xjp1o6,xjp1o2,xjp1o6)
wjph_v = Cubic1(x,xjm1o2,xjm1o6,xjp1o6,xjp1o2,xjp1o2)

dwjmh =  diff(wjmh,x)
dwjms =  diff(wjms,x)
dwjps =  diff(wjps,x)
dwjph =  diff(wjph,x)

H = Matrix([[hjmh],[hjms],[hjps],[hjph]])
U = Matrix([[ujmh],[ujms],[ujps],[ujph]])
G = Matrix([[Gjmh],[Gjms],[Gjps],[Gjph]])
W = Matrix([[wjmh],[wjms],[wjps],[wjph]])
DW = diff(W, x) #Matrix([[dwjm1o2],[dwjm1o6],[dwjp1o6],[dwjp1o2]])

#G term
# # (G . W)  W = (W^T G)W = W(W^T G) = WW^T G
# GM = W * W.T
# GMint = integrate(GM.replace(wjmh,wjmh_v ).replace(wjms,wjms_v).replace(wjps,wjps_v ).replace(wjph,wjph_v )  , (x,-1,1))

# PrintElements(GMint,'Ge','dx/2*')


#UH term
# # (H . W) (U . W)  W =  (H^T W) (W^T U) W = (H^T W) W (W^T U) = = (H^T W W W^T) U
# UHM = (H.T * W)[0] *(W * W.T)
# UHMint = integrate(UHM.replace(wjmh,wjmh_v ).replace(wjms,wjms_v).replace(wjps,wjps_v ).replace(wjph,wjph_v )  , (x,-1,1))

# PrintElements(UHMint,'uhe','dx/2*')

# h^3/3 ux wx term
# # (H . W)  (H . W) (H . W) (U . DW)  DW =  (H . W)  (H . W) (H . W) (DW DW^T) U
H3DUM = (H.T * W)[0]*(H.T * W)[0]*(H.T * W)[0] *(DW * DW.T)
H3DUMRep = H3DUM.replace(wjmh,wjmh_v ).replace(wjms,wjms_v).replace(wjps,wjps_v ).replace(wjph,wjph_v ).doit()
H3DUMint = integrate(H3DUMRep, (x,-1,1))

PrintElements(H3DUMint,'h3uxe','(1.0/3.0)*2/dx*')