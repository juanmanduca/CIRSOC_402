'''Module with the functions necessary to compute the bearing
capacity of piles
'''
import numpy as np

def tip():
    pass

def bolton(phi_cv,K0,Dr,qu,Q=10):
    p = qu*(1+ 2*K0)/3
    phi = phi_cv + 3*Dr*(Q - np.log(p)) - 3
    return phi


def young_modulus(Dr,K0,gamma_Nq,Depth, p_atm = 100):
    
    C = 100 + 1000*(Dr + Dr**2)
    m = 2 - 0.5*np.log(C)    
    sigma_eff_vert = gamma_Nq*Depth   
    E = p_atm*C*(sigma_eff_vert*K0/p_atm)**m 

    return E


def I_rr(phi,Dr,Depth,K0,gamma_Nq,cohesion,poisson):
    