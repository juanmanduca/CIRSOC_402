'''Module with the functions necessary to compute the bearing
capacity of piles
'''
import numpy as np



def bolton(phi_cv,K0,Dr,qu,Q=10):
    """Bolton equation use to calculate critic friction angle 
     according to reference [xxx]     
    Parameters
     ----------
    phi_cv : float,int
        Friction angle in a constant volume         
    K0 : float , int
        Coefficient of earth pressure at rest        
    Dr : float, int
        Relative density of soil [-]
    qu : float, int
        Ultimate capacity of the tip in piles, qu is calculated in the
        function of bearing capacity [kN/m2]
    Q : float, int
        parameter of the type of grain [-]   
    Returns
    -------
    float
        Critic friction angle [°]

    """
    p = qu*(1+ 2*K0)/3
    phi = phi_cv + 3*Dr*(Q - np.log(p)) - 3
    return phi
  

def young_modulus(Dr,K0,gamma_Nq,Depth, p_atm = 100):
    """Use to calculate the young modulus of the soil using the relativ density 
     according to reference [xxx]     
    Parameters
     ----------
    Dr : float, int
        Relative density of soil [-]       
    K0 : float , int
        Coefficient of earth pressure at rest        
    gamma_nq : float, int
        Unit weight used in the surcharge term in the bearing capacity
        equation [kN/m3]
    Depth: float, int
        Depth of the foundation level [m]
    p_atm: float, int
        Value of the atmosferic pressure [KN/m2]   
    Returns
    -------
    float
        Young modulus [KN/m2]
    """    
    C = 100 + 1000*(Dr + Dr**2)
    m = 2 - 0.5*np.log(C)    
    sigma_eff_vert = gamma_Nq*Depth   
    E = p_atm*C*(sigma_eff_vert*K0/p_atm)**m 

    return E


def Irr(phi,Dr,Depth,K0,gamma_Nq,cohesion,poisson):
    """Use to calculate the Irr of the vesic formule  
     according to reference [xxx]     
    Parameters
     ----------
    phi : float, int
        Critic friction angle [°] 
    Dr : float, int
        Relative density of soil [-]  
    Depth: float, int
        Depth of the foundation level [m]     
    K0 : float , int
        Coefficient of earth pressure at rest        
    gamma_nq : float, int
        Unit weight used in the surcharge term in the bearing capacity
        equation [kN/m3]
    cohesion : float, int
        Soil cohesion [kPa]
    poisson: float, int
        Poisson's ratio  [-]

      Returns
    -------
    float
        Irr [KN/m2]
    """       
    phi_rad = np.radians(phi)
    TanPhi = np.tan(phi_rad)
    sigma_eff_vert = gamma_Nq*Depth    
    I_r = young_modulus(Dr,K0,gamma_Nq,Depth , p_atm=100) / ((cohesion + sigma_eff_vert *TanPhi)*(1 + poisson)*2)    
    Delta = 50*(I_r)**(-1.8)    
    I_rr = I_r / (1 + I_r*Delta)
    
    return I_rr

def N_sigma(phi,Depth,Dr,K0, gamma_Nq, cohesion,poisson):
    """Use to calculate the young modulus of the soil using the relativ density 
     according to reference [xxx]     
    Parameters
     ----------
    phi : float, int
        Critic friction angle [°]  
    Depth: float, int
        Depth of the foundation level [m]   
    Dr : float, int
        Relative density of soil [-]       
    K0 : float , int
        Coefficient of earth pressure at rest        
    gamma_nq : float, int
        Unit weight used in the surcharge term in the bearing capacity
        equation [kN/m3]
    cohesion : float, int
        Soil cohesion [kPa]
    poisson: float, int
        Poisson's ratio  [-]   
    Returns
    -------
    float
        N_sigma [-]
    """    
    phi_rad = np.radians(phi)
    SinPhi = np.sin(phi_rad)
    TanPhi = np.tan(phi_rad)    
    c1 = (3 / (3 - SinPhi)) * np.exp((np.pi/2 - phi_rad)*TanPhi)
    c2 = np.tan(np.pi/4 + phi_rad/2)**2 
    c3 = (4*SinPhi) / (3*(1 + SinPhi))   
    I_rr = Irr(phi,Dr,Depth,K0,gamma_Nq,cohesion,poisson) 
    N_sigma = c1*c2*I_rr**c3

    return N_sigma

def Nc_vesic(phi,Depth,Dr,gamma_Nq,cohesion,poisson,K0):
    """Parameter use to calculate the bearing capacity by Vesic [xxx]     
    Parameters
     ----------
    phi : float, int
        Critic friction angle [°]  
    Depth: float, int
        Depth of the foundation level [m]    
    Dr : float, int
        Relative density of soil [-]    
    gamma_nq : float, int
        Unit weight used in the surcharge term in the bearing capacity
        equation [kN/m3] 
    cohesion : float, int
        Soil cohesion [kPa]       
    poisson: float, int
        Poisson's ratio  [-]
    K0 : float , int
        Coefficient of earth pressure at rest 
    Returns
    -------
    float
        N_c [-]
    """    
    phi_rad = np.radians(phi)    
    TanPhi = np.tan(phi_rad)    
    sigma_eff_vert =gamma_Nq*Depth    
    sigma_0 = ((1 + 2*K0)/3 )* sigma_eff_vert    
    N_q = (sigma_0*N_sigma(phi,Depth,Dr,K0,gamma_Nq,cohesion,poisson)) / sigma_eff_vert
    I_rr = Irr(phi,Dr,Depth,K0,gamma_Nq,cohesion,poisson)     
    if phi == 0:         
        N_c = (4/3) * (np.log(I_rr) + 1) + np.pi/2 + 1

    elif phi > 0: 
        N_c = (N_q - 1) / TanPhi
        
    return N_c    

def bearingcapacity(poisson,K0,Depth,Dr,shape,foundation, method, gamma_Ng, gamma_Nq,phi, phi_cv, cohesion, effective_width,  effective_length,surface_load ):   
    """function which estimates the ultimate bearing capacity of the pile tip
     as the minimum between the method of vesic and Brinch Hanse [xxx]     
    Parameters
     ----------
    poisson: float, int
        Poisson's ratio  [-]
    K0 : float , int
        Coefficient of earth pressure at rest
    Depth: float, int
        Depth of the foundation level [m] 
    Dr : float, int
        Relative density of soil [-] 
    shape : str
        Shape of the foundation. The supported shapes can be seen with
        cirsoc_402.constants.BEARINGSHAPE.
    foundation : str
        type of the foundation. The supported types can be seen with
        cirsoc_402.constants.BEARINGSHAPE.   
    method : str
        Calculation method for the soil weight bearing capacity factor.
        The supported methods can be seen with
        cirsoc_402.constants.BEARINGMETHOD. 
    gamma_ng : float, int
        Unit weight used in the soil weight term in the bearing capacity
        equation [kN/m3]
    gamma_nq : float, int
        Unit weight used in the surcharge term in the bearing capacity
        equation [kN/m3]
    phi : float, int
        Critic friction angle [°]  
    phi_cv : float,int
        Friction angle in a constant volume
    cohesion : float, int
        Soil cohesion [kPa]
    poisson: float, int
        Poisson's ratio  [-]
    effective_width : float, int
        Effective width of the equivalent rectangular load area [m] 
    effective_length : float, int
        Effective length of the equivalent rectangular load area [m]
    surface_load : float, int, optional
        Load acting on the ground surface considered in the surcharge
        term of the bearing capacity equation. by default 0 [kPa]              
    Returns
    -------
    float
        qu [KN/m2]
    """    
    qu_c = cohesion * Nc(phi) * c_cs(shape, phi, B=B, L=L) * c_cd(Depth, B, phi)
    qu_p = (q + Depth * gamma_Nq) * Nq(phi) * c_qs(shape, phi, B=B, L=L) * c_qd(Depth, B, phi)
    sigma_0 = ((1 + 2*K0)/3 )* (q + Depth*gamma_Nq)

    if foundation == 'shallow':
        if shape == ['rectangulo','rectangulo']:
            qu_g_B = 0.5 * gamma_Ng * B * Ng(phi, method) * c_gs(shape, B=B, L=L) * c_gd(phi)
            qu_g_L = 0.5 * gamma_Ng * L * Ng(phi, method) * c_gs(shape, B=B, L=L) * c_gd(phi)
            qu = np.minimum(qu_c + qu_p + qu_g_B , qu_c + qu_p + qu_g_L)
        else:
            qu = qu_c + qu_p + 0.5 * gamma_Ng * B * Ng(phi, method) * c_gs(shape, B=B, L=L) * c_gd(phi)
        
    elif foundation == 'pile':
        
        qu_BH = 0.5 * gamma_Ng * B * Ng(phi, method) * c_gs(shape, B=B, L=L) * c_gd(phi)
        qu_vesic = N_sigma(phi,Depth,Dr,K0, gamma_Nq, cohesion,poisson)*sigma_0 + Nc_vesic(phi,Depth,Dr,gamma_Nq,cohesion,poisson,K0)*cohesion
        qu = np.minimum(qu_BH , qu_vesic)        
    phi = bolton(phi_cv,K0,Dr,qu,Q=10)       
                
    return(phi,qu)

def tip(poisson,K0,D,Dr,shape,foundation, method, gamma_Ng, gamma_Nq,phi , phi_cv, cohesion, effective_width,  effective_length,surface_load ):
    """function which iterates the values of phi in bolton equation using the ultimate load of the piles
    and returns the ultimate load of the tip
     ----------
    poisson: float, int
        Poisson's ratio  [-]
    K0 : float , int
        Coefficient of earth pressure at rest
    Depth: float, int
        Depth of the foundation level [m] 
    Dr : float, int
        Relative density of soil [-] 
    shape : str
        Shape of the foundation. The supported shapes can be seen with
        cirsoc_402.constants.BEARINGSHAPE.
    foundation : str
        type of the foundation. The supported types can be seen with
        cirsoc_402.constants.BEARINGSHAPE.   
    method : str
        Calculation method for the soil weight bearing capacity factor.
        The supported methods can be seen with
        cirsoc_402.constants.BEARINGMETHOD. 
    gamma_ng : float, int
        Unit weight used in the soil weight term in the bearing capacity
        equation [kN/m3]
    gamma_nq : float, int
        Unit weight used in the surcharge term in the bearing capacity
        equation [kN/m3]
    phi : float, int
        Critic friction angle [°]  
    phi_cv : float,int
        Friction angle in a constant volume
    cohesion : float, int
        Soil cohesion [kPa]
    poisson: float, int
        Poisson's ratio  [-]
    effective_width : float, int
        Effective width of the equivalent rectangular load area [m] 
    effective_length : float, int
        Effective length of the equivalent rectangular load area [m]
    surface_load : float, int, optional
        Load acting on the ground surface considered in the surcharge
        term of the bearing capacity equation. by default 0 [kPa]              
    Returns
    -------
    float
        Qu [KN/m2]
    """
    # Definición de lista para valores de  numero de  iteraciones 
    Cont =[]
    # Definición de lista para valores iterativos de phi
    Phi = []
    # Definición de lista para valores iterativos de qu
    Qu = []    
    contador = 15   
    for i in range(contador):
        Phi_i = bearingcapacity(poisson,K0,D,Dr,shape,foundation, method, gamma_Ng, gamma_Nq,phi, phi_cv, c, B, q = 0, L =np.nan )[0]
        Qu_i = bearingcapacity(poisson,K0,D,Dr,shape,foundation, method, gamma_Ng, gamma_Nq,phi, phi_cv, c, B, q = 0, L =np.nan )[1]
        Phi.append(Phi_i)
        Qu.append(Qu_i)
        Cont.append(i + 1)
        phi = Phi_i
       
    return  phi[15] , Qu[15]
        




def shaft(phi,dilatancy,young_modulus,diameter,gamma,depth,v,poisson=0.2,K0=0.5):
    
    phi_rad = np.radians(phi)
    psi_rad = np.radians(dilatancy)
    
    shear_modulus = young_modulus/(2*(1 + poisson))
    
    Delta_sigma_h = (shear_modulus/(0.5*diameter))*np.tan(psi_rad)*v

    sigma_h = K0*gamma*depth
    
    fs = sigma_h + Delta_sigma_h
    