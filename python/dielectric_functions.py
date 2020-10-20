from numpy import pi
c_speed = 3e17 # nm

def lorentzAg(wavelength):
    """
    wavelength: Wavelength in nanometers

    Outputs dielectric constant (relative permittivity) for the given wavelength
    """
    freq = 2*pi*c_speed/wavelength
    eps_inf_Ag = 1
    f_Ag1 = 0.845
    f_Ag2 = 0.065
    f_Ag3 = 0.124
    f_Ag4 = 0.011
    f_Ag5 = 0.84
    f_Ag6 = 5.646

    wo_Ag1 = 0
    wo_Ag2 = 1.2397e+15
    wo_Ag3 = 6.8078e+15
    wo_Ag4 = 1.2435e+16
    wo_Ag5 = 1.38e+16
    wo_Ag6 = 3.0826e+16
    wp_Ag = 1.3689e+16

    y_Ag1 = 7.2925e+13
    y_Ag2 = 5.9039e+15
    y_Ag3 = 6.8671e+14
    y_Ag4 = 9.8752e+13
    y_Ag5 = 1.3916e+15
    y_Ag6 = 3.6751e+15

    return eps_inf_Ag + f_Ag1*wp_Ag**2/(wo_Ag1**2-freq**2-complex(0,1)*y_Ag1*freq) + \
    f_Ag2*wp_Ag**2/(wo_Ag2**2-freq**2-complex(0,1)*y_Ag2*freq) + \
    f_Ag3*wp_Ag**2/(wo_Ag3**2-freq**2-complex(0,1)*y_Ag3*freq) + \
    f_Ag4*wp_Ag**2/(wo_Ag4**2-freq**2-complex(0,1)*y_Ag4*freq) + \
    f_Ag5*wp_Ag**2/(wo_Ag5**2-freq**2-complex(0,1)*y_Ag5*freq) + \
    f_Ag6*wp_Ag**2/(wo_Ag6**2-freq**2-complex(0,1)*y_Ag6*freq)

def lorentzTDBC(wavelength):
    """
    wavelength: Wavelength in nanometers

    Outputs dielectric constant (relative permittivity) for the given wavelength
    """
    freq = 2*pi*c_speed/wavelength
    eps_b_organic3 = 2
    fo_organic3 = 0.2
    f1_organic3 = 0.03

    wo_organic3 = 3.236098450319052e+15
    w1_organic3 = 3.5706e15

    y1_organic3 = 3.038590094196293e+14
    yo_organic3 = 7.596475235490733e+13
    
    return eps_b_organic3 + \
           fo_organic3*wo_organic3**2/(wo_organic3**2 - freq**2 - complex(0,1)*yo_organic3*freq) + \
           f1_organic3*w1_organic3**2/(w1_organic3**2 - freq**2 - complex(0,1)*y1_organic3*freq)

def sellmeier(wavelength):
    """
    wavelength: Wavelength in nanometers 

    Outputs dielectric constant (relative permittivity) for the given wavelength
    """  
    B1 = 0.6961663
    B2 = 0.4079426
    B3 = 0.8974794
    C1 = 0.0684043
    C2 = 0.1162414
    C3 = 9.896161
    return 1+(B1*wavelength**2/(wavelength**2-C1**2))+\
           (B2*wavelength**2/(wavelength**2-C2**2))+\
           (B3*wavelength**2/(wavelength**2-C3**2))