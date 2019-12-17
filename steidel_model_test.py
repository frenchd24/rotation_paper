import numpy as np



def steidel_model(vc, inclination, hv, impact, azimuth, Dlos):
    """
        Calculate v_los based on the model of Steidel et al. (2002). Also used by
        Kacprzak et al. (2019)
        
        Parameters
        ----------
        vmax     : float
                maximum rotation velocity (inclination corrected)
                
        hv       : float
                scale height of the disk
                
        impact   : float
                impact parameter
                
        azimuth  : float
                azimuth angle
                
        y        : float
                the y component
        
        Returns
        -------
        vlos    : float
                the vlos at Dlos position along the sightline
        
    """
    
    # --- convert to radians
    to_radians = np.pi/180.
    az_radians = azimuth * to_radians
    inc_radians = inclination * to_radians
    
    y0 = impact * np.sin(az_radians) / np.cos(inc_radians)
    y = (Dlos * np.sin(inc_radians)) + y0
    p = impact * np.cos(az_radians)
    
    exp_term = np.exp(-(abs(y - y0)/(hv * np.tan(inc_radians))))
    
    sqrt_term = np.sqrt(1 + (y/p)**2)
    
    vlos = -(vc  * np.sin(inc_radians) / sqrt_term) * exp_term
    
    return vlos
    

vc = 232.286425223
inclination = 85.
hv = 0.1
impact = 26.109
azimuth = 0
Dlos = 0

print
print 'steidel_model: ',steidel_model(vc, inclination, hv, impact, azimuth, Dlos)




