import constants as costs
import code_units as cdu
import numpy as np

class uniform:
    """
        This class simply give the value of pressure, internal energy, density and velocity of
        a uniform distribution in code units. Code units are set in code_units.py
    """
    def __init__(self, n = 1., T = 1e4, gamma = 5./3 ) -> None:
        """

        Args:
            n     (float, cm**-3)  : Particle dentity of the gas. Defaults to 1.
            T     (float, K)       : Temperature of the gas. Defaults to 1e4.
            gamma (float, dimensionsless): Gas constant. Defaults to 5./3.
        """
        
        self.gamma = gamma
        self.costs = costs.kB / (costs.mu_mol * costs.mH * (self.gamma - 1))
        self.n = n 
        self.T = T

    def d(self) -> float :
        """
        Density.

        Returns:
            float: return the density of the gas in code units.
        """
        rho_cgs = self.n * costs.mu_mol * costs.mH
        
        return rho_cgs / cdu.scale_d
    
    def e(self) -> float: 
        """
        Energy

        Returns:
            float: return the internal energy of the gas in code units [velocity**2]
        """
        e_cgs = self.T * self.costs
        return e_cgs / cdu.scale_v**2
    
    def p(self):
        """
        Pressure.

        Returns:
            float: pressure of the gas in code units.
        """
        return self.d() * self.e() * (self.gamma - 1)
    
    def v(self,v):
        """the
        Velocity.

        Args:
            v (float, np.array): velocity of the gasin km/s.

        Returns:
            float: velocity of the gas in code units.
        """
        return v * cdu.kms / cdu.scale_v



class exponential:
    """
    This class generates the potential and density profiles associated to an exponential distribution of gas.
    Common ARGS:
        z       (float/ 1-2-3D array)        : z-coordinates.
        y       (float/ 2-3Darray , optional): y-coordinates. Defaults to 0. When this is different from zero the result should become 2D.
                                                                        If x is an array, that should turn 3D too. 
        x      (float/ 3Darray   , optional): x-coordinates. Defaults to 0. When this is different from a float the result should become 3D
        
    """
    def __init__(self, 
                midplane_density = 0.7,
                midplane_temp    = 8e3,
                B0               = 0. ,
                H                = 150.,
                R0               = 8.e0 * costs.kpc,
                vscale           = 1.,
                ) -> None:
        """
        Outputs depends on the units used to make the computations.
        Args:
            midplane_density (float, optional): density of the gas in the midplane. Defaults to 0.7 cm**-3.
            midplane_temp    (float, optional): density of the gas in the midplane. Defaults to 8e3 K.          
            B0               (float, optional): Magnetic field magnitude. Defaults to 0..
            H                (float, optional): scale height of the disk. Defaults to 150.
            vscale           (floar, optional): velocity units. Needed to convert relevant constants 
                                                in case the parameters are not in cgs. Defaults to 1.
            R0               (float, optional): scale parameter for the gas distribution. Not used. Defaults to 8.e0*kpc.
        """
        self.d_mp             = midplane_density * costs.mu_mol * costs.mH 
        self.T_mp             = midplane_temp
        self.H                = H        
        self.B0               = B0 
        self.R0               = R0
        # computing useful constants
        self.Rgas          = costs.kB * self.T_mp  / (costs.mu_mol * costs.mH) / vscale
        self.RMgas         = B0**2 / ( 8. * costs.pi * self.d_mp ) 
        self.g             = self.Rgas / self.H
        
        pass
    def phi(self, z):
        """
        External gravitational potent(z)ial for the gas  distribution.

        Args:
            z (float/array): coordinates

        Returns:
            float/array: potential at the corresponding coordinates.
        """
        return self.g*z
    
    def profile(self, z, z_mid = 0. ):
        """
        Isothermal density distribution associated to the potential.

        Args:
            z      (float/array): coordinates
            z_mid  (float)      : where to compute phi_0. 

        Returns:
            _type_: _description_
        """
        delta_phi = -(self.phi(z) - self.phi(z_mid)) 
        return self.d_mp * np.exp(delta_phi / (self.Rgas + self.RMgas))
    


class BetaModel:
    """
        Beta model. This is supposed to be a spherical model.  Only implemented for beta = 2/3. 
    """
    def __init__(self,  
                 beta = 2./3.,
                 Rc = 8 * costs.kpc2cm,
                 rho_c = 1):
        """

        Args:
            beta  (float, optional) : Beta parameter. Defaults to 2./3..
            Rc    (_type_, optional): Radius of the distribution core. Defaults to 8*costs.kpc2cm.
            rho_c (int, optional)   : Core density of the model. sDefaults to 1.
        """
        self.beta  = beta
        self.Rc    = Rc
        self.rho_c = rho_c
        pass

    def profile(self, z, y = 0., x = 0.):
        """
        Generate the density profile.
        
        Returns:
            float/array: density distribution of the beta model
        """

        
        r = np.sqrt( x** 2 + y** 2 + z**2 )
    
        return self.rho_c * (1 + ( r / self.Rc )**2.) ** (-3./2. * self.beta)
    
    def profile_iso(self, 
                z, y = 0., x = 0.,  
                midplane_temp    = 8e3,
                B0               = 0. ,
                vscale           = 1.):
        """
        Isothermal profile corresponding to the potential of a beta-model. Multidimensionality need to be tested. 

        Args:
            midplane_temp (float, optional): temperature of the gas. Defaults to 8e3K.
            B0            (float, optional): magnetic field magnitude. Defaults to 0.
            vscale        (float, optional): velocity units. Needed to convert relevant constants in case the parameters are not in cgs. Defaults to 1.
            
        Returns:
            float/array: density distribution of the isothermal beta model
        """
        d_mp          = self.rho_c
        Rgas          = costs.kB * midplane_temp/(costs.mu_mol * costs.mH) / vscale
        RMgas         = B0**2 / ( 8. * costs.pi * d_mp ) 
        
        r = np.sqrt( x** 2 + y**2 + z**2 )
        delta_phi     = (self.phi(0.) - self.phi(r)) 
        
        return d_mp * np.exp(delta_phi / (Rgas + RMgas))
    
    # potential
    def phi(self, z, y = 0., x = 0.):
        """
            Potential coming from the integration of the poisson equation... not sure it makes sense 
        """
        if self.beta != 2./3.: raise Exception("Only beta allowed for phi is 2/3")
        
        r = np.sqrt( x** 2 + y**2 + z**2 )
        
        return - 2. * np.pi * costs.factG_in_cgs * self.rho_c * self.Rc *(
            self.Rc * np.log( self.Rc**2 + r**2 ) +
            2 * z * np.arctan2( self.Rc, r)  
            )


    def dphi(self, z, y = 0., x = 0.):
        """Gravitational acceleration associated to self.phi.

        Args:

        Raises:
            Exception: if beta != 2/3

        Returns:
            float/array: gravitational force
        """
        if self.beta != 2./3.: raise Exception("Only beta allowed for phi is 2/3")
        
        r = np.sqrt( x** 2 + y**2 + z**2 )
        
        return -4. * np.pi * costs.factG_in_cgs * self.rho_c * self.Rc * np.arctan(self.Rc / r)
    
    def drho(self, z, y = 0., x = 0.):
        """
        Density gradient.

        Returns:
            float/array: density gradient
        """
        r = np.sqrt( x** 2 + y**2 + z**2 )
        return - (2. * r * self.rho_c * self.Rc**2) / ((self.Rc**2 + r**2)**2)

class Barros:
    """This class computes the potentials associated to the various disks described in Barros et al. 2015. 
        It also allows to compute the density profiles associated to it. 3D not implemented.
        Common Args:
            z (float/array): z-coordinates.

    """
    def __init__(self, 
                
                midplane_density = 0.7,
                midplane_temp    = 8e3,
                midplane_B0      = 0.,
                R0               = 8.e0 * costs.kpc,
                ) -> None:
        """

        Args:
            
            midplane_density(float, optional): Density of the midplane.     Defaults to 0.7.
            midplane_temp   (float, optional): Temperature of the midplane. Defaults to 8e3.
            midplane_B0     (float, optional): Magnetic field magnitude at the disk midplane. Defaults to 0..
            R0              (float, optional): Radial distance from the center of the galaxy. Defaults to 8.e0*costs.kpc.
        """
        
        self.d_mp = midplane_density * costs.mu_mol * costs.mH 

        self.B_mp          = midplane_B0 
        self.T_mp          = midplane_temp
        self.R0            = R0
        self.Rgas          = costs.kB * self.T_mp / (costs.mu_mol * costs.mH)
        self.RMgas         = self.B_mp**2 / (8. * costs.pi * self.d_mp )


        pass

    
        

    def profile(self, z):
        """
        Isothermal distribution of the gas.

        Returns:
            float/array: density distribution.
        """

        phi_tot_0 = self.phi( 0.)#phi_d_0 + self.phi_b(0.e0,2.61e10*M_sun,0.44e0*kpc2cm) + self.phi_h(0.e0,5.4e0*kpc2cm,166.e0*kms)  
        phi_tot = self.phi( z )#phi_d_0 + self.phi_b(0.e0,2.61e10*M_sun,0.44e0*kpc2cm) + self.phi_h(0.e0,5.4e0*kpc2cm,166.e0*kms)  
        
        d    = self.d_mp*np.exp((phi_tot_0-phi_tot)/((self.Rgas +self.RMgas  )))#*scale_v**2))
        #print(phi_tot, phi_tot_0)
        return d
    def phi(self, z):
        """Compute total potential of the disk at z.


        Returns:
            float/array: potentials
        """
        zz = z
        phi_d =  self.phi_d_thin( zz, 2.106e10*costs.M_sun, 3.859e0*costs.kpc2cm, 2.162e10*costs.M_sun, 9.052e0*costs.kpc2cm,-1.704e10*costs.M_sun, 3.107e0*costs.kpc2cm, 0.243e0*costs.kpc2cm)
        phi_d += self.phi_d_thick(zz, 0.056e10*costs.M_sun, 0.993e0*costs.kpc2cm, 3.766e10*costs.M_sun, 6.555e0*costs.kpc2cm,-3.250e10*costs.M_sun, 7.651e0*costs.kpc2cm, 0.776e0*costs.kpc2cm)
        phi_d += self.phi_d_HI(   zz, 2.046e10*costs.M_sun, 9.021e0*costs.kpc2cm, 2.169e10*costs.M_sun, 9.143e0*costs.kpc2cm,-3.049e10*costs.M_sun, 7.758e0*costs.kpc2cm, 0.168e0*costs.kpc2cm)
        phi_d += self.phi_d_H2(   zz, 0.928e10*costs.M_sun, 6.062e0*costs.kpc2cm, 0.163e10*costs.M_sun, 3.141e0*costs.kpc2cm,-0.837e10*costs.M_sun, 4.485e0*costs.kpc2cm, 0.128e0*costs.kpc2cm)
        phi_tot = phi_d+(self.phi_b(zz,2.61e10*costs.M_sun,0.44e0*costs.kpc2cm) +
                         self.phi_h(zz,5.40e0 *costs.kpc2cm,166.e0*costs.kms))  
        return phi_tot




    def grav_profile(self, z):
        """Compute total gravity of the disk at z.


        Returns:
            float/array: potentials
        """

        
        zz = z
        
        dphi_d = (  self.dphi_d_thin( zz, 2.106e10*costs.M_sun, 3.859e0*costs.kpc2cm, 2.162e10*costs.M_sun, 9.052e0*costs.kpc2cm,-1.704e10*costs.M_sun, 3.107e0*costs.kpc2cm, 0.243e0*costs.kpc2cm) + 
                    self.dphi_d_thick(zz, 0.056e10*costs.M_sun, 0.993e0*costs.kpc2cm, 3.766e10*costs.M_sun, 6.555e0*costs.kpc2cm,-3.250e10*costs.M_sun, 7.651e0*costs.kpc2cm, 0.776e0*costs.kpc2cm) + 
                    self.dphi_d_HI(   zz, 2.046e10*costs.M_sun, 9.021e0*costs.kpc2cm, 2.169e10*costs.M_sun, 9.143e0*costs.kpc2cm,-3.049e10*costs.M_sun, 7.758e0*costs.kpc2cm, 0.168e0*costs.kpc2cm) + 
                    self.dphi_d_H2(   zz, 0.928e10*costs.M_sun, 6.062e0*costs.kpc2cm, 0.163e10*costs.M_sun, 3.141e0*costs.kpc2cm,-0.837e10*costs.M_sun, 4.485e0*costs.kpc2cm, 0.128e0*costs.kpc2cm))
        dphi_tot = dphi_d +(self.dphi_b(zz,2.61e10*costs.M_sun,  0.44e0*costs.kpc2cm) + 
                            self.dphi_h(zz,5.40e0 *costs.kpc2cm,166.e0*costs.kms))
        return -dphi_tot

    # potentials
    @staticmethod
    def zeta(z,b):
        return np.sqrt(z**2 + b**2)
    
    def phi_MN_1(self, z, M, a, b, R0 = None):
        if R0 is None: R0 = self.R0
        return (-costs.factG_in_cgs*M)/(np.sqrt(R0**2+(a+self.zeta(z,b))**2))
    
    def phi_MN_2(self, z, M, a, b, R0 = None):
        if R0 is None: R0 = self.R0
        return self.phi_MN_1(z,M,a,b)*(1.e0+a*(a+self.zeta(z,b))/(R0**2+(a+self.zeta(z,b))**2))
    
    def phi_MN_3(self, z, M, a, b, R0 = None):
        if R0 is None: R0 = self.R0
        return  self.phi_MN_1(z,M,a,b)*(1.e0+a*(a+self.zeta(z,b))/(R0**2+(a+self.zeta(z,b))**2)-1.e0/3.e0*a**2*(R0**2-2.e0*(a+self.zeta(z,b))**2)/(R0**2+(a+self.zeta(z,b))**2)**2)
    
    def phi_d_thin(self, z, M1, a1, M2, a2, M3, a3, b):
        return self.phi_MN_3(z,M1,a1,b) + self.phi_MN_3(z,M2,a2,b) + self.phi_MN_3(z,M3,a3,b)
    
    def phi_d_thick(self, z, M1, a1, M2, a2, M3, a3, b):
        return self.phi_MN_1(z,M1,a1,b) + self.phi_MN_1(z,M2,a2,b) + self.phi_MN_1(z,M3,a3,b)
    
    def phi_d_HI(self, z, M1, a1, M2, a2, M3, a3, b):
        return self.phi_MN_2(z,M1,a1,b) + self.phi_MN_2(z,M2,a2,b) + self.phi_MN_2(z,M3,a3,b)
    
    def phi_d_H2(self, z, M1, a1, M2, a2, M3, a3, b):
        return self.phi_MN_3(z,M1,a1,b) + self.phi_MN_3(z,M2,a2,b) + self.phi_MN_3(z,M3,a3,b)
    
    def phi_b(self, z, M_b, a_b, R0 = None):
        if R0 is None: R0 = self.R0
        return (-costs.factG_in_cgs*M_b)/(np.sqrt(R0**2+z**2)+a_b)
    
    def phi_h(self, z, r_h, v_h, R0 = None):
        if R0 is None: R0 = self.R0
        return 0.5e0*v_h**2*np.log((R0/r_h)**2+(z/r_h)**2+1.e0)
    
    # gravitational accelerations

    def phi_MN_1(self, z,M,a,b, R0 = None):
        if R0 is None: R0 = self.R0
        return (-costs.factG_in_cgs*M)/(np.sqrt(R0**2+(a+self.zeta(z,b))**2))
   
    def dphi_MN_1(self, z,M,a,b, R0 = None):
        if R0 is None: R0 = self.R0
   
        return costs.factG_in_cgs*M*z*(a+self.zeta(z,b))/(self.zeta(z,b)*(np.sqrt(R0**2+(a+self.zeta(z,b))**2))**3)
    
    def dphi_MN_2(self, z,M,a,b, R0 = None):
        if R0 is None: R0 = self.R0
        return (
            self.dphi_MN_1(z,M,a,b)*(1.e0+a*(a+self.zeta(z,b))/(R0**2+(a+self.zeta(z,b))**2))+
        self.phi_MN_1(z,M,a,b)*(a*z/self.zeta(z,b)*(R0**2-(a+self.zeta(z,b))**2)/(R0**2+(a+self.zeta(z,b))**2)**2)
        )
   
    def dphi_MN_3(self, z,M,a,b, R0 = None):
        if R0 is None: R0 = self.R0
        return (self.dphi_MN_1(z,M,a,b)*(1.e0+a*(a+self.zeta(z,b))/(R0**2+(a+self.zeta(z,b))**2)-
                1.e0/3.e0*a**2*(R0**2-2.e0*(a+self.zeta(z,b))**2)/(R0**2+(a+self.zeta(z,b))**2)**2)+
                self.phi_MN_1(z,M,a,b)*(a*z/self.zeta(z,b)*(R0**2-(a+self.zeta(z,b))**2)/(R0**2+(a+self.zeta(z,b))**2)**2+
                4.e0*a**2*(a+self.zeta(z,b))*z/(3.e0*self.zeta(z,b))*
                (2.e0*R0**2-(a+self.zeta(z,b))**2)/(R0**2+(a+self.zeta(z,b))**2)**3))


    # disks 
    def dphi_d_thin(self, z,M1,a1,M2,a2,M3,a3,b):
        return self.dphi_MN_3(z,M1,a1,b) + self.dphi_MN_3(z,M2,a2,b) + self.dphi_MN_3(z,M3,a3,b)
   
    def dphi_d_thick(self, z,M1,a1,M2,a2,M3,a3,b):
        return self.dphi_MN_1(z,M1,a1,b) + self.dphi_MN_1(z,M2,a2,b) + self.dphi_MN_1(z,M3,a3,b)
    
    def dphi_d_HI(self, z,M1,a1,M2,a2,M3,a3,b):
        return self.dphi_MN_2(z,M1,a1,b) + self.dphi_MN_2(z,M2,a2,b) + self.dphi_MN_2(z,M3,a3,b)
    
    def dphi_d_H2(self, z,M1,a1,M2,a2,M3,a3,b):
        return self.dphi_MN_3(z,M1,a1,b) + self.dphi_MN_3(z,M2,a2,b) + self.dphi_MN_3(z,M3,a3,b)
    
    def dphi_b(self, z,M_b,a_b, R0 = None):
        if R0 is None: R0 = self.R0
        return (costs.factG_in_cgs*M_b*z)/(np.sqrt(R0**2+z**2)*(np.sqrt(R0**2+z**2)+a_b)**2)
    
    def dphi_h(self, z,r_h,v_h, R0 = None):
        if R0 is None: R0 = self.R0
        return v_h**2*z/(R0**2+z**2+r_h**2)
   