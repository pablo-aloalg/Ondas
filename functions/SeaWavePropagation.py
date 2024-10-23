import numpy as np
from scipy.optimize import newton
from scipy.constants import g

LIMIT_KH = 700.
LIMIT_ANGULO = 80.
LIMIT_HS = .2
LIMIT_T = 2.5

class Dispersion(object):
    def __init__(self, T, h='ind'):
        self.T = T
        self.omega = 2.*np.pi/self.T
        self.Lo = g*T**2./2./np.pi
        self.ko = 2.*np.pi/self.Lo
        if isinstance(h, str):
            if h.strip().lower()=='ind':
                # Estamos en indefinidas
                self.L = self.Lo
                self.k = 2.*np.pi/self.L
                self.n = .5
                self.h = 1e6
            else:
                raise Exception("Cadena de profundidad no reconocida")
        else:
            self.h = h
            self.k = (newton(self.dispersion, x0=(self.ko), args=(self.omega, h),
                fprime=self.deriv_disper, maxiter=450))
            # print("El valor resuelto para k es {0:.3f}".format(self.k))
            if self.k == 0:
                raise Exception("Número de onda nulo. Abortando...")
            self.L = 2.*np.pi/self.k
            if self.k*self.h > LIMIT_KH:
                self.n = .5
            else:
                self.n = .5 * (1 + (2 * self.k * h) / np.sinh(2 * self.k * h))
            
        self.c = self.L / self.T
        self.cg = self.n * self.c

    def __str__(self):
        cadena = "T = {0:5.2f}, L = {1:6.2f}, h = {2:8.2f}, c = {3:5.2f}"
        cadena += "\n"
        cadena += "O = {4:5.2f}, k = {5:6.2f}, n = {6:5.2f}, cg = {7:5.2f}"
        return cadena.format(self.T, self.L, self.h, self.c, self.omega, 
            self.k, self.n, self.cg)

    @staticmethod
    def dispersion(k, omega, h):
        return omega**2. - g*k*np.tanh(k*h)

    @staticmethod
    def dispersion_logk(logk, omega, h):
        return omega**2. - g*np.exp(logk)*np.tanh(np.exp(logk)*h)

    @staticmethod
    def deriv_disper(k, omega, h):
        kh = k*h
        if kh > LIMIT_KH:
            return -g * np.tanh(kh)
        else:
            return -g*kh/np.cosh(kh)**2. -g*np.tanh(kh)

    @staticmethod
    def deriv_disper_logk(logk, omega, h):
        kh = np.exp(logk)*h
        if kh > LIMIT_KH:
            return -g * np.tanh(kh)
        else:
            return -g*kh/np.cosh(kh)**2. -g*np.tanh(kh)

class PuntoOla(object):
    def __init__(self, H0, Dir, Tp, h='ind'):
        self.H0 = H0
        self.T = Tp
        self.theta0_d = Dir
        self.theta0_r = np.deg2rad(Dir)
        self.onda = Dispersion(Tp, h)

    def propaga(self, hd, DirCosta):
        # self.DifAngulo_d = 180. - abs(abs(DirCosta - self.theta0_d) - 180.)
        self.DifAngulo_d = DirCosta - self.theta0_d
        self.DifAngulo_r = np.deg2rad(self.DifAngulo_d)
        if np.abs(self.DifAngulo_d) > LIMIT_ANGULO:
            return (0, 0)
        destino = Dispersion(self.T, hd)
        thetad_r = np.arcsin(destino.c*np.sin(self.DifAngulo_r)/self.onda.c)
        Hd  = self.H0 * np.sqrt(self.onda.cg / destino.cg) * np.sqrt(np.cos(self.DifAngulo_r) / np.cos(thetad_r))
        return(Hd, np.rad2deg(thetad_r))


    def Goda(self, hb, tanbeta):
        return 0.17 * self.onda.Lo * (1. - np.exp(-1.5*np.pi*hb*(1.+15.*(tanbeta)**(4./3.))/self.onda.Lo))

    def condicion_rotura(self, hd, tanbeta):
        # HGoda = 0.17 * self.onda.Lo * (1. - np.exp(-1.5*np.pi*hd*(1.+15.*(tanbeta)**(4./3.))/self.onda.Lo))
        Hlimite = .8 * hd
        destino = Dispersion(self.T, hd)
        thetad_r = np.arcsin(destino.c*np.sin(self.DifAngulo_r)/self.onda.c)
        Hd = self.H0 * np.sqrt(self.onda.cg / destino.cg) * np.sqrt(np.cos(self.DifAngulo_r) / np.cos(thetad_r))
        return Hlimite - Hd

    def condicion_rotura_logh(self, logh, tanbeta):
        h = np.exp(logh)
        Hlimite = .8 * h
        # HGoda = 0.17 * self.onda.Lo * (1. - np.exp(-1.5*np.pi*h*(1.+15.*(tanbeta)**(4./3.))/self.onda.Lo))
        destino = Dispersion(self.T, h)
        thetad_r = np.arcsin(destino.c*np.sin(self.DifAngulo_r)/self.onda.c)
        Hd = self.H0 * np.sqrt(self.onda.cg / destino.cg) * np.sqrt(np.cos(self.DifAngulo_r) / np.cos(thetad_r))
        return Hlimite - Hd

    def resuelve_rotura(self, tanbeta, DirCosta):
        self.DifAngulo_d = DirCosta - self.theta0_d
        self.DifAngulo_r = np.deg2rad(self.DifAngulo_d)
        if np.abs(self.DifAngulo_d) > LIMIT_ANGULO:
            self.Hb, self.hb, self.thetab_d = (0., 0., 0.)
            self.thetab_r = 0.
            return (0., 0., 0.)
        elif (self.H0 < LIMIT_HS) or (self.onda.T < LIMIT_T):
            # Si no hay rotura en un palmo de agua es que no hay rotura
            self.Hb, self.hb, self.thetab_d = (0., 0., 0.)
            self.thetab_r = 0.
            return (0., 0., 0.)
        else:
            self.hb = (newton(self.condicion_rotura,
                x0 = .25, 
                args=(tanbeta, )))
            if self.hb > 10. * self.H0:
                self.hb = 1.25 * self.H0
            self.Hb, self.thetab_d = self.propaga(self.hb, DirCosta)
            self.Hrmsb = self.Hb / 1.4
        self.thetab_r = np.deg2rad(self.thetab_d)
        return (self.hb, self.Hb, self.thetab_d)

    def transporte_CERC(self, D50=0.3, p=0.4, rhow=1025, rhos=2650):
        #K = 1.4 * np.exp(-2.5*D50)
        K = 0.39
        if (self.hb == 0) or (self.Hb == 0):
            return 0.
        gamma = self.Hb / self.hb
        # return (K*rhow*g**.5/16*(rhow-rhos)*(1.-p)*np.sqrt(sigma))*(self.Hb**(5./2.)*np.sin(2.*self.thetab_r))
        return (K * rhow * g**.5 * self.Hrmsb**2.5 * np.sin(2. * self.thetab_r) 
            / (16. * (rhos - rhow) * (1. - p) * gamma**.5) )
