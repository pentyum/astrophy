'''
Created on 2018年11月6日

@author: yaoyao
'''
import numpy as np
from scipy.interpolate import interp1d
import sfdmap
import matplotlib.pyplot as plt
import cn.edu.ustc.zesenlin.fit_line as fl
from cn.edu.ustc.zesenlin.my_ccm89 import ccm89

class Spectra:
    objid = None
    hdulist = None
    flux_list = None
    wavelen_list = None
    flux_e_list = None
    __t3 = 0
    __t2 = 0
    
    def __init__(self,objid,hdulist):
        self.hdulist=hdulist
        self.objid=objid
        self.wavelen_list=np.power(10,hdulist[1].data['loglam'])
        self.flux_list=hdulist[1].data['flux']
        self.flux_e_list=1./np.sqrt(hdulist[1].data['ivar'] + 1.e-8)
    
    def interp1d(self):
        '''
        插值
        '''
        flux_list = self.flux_list
        flux_e_list = self.flux_e_list
        wave = self.wavelen_list
        self.wavelen_list = np.arange(min(wave), max(wave), 1.)
        f_interp_flux = interp1d(wave, flux_list)
        self.flux_list = f_interp_flux(self.wavelen_list)
        f_interp_err = interp1d(wave, flux_e_list)
        self.flux_e_list = f_interp_err(self.wavelen_list)
        return self
        
    def reduce_mw_extinct(self):
        '''
        银河消光改正
        '''
        wav = self.wavelen_list
        flux = self.flux_list
        flux_err = self.flux_e_list
        m = sfdmap.SFDMap('./sfddata-master')
        ebv = m.ebv(self.get_ra(), self.get_dec())
        R_V = 3.1
        A_V = ebv * R_V
        A_lam = ccm89(wav, A_V, R_V)
        fac = 10**(-0.4 * A_lam)
        flux_deredden = flux / fac
        flux_err_deredden = flux_err / fac
        self.flux_list = flux_deredden
        self.flux_e_list = flux_err_deredden
        return self
    
    def reduce_redshift(self):
        '''
        修正红移
        '''
        z = self.get_redshift()
        
        self.wavelen_list = self.wavelen_list / (1. + z)
        self.flux_list = self.flux_list * (1. + z)
        self.flux_e_list = self.flux_e_list * (1. + z)
        return self
    
    def get_line(self,line_name):
        '''
        获得发射线流量
        '''
        res = fl.fit_line(self.wavelen_list, self.flux_list, self.flux_e_list, lines=[0, 18, 25, 29, 46, 52], isplot=True)
        for line in res:
            if line.name == line_name:
                return line
        return None
    
    def __get_t3(self):
        if self.__t3 == 0:
            l4363 = self.get_line("[OIII]4363")
            l4959 = self.get_line("[OIII]4959")
            l5007 = self.get_line("[OIII]5007")
            lg = np.log10((l4959.flux+l5007.flux)/l4363.flux)
            t3 = 1
            i = 0
            while i<=10:
                CT = 8.44-1.09*t3+0.5*t3**2-0.08*t3**3
                t3 = 1.432/(lg-np.log10(CT))
                i=i+1
            self.__t3 = t3
            return t3;
        else:
            return self.__t3
    
    def get_z3(self):
        t3 = self.__get_t3()
        lHb = self.get_line("Hb")
        l4959 = self.get_line("[OIII]4959")
        l5007 = self.get_line("[OIII]5007")
        return np.log10((l4959.flux+l5007.flux)/lHb.flux)+6.2+1.251/t3-0.55*np.log10(t3)-0.014*t3
    
    def __get_t2(self):
        if self.__t2 == 0:
            t3 = self.__get_t3()
            z3 = self.get_z3()
            if z3<7.2:
                t2 = -0.577+t3*(2.065-0.489*t3)
            elif z3>8.2:
                t2 = 2.967+t3*(-4.797+2.827*t3)
            else:
                t2 = -0.744+t3*(2.338-0.61*t3)
            self.__t2 = t2
            return t2
        else:
            return self.__t2
    
    def get_z2(self):
        t2 = self.__get_t2()
        lHb = self.get_line("Hb")
        l7325 = self.get_line("[OII]7325")
        if l7325.flux/l7325.flux_err < 5:
            print("Warning: l7325.flux/l7325.flux_e < 5 !")
        return np.log10(2*l7325.flux/lHb.flux)+6.901+2.487/t2-0.483*np.log10(t2)-0.013*t2
    
    def get_z(self):
        '''
        电子温度法求金属丰度
        '''
        return 12+np.log10(10**(self.get_z2()-12)+10**(self.get_z3()-12))
        
    def get_ra(self):
        return self.hdulist[2].data['PLUG_RA'][0]
    
    def get_dec(self):
        return self.hdulist[2].data['PLUG_DEC'][0]
    
    def get_redshift(self):
        return self.hdulist[2].data['Z'][0]
    
    def plot(self,start=0,end=0):
        plt.plot(self.wavelen_list,self.flux_list)
        if start != 0 and end != 0:
            plt.xlim((start, end))