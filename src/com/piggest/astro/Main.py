'''
Created on 2018年11月6日

@author: yaoyao
'''
from com.piggest.astro.Galaxy import Galaxy

if __name__ == '__main__':
    gal = Galaxy("I Zw 18")
    gal.download_sdss_spec()
    spec = gal.get_spec()
    spec.interp1d().reduce_mw_extinct().reduce_redshift()
    print("金属丰度=",spec.get_z())