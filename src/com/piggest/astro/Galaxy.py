'''
Created on 2018年10月24日

@author: yaoyao
'''
import astropy.coordinates as coord
import astropy.units as u
from urllib.parse import urlencode
from urllib.request import urlopen
from astroquery.sdss import SDSS
import numpy as np
import cv2
import IPython.display as ipydp
from com.piggest.astro.Spectra import Spectra

class Galaxy:
    loc = None
    name = None
    photo_data_table = None
    spec_data_table = None
    spec = None
    info_table = None
    objid = 0
    r = 2
    
    def __init__(self,name):
        if isinstance(name,coord.sky_coordinate.SkyCoord):
            self.loc=name
        else:
            self.loc = coord.SkyCoord.from_name(name)
            self.name = name
        self.download_sdss_data()
        
    def download_sdss_data(self):
        sdss = SDSS()
        self.info_table = sdss.query_region(self.loc,radius=self.r*u.arcsec)
        self.photo_data_table = sdss.query_region(self.loc,radius=self.r*u.arcsec,photoobj_fields=['u','g','r','i','z'])
        while self.info_table is None and self.r<=10:
            self.r=self.r+0.4
            self.info_table = sdss.query_region(self.loc,radius=self.r*u.arcsec)
            self.photo_data_table = sdss.query_region(self.loc,radius=self.r*u.arcsec,photoobj_fields=['u','g','r','i','z'])
        self.objid = self.info_table['objid'][0]
        #print("r=%f"%(self.r))
    
    def download_sdss_spec(self,add_r=0):
        sdss = SDSS()
        self.spec_data_table = sdss.query_region(self.loc,radius=(self.r+add_r)*u.arcsec,spectro=True)
        if self.has_spec() == True:
            print("在r=%f角秒内找到了%d个光谱"%(self.r+add_r,len(self.spec_data_table)))
            self.spec = SDSS.get_spectra(matches=self.spec_data_table)
            i=0
            for match in self.spec_data_table:
                self.spec[i] = Spectra(match['objid'],self.spec[i])
                i=i+1
        else:
            print("在r=%f角秒内没有找到光谱"%(self.r+add_r))
            
    def has_spec(self):
        if self.spec_data_table == None:
            return False
        else:
            return True
        
    def get_redshift(self,n=0):
        if self.has_spec() == False:
            return None
        else:
            return self.spec_data_table["z"][n]
        
    def get_spec(self,n=0):
        return self.spec[n]

    def print_spec(self,n=0):
        self.spec[n].plot()
        
    def get_sdss_img_data(self,field_arcmin,pixels): #返回rgb三维数组
        im_size = field_arcmin*u.arcmin # 视场
        im_pixels = pixels
        cutoutbaseurl = 'http://skyserver.sdss.org/dr14/SkyServerWS/ImgCutout/getjpeg'
        query_string = urlencode(dict(ra=self.loc.ra.deg,
                              dec=self.loc.dec.deg,
                              width=im_pixels, height=im_pixels,
                              scale=im_size.to(u.arcsec).value/im_pixels))
        url = cutoutbaseurl + '?' + query_string
        return urlopen(url).read()
    
    def get_sdss_img_rgb(self,field_arcmin,pixels): #输入url返回rgb三维数组
        data = self.get_sdss_img_data(field_arcmin,pixels)
        data = np.frombuffer(data, np.uint8)
        img_data = cv2.imdecode(data, cv2.IMREAD_COLOR)
        img_data = cv2.cvtColor(img_data, cv2.COLOR_BGR2RGB)  #opencv使用的格式是bgr，要转化一下
        return img_data
    
    def show_sdss_img(self,field_arcmin=3,pixels=360):
        ipydp.display(ipydp.Image(self.get_sdss_img_data(field_arcmin,pixels)))
    
    def get_photo(self,wave,n=0):
        return self.photo_data_table[wave][n];