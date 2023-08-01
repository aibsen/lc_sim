import numpy as np
import sncosmo
import os
from scipy.stats import truncnorm
import csv
from astropy.cosmology import Planck18


# model, z_range,parameters,rate,flux_scale

def SFRate(z,mode):
    # Madau & Dickinson 2014
    if mode == 'md14':
        a, b, c, d = 0.015, 2.9, 2.7, 5.6
    # Strogler 2015
    elif mode == 's15':
        a, b, c, d = 0.015, 1.50, 5.0, 6.10
    sfr = (a*(z+1)**c)/(1+((1+z)/b)**d)
    return sfr

def SNRate(z,R0,mode):
    if mode == 'md14':
        return R0*SFRate(z,mode)/SFRate(0,mode)
    elif mode == 's15':
        h = Planck18.H0.value/100.0
        k = R0/SFRate(0,mode)
        return k*(h*h)*SFRate(z,mode)

def Iax_mag(template_index,sed, **kwargs):
    # print(template_index)
    required_sed = sed[template_index]
    f = open('../SED/SIMSED.Iax/SED.INFO')
    sed_data = list(csv.reader(f))[12:]
    data = [d[0] for d in sed_data]
    data = list(filter(lambda d: required_sed in d ,data))
    mean_mag = float(data[0].split()[3])
    a = (-18-mean_mag)/0.5
    b = (-13-mean_mag)/0.5
    mabs = truncnorm.rvs(a=a, b=b, loc=mean_mag, scale=0.5)
    return mabs

transients = {
    'Ia':  {'rate': lambda z: (2.5e-5)*(1+z)**1.5,
            'z_range':(0.011,0.2),
            'built-in':True,
            'template':'salt2'}, 
    'Ibc': {'rate': lambda z: SNRate(z,1.9e-5,'s15'), 
            'z_range':(0.011,0.2),
            'built-in':True,
            'template':'nugent'},  
    'IIP': {'rate': lambda z: SNRate(z,1.4e-4,'s15'), 
            'z_range':(0.011,0.2),
            'built-in':True,
            'template':'nugent'},
    'IIn': {'rate': lambda z: SNRate(z,7.5e-6,'s15'), 
        'z_range':(0.011,0.2),
        'built-in':True,
        'template':'nugent'},
    # 'Iax': {'rate': lambda z: SNRate(z,6.00e-06,'md14'),
    #         'z_range':(0.011,0.2),
    #         'parameters':{'mag':Iax_mag},
    #         'built-in':False,
    #         'flux_scale': 1},
    # 'Ia-91bg':{
    #         'rate': lambda z: (3e-6)*(1+z)**1.5,
    #         'z_range':(0.011,0.2),
    #         'parameters':{'mag': lambda **m : truncnorm.rvs(a=-2.5, b=np.inf, loc=-17.5, scale=0.2)},
    #         'built-in':False,
    #         'flux_scale':1 },
    'SLSN':{'rate': lambda z: SNRate(z,2.00e-08,'md14'),
            'z_range':(0.011,0.2),
            'parameters':{'mag': lambda **m: truncnorm.rvs(a=-0.75, b=3.52, loc=-22.47, scale=0.7)},
            'built-in':False,
            'flux_scale':8.358E-41} 
    }

def load_fields_ccd(fields_ccd_dir=''):

    fields_file=os.path.join(fields_ccd_dir,'ztf_fields.txt')
    ccd_file=os.path.join(fields_ccd_dir,'ztf_ccd_corners.txt')
    fields_raw = np.genfromtxt(fields_file, comments='%')

    fields = {'field_id': np.array(fields_raw[:,0], dtype=int),
              'ra': fields_raw[:,1],
              'dec': fields_raw[:,2]}

    ccd_corners = np.genfromtxt(ccd_file, skip_header=1)
    ccds = [ccd_corners[np.array([0,1,3,2])+4*k, :2] for k in range(16)]

    return fields, ccds

def load_ztf_bands(bandpass_dir=''):
    bands = {
        'ztfg' : 'ztfg_eff.txt',
        'ztfr' : 'ztfr_eff.txt'
        # 'ztfi' : 'ztfi_eff.txt',
    }

    for bandname in bands.keys() :
        fname = bands[bandname]
        b = np.loadtxt(os.path.join(bandpass_dir,fname))
        band = sncosmo.Bandpass(b[:,0], b[:,1], name=bandname)
        sncosmo.registry.register(band, force=True)


