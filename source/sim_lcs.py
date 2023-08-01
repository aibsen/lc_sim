from statistics import mean
import numpy as np
import pandas as pd
import os
import sncosmo
import simsurvey
from astropy.cosmology import Planck18
from astropy.table import Table
from sim_utils import *
import simsurvey.utils as utils



sfd98_dir = '../sfddata-master'

def random_parameters(redshifts, model, n_templates, mag,
                      r_v=2., ebv_rate=0.11,
                      cosmo=Planck18,               
                      **kwargs):
    out = {}
    amp = []
    # print(n_templates)
    if n_templates > 1:
        out['template_index'] = np.random.randint(0, n_templates, len(redshifts))
        model_kw = [{'template_index': k} for k in out['template_index']]
    elif n_templates == 1:
        model_kw = [{} for z in redshifts]

    # Amplitude
    for z, model_kw_ in zip(redshifts, model_kw):
        mabs = mag(template_index=model_kw_['template_index'],**kwargs)
        model.set(z=z, **model_kw_)
        model.set_source_peakabsmag(mabs, 'bessellb', 'vega', cosmo=cosmo)
        amp.append(model.get('amplitude'))

    out['amplitude'] = np.array(amp)
    out['hostr_v'] = r_v * np.ones(len(redshifts))
    out['hostebv'] =  np.random.exponential(ebv_rate, len(redshifts))
    return out



def generate_lightcurves(transient_name,transient,mjd_range,plan,real_nights,**kwargs):
# def generate_lightcurves(transient_name,transient,mjd_range,plan,**kwargs):

    if transient["built-in"]:
        print(transient)
        transientprop = None
        template = transient['template']
    else:
        template = None
        transientprop = {
            'lcmodel': transient["model"],
            'lcsimul_func': random_parameters,
            'lcsimul_prop': transient["parameters"]
        }

    redshift_range = transient['z_range']

    tr = simsurvey.get_transient_generator(redshift_range,
                                        ratefunc=transient['rate'],
                                        ra_range=(0,360),
                                        dec_range=(-30,90),
                                        mjd_range=(mjd_range[0],mjd_range[1]),
                                        transient=transient_name,
                                        transientprop=transientprop,
                                        template=template,
                                        sfd98_dir=sfd98_dir,
                                        cosmo=Planck18
                                        )

    survey = simsurvey.SimulSurvey(generator=tr, plan=plan)
    lcs = survey.get_lightcurves(progress_bar=True)

    def filterfunc(lc):
        mask = np.array([(int(t_) in real_nights) and (t_ < lc.meta['t0'] + 100) for t_ in lc['time']])
        return lc[mask]

    lcs = lcs.filter(filterfunc)
    
    lcs.save("lcs/lc_{}_test_3.pkl".format(transient_name))


def run_simulation():

    fields, ccds = load_fields_ccd('../ztf_data')
    load_ztf_bands('../ztf_data')
    # survey_file = 'test_schedule_v8_msip.db'
    # survey_file = 'test_schedule_v6.db'
    
    # plan = simsurvey.SurveyPlan(load_opsim=survey_file, band_dict={'g': 'ztfg', 'r': 'ztfr'}, ccds=ccds)
    # plan = simsurvey.SurveyPlan(load_opsim=survey_file, band_dict={'g': 'ztfg', 'r': 'ztfr', 'i': 'desi'}, ccds=ccds)

    # print(obs['band'].shape)
    # obs = Table.read("plan_sim_paper.csv", format="ascii.csv")

    survey_file = 'test_schedule_v8_msip.db'
    plan = simsurvey.SurveyPlan(load_opsim=survey_file, band_dict={'g': 'ztfg', 'r': 'ztfr', 'i': 'desi'}, ccds=ccds)
    # plan = simsurvey.SurveyPlan(time=obs['time'], band=obs['band'], obs_field=obs['field'],
    #                         skynoise=obs['skynoise'], comment=obs['comment'],
    #                         fields={k: v for k, v in fields.items()
    #                                 if k in ['ra', 'dec', 'field_id',
    #                                             'width', 'height']},
    #                         ccds=ccds)



    mjd_range = (plan.cadence['time'].min() - 30, plan.cadence['time'].max() + 30)
    # mjd_range = (plan.pointings['time'].min(), plan.pointings['time'].max())

    real_nights = np.genfromtxt('hours_per_obsnight.dat')
    idx = np.concatenate((range(31,365), range(31)))
    # rn_2016 = np.where(real_nights[idx, -3] > 3.5)[0] +2458151  
    rn_2016 = np.where(real_nights[idx, -3] > 3.5)[0] + int(plan.pointings["time"].min())
    

    for transient_name, transient in utils.transients.items():

        print(transient_name)

        if not transient['built-in']:
            sed_dir =  '../SED/SIMSED.{}/'.format(transient_name)
            files =  next(os.walk(sed_dir))[2]
            n_templates = len(files)-1
            data_sed = {'p': [], 'w': [], 'f': [], 'sed_name':[]}

            scale = transient['flux_scale']
            for sed in files:
                if sed != "SED.INFO" :
                    print(sed)
                    p, w, f = sncosmo.read_griddata_ascii(sed_dir+sed)
                    data_sed['p'].append(p)
                    data_sed['w'].append(w)
                    data_sed['f'].append(f*scale)
                    data_sed['sed_name'].append(sed)

            transient["parameters"]['sed'] = data_sed['sed_name']
            transient["parameters"]['n_templates'] = n_templates
            # print(transient['parameters'])
            # print(transient['parameters']['mag'])

            source = simsurvey.MultiSource(data_sed['p'], data_sed['w'], data_sed['f'])

            transient["model"] = sncosmo.Model(
                source=source,
                effects=[sncosmo.CCM89Dust(),sncosmo.CCM89Dust()],
                effect_names=['host','mw'],
                effect_frames=['rest', 'obs']
            )

        generate_lightcurves(transient_name,transient,mjd_range,plan,rn_2016)


run_simulation()

