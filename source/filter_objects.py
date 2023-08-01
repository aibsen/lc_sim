import numpy as np
import pandas as pd
import math
import sys


data_dir = "lcs/"

p = lambda f: 0 if f=="ztfr" else 1 #r=0, g=1

ztf_type_dict = {
    '0':'Ia',
    '1':'Ia-91bg', 
    '2':'Iax',  
    '3':'IIP', 
    '4':'Ibc',
    '5':'SLSN',
    '6':'IIn'
}

def filter_lightcurves(input,run_number,class_code):
    sn_data = pd.read_pickle(input)
    lcs = sn_data['lcs']
    meta = sn_data['meta']
    stats = sn_data['stats']

    metadata_keys = ['object_id','z', 't0', 'mwebv','ra','dec','mwebv_sfd98',
        'max_mag_r', 'max_mag_g','t_peak']

    filtered_meta={k:[] for k in metadata_keys}

    filtered_data = {
        'object_id':[],
        'mjd':[],
        'passband':[],
        'magpsf':[],
        'sigmagpsf':[],
        'flux':[],
        'flux_err':[]
    }

    for l in range(len(lcs)):

        z = meta['z'][l]
        t0 = meta['t0'][l]
        mwebv = meta['mwebv'][l]
        ra = meta['ra'][l]
        dec = meta['dec'][l]
        mwebv_sfd98 = meta['mwebv_sfd98'][l]
        object_id = meta['idx_orig'][l]
        max_mag_r = stats['mag_max']['ztfr'][l]
        max_mag_g = stats['mag_max']['ztfg'][l]

        # filter galactic plane and magnitude limits
        if np.abs(dec) > 7 and (max_mag_r <=19 and max_mag_g<=19):
            lc = lcs[l]
            #limit points to -+100 days before/after explosion
            mjds = np.array([point[0] for point in lc])
            time_mask = np.abs(mjds-t0)<100
            
            mjds = mjds[time_mask]
            passbands = np.array([point[1] for point in lc])[time_mask]
            flux_errs = np.array([point[3] for point in lc])[time_mask]
            fluxes = np.array([point[2] for point in lc])[time_mask]

            #points in g and r band only and convert to numeric code r=0, g=1
            passband_mask = np.where((passbands == 'ztfg') | (passbands == 'ztfr'))
            mjds = mjds[passband_mask]
            passbands = passbands[passband_mask]
            flux_errs = flux_errs[passband_mask]
            fluxes = fluxes[passband_mask]
            passbands = np.array([p(point) for point in passbands])

            #remove points with high uncertainty
            for c in range(5): #sigma clipping 5 times
                if len(flux_errs>0):
                    std_fluxerr = flux_errs.std()
                    sigma_mask = np.abs(flux_errs-flux_errs.mean()) <= 3*std_fluxerr
                    
                    mjds = mjds[sigma_mask]
                    passbands = passbands[sigma_mask]
                    flux_errs = flux_errs[sigma_mask]
                    fluxes = fluxes[sigma_mask]

            #keep positive fluxes only
            positive_mask = fluxes>0
            mjds = mjds[positive_mask]
            passbands = passbands[positive_mask]
            flux_errs = flux_errs[positive_mask]
            fluxes = fluxes[positive_mask]

            #convert fluxes to mags
            flux_to_mag = lambda f: 30-2.5*math.log10(f)
            fluxerr_to_sigmag = lambda ferr,f: np.sqrt(np.abs(2.5/math.log(10)*(ferr/f)))

            magpsfs = np.array([flux_to_mag(f) for f in fluxes])
            sigmagpsfs = np.array([fluxerr_to_sigmag(f_err,f) for f_err,f in np.stack((flux_errs,fluxes), axis = 1)])
            
            #if we still have points:
            if magpsfs.shape[0]>0:
                #find time of max_mag
                max_magpsfs = magpsfs.min()
                time_max_mag = mjds[magpsfs==max_magpsfs][0]
                
                # require at least 3 detections per band before that
                early_points = passbands[mjds <= time_max_mag]
                enough_early_points = np.unique(early_points, return_counts=True)[1]
                
                if enough_early_points.shape[0]==2 and enough_early_points[0]>=3 and enough_early_points[1]>=3:
                    object_id = (object_id+1)*10000+10*run_number+class_code
                    obj_ids = np.full(mjds.shape,object_id)
                    
                    filtered_data['object_id']+=list(obj_ids)
                    filtered_data['mjd']+=list(mjds)
                    filtered_data['passband']+=list(passbands)
                    filtered_data['magpsf']+=list(magpsfs)
                    filtered_data['sigmagpsf']+=list(sigmagpsfs)
                    filtered_data['flux']+=list(fluxes)
                    filtered_data['flux_err']+=list(flux_errs)
                    
                    array_lens = len(filtered_data['object_id'])
                    for k,v in filtered_data.items():
                        assert(len(v)==array_lens)

                    filtered_meta['object_id'].append(object_id)
                    filtered_meta['t0'].append(t0)
                    filtered_meta['z'].append(z)
                    filtered_meta['mwebv'].append(mwebv)
                    filtered_meta['ra'].append(ra)
                    filtered_meta['dec'].append(dec)
                    filtered_meta['mwebv_sfd98'].append(mwebv_sfd98)
                    filtered_meta['max_mag_g'].append(max_mag_g)
                    filtered_meta['max_mag_r'].append(max_mag_r)
                    filtered_meta['t_peak'].append(time_max_mag)

    return filtered_meta, filtered_data


def main(class_codes=[int(k) for k in ztf_type_dict], run_numbers=[0]):
    run_numbers = [0,1,2,3]
    class_codes = [0, 3, 4, 5, 6]
    for run_number in run_numbers:
        for class_code in class_codes:
            class_name = ztf_type_dict[str(class_code)]
            try: 
                input_fn = data_dir+"lc_{}_test_{}.pkl".format(class_name,run_number)
                # input_fn = data_dir+"lc_{}_{}.pkl".format(class_name,run_number)
                print(input_fn)
                meta, data = filter_lightcurves(input_fn,run_number,class_code)
                df_meta = pd.DataFrame(meta)
                df_data = pd.DataFrame(data)
                if df_meta.shape[0] > 0:
                    print(df_meta.shape)
                    # print(df_data.shape)
                    df_meta.to_csv(data_dir+'meta_{}_test_{}.csv'.format(class_name,run_number),index=False)
                    # df_meta.to_csv(data_dir+'meta_{}_{}.csv'.format(class_name,run_number),index=False)
                    df_data.to_csv(data_dir+'lcs_{}_test_{}.csv'.format(class_name,run_number), index=False)
                    # df_data.to_csv(data_dir+'lcs_{}_{}.csv'.format(class_name,run_number), index=False)
                else:
                    ("No objects left for run{} of class{}".format(run_number,class_code))
                print("")
            except FileNotFoundError:
                print("Light curve file for class {} in run {} does not exist.".format(ztf_type_dict[str(class_code)],run_number))
                print("File should be in lcs/ and have the format lc_classcode_runnumer.pkl")
                sys.exit(1)

if __name__ == "__main__":
    main()