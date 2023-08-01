import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import matplotlib.patheffects as PathEffects
from matplotlib.colors import LinearSegmentedColormap


# --------------------------- #
# - Plots Methods           - #
# --------------------------- #\

data_dir = "lcs/"
_d2r = np.pi / 180 

def convert_radec_azel(ra, dec, edge=0):
    if edge > 0:
        if type(ra) == float:
            if ra < -180 + edge:
                ra = -180 + edge
            elif ra > 180 - edge:
                ra = 180 - edge
        else:
            ra[ra < -180 + edge] = -180 + edge
            ra[ra > 180 - edge] = 180 - edge

    ra = ((ra + 180) % 360) - 180
    ra *= -1

    az = _d2r * ra
    el = _d2r * dec

    return az, el

def plot_transients(ra,dec):
    fig_size = (12,8)
    rect = [0.1, 0.1, 0.8, 0.8]
    fig = plt.figure(figsize=fig_size)
    ax = ax = fig.add_axes(rect, projection='mollweide')
    ax.grid(True)

    # ra_range=(0,360) 
    # ra = np.linspace(ra_range[0],ra_range[1],1000)
    # dec_range=(-90,90)
    # dec = np.linspace(dec_range[0],dec_range[1],1000)
    az, el = convert_radec_azel(ra, dec)
    pl = ax.scatter(az, el)
    plt.show()
    


def show_sky_coverage(plan):
    plan_public = plan[(plan.comment=='all_sky')|(plan.comment=='nightly_plane')]
    plan_hc = plan[plan.comment=='high_cadence']
    rect = [0.1, 0.1, 0.8, 0.8]
    fig = plt.figure()
    ax = fig.add_axes(rect, projection='mollweide')

    
    # az, el = convert_radec_azel(ra, dec)
    # axcar = ax.insert_ax('bottom',shrunk=0.85,space=0,axspace=0.08)
    bsize = 26
    cup = 300
    ra_bins = np.linspace(-180,180,bsize)
    dec_bins = np.linspace(-90,90,bsize)
    c = lambda r: -r if r<= 180 else 360-r

    ra = plan_public.ra.values
    dec = plan_public.dec.values
 

    ra_conv = np.array([c(r) for r in ra])
    hist, ra_edge, dec_edge = np.histogram2d(ra_conv,dec,bins=[ra_bins,dec_bins])
    # print(ra_edge)
    # ra_edge_neg = np.linspace(-360,360,bsize)
    # az_edge, el_edge = convert_radec_azel(ra_edge_neg, dec_edge)
    # print(hist.shape)
    # print(hist)
    hist=np.swapaxes(hist,1,0)
    hist = np.where(hist<cup,hist,cup)
    cmap = plt.cm.Blues

    if not plan_hc.empty:
        ra_hc = plan_hc.ra.values
        dec_hc = plan_hc.dec.values
        ra_hc_conv = np.array([c(r) for r in ra_hc])
        hist_hc, _, _ = np.histogram2d(ra_hc_conv,dec_hc,bins=[ra_bins,dec_bins])
        hist_hc= np.swapaxes(hist_hc,1,0) 

        hist_hc_nonzero=hist_hc.nonzero()
        min_lat = hist_hc_nonzero[0].min()
        min_lon = hist_hc_nonzero[1].min()
        lat_shape = hist_hc_nonzero[0].max() - hist_hc_nonzero[0].min() +1
        lon_shape = hist_hc_nonzero[1].max() - hist_hc_nonzero[1].min() +1
        hist_hc_0 = np.zeros((lat_shape,lon_shape))
        hc_points = [(hist_hc_nonzero[0][i],hist_hc_nonzero[1][i]) for i in range(hist_hc_nonzero[0].shape[0])]
        min_point = min_lat if min_lat<min_lon else min_lon
        for x,y in hc_points:
            hist_hc_0[x-min_lat][y-min_lon] = hist_hc[x][y]

        # get colormap
        ncolors = 150
        color_array = plt.get_cmap('Oranges')(range(ncolors))
        # change alpha values
        color_array[:,-1] = np.linspace(0.0,1.0,ncolors)
        # create a colormap object
        map_object = LinearSegmentedColormap.from_list(name='cmap_alpha',colors=color_array)
        # register this new colormap with matplotlib
        plt.register_cmap(cmap=map_object)
        lon_hc = np.linspace(-np.pi, np.pi,bsize)[min_point:min_point+hist_hc_0.shape[1]+1]
        lat_hc = np.linspace(-np.pi/2., np.pi/2.,bsize)[min_point:min_point+hist_hc_0.shape[0]+1]



    lon = np.linspace(-np.pi, np.pi,bsize)
    lat = np.linspace(-np.pi/2., np.pi/2.,bsize)

    im=ax.pcolormesh(lon,lat,hist, cmap=cmap)
    tick_labels = np.array([150, 120, 90, 60, 30, 0, 330, 300, 270, 240, 210])
    ax.set_xticklabels(tick_labels, color='black',
        size=8,path_effects=[PathEffects.withStroke(linewidth=2.5, foreground='w')])
    
    if not plan_hc.empty:
        im2=ax.pcolormesh(lon_hc,lat_hc,hist_hc_0, cmap='cmap_alpha')
        plt.colorbar(im,ax=ax,location='bottom',shrink=0.5,anchor=(0.5,2.7) )
        plt.colorbar(im2,ax=ax,location='bottom',shrink=0.5,anchor=(0.5,1.5))

    else :
        plt.colorbar(im,ax=ax,location='bottom',shrink=0.5,anchor=(0.5,1.5) )

    ax.grid(True)
    plt.show()

plan = pd.read_csv('plan_pointings_v8.csv')
# plan = pd.read_csv('surveyplan_msip.csv', sep=" ")
fields = pd.read_csv('../ztf_data/ztf_fields.txt')
# print(plan)
# print(fields)
# fields.apply(lambda r: r.values.astype('str')[0].split(" "))

# plan

# plan = plan[plan.comment=='public']
print(plan.comment.unique())

plan = plan[(plan.band=='ztfg') | (plan.band=='ztfr') | (plan.band==0) | (plan.band==1)]
print(plan.band.unique())

field_dict = {"ra":[], "dec":[], "field":[]}
for field in fields.iterrows():
    items = np.array(field[1].values.astype('str')[0].split("  "))
    blank_mask = (items != "")
    items = items[blank_mask]
    field_dict["ra"].append(float((items[1])))
    field_dict["dec"].append(float(items[2]))
    field_dict["field"].append(items[-1])
    
field_df = pd.DataFrame(field_dict)


x = plan.join(field_df, on='field', how='left',lsuffix="l")
# print(x.shape)
# ra = x.ra.values
# dec = x.dec.values

show_sky_coverage(x)

