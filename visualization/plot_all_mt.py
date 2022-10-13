import pygmt
import sys
import json
import glob

import numpy as np

import faulthandler
faulthandler.enable()


# Change to where you store results
data_dir = '/home/ammcpherson/REPOSITORIES/mtuq/'
cutoff1 = len(data_dir)

# Controls default beachball size
ballsize = '8p'

# Region of interest
ax2 = [-164, -128, 54, 72]

# Some GMT stuff
Jbscale = 15500000

xmin = ax2[0] - 1
xmax = ax2[1] + 1
ymin = ax2[2] - 1
ymax = ax2[3]
        
xcen = (xmax + xmin) / 2
ycen = (ymax + ymin) / 2

J = 'b%.1f/%.1f/%.1f/%.1f/1:%i' % (xcen,ycen,ymin,ymax,Jbscale)
R = '%.1f/%.1f/%.1f/%.1f' % (xmin,xmax,ymin,ymax)
xran = xmax - xmin
yran = ymax - ymin

xL = xmin + 0.9 * xran
yL = ymin + 0.1 * yran

yL = ymin + 0.15 * yran
sbarinfo = 'g%.1f/%.1f+w100+ab+u' % (xL,yL)

# Find all mt json files
mt_file_ext = 'DC_mt.json'
cutoff2 = len(mt_file_ext)
mt_path = data_dir + '*' + mt_file_ext
loc_file_ext = 'DC_best_origin.json'

mt_files = glob.glob(mt_path)


# Grab catalog informations
xvals, yvals, depths, mags = [], [], [], []
mrr, mtt, mff, mrt, mrf, mtf, exp = [], [], [], [], [], [], []

# Need to organize all the hypocenter and mt information
eids = []
for file in mt_files:
    eid = file[cutoff1:-cutoff2]
    eids.append(eid)

for eid in eids:
    file_mt = data_dir + eid + mt_file_ext
    file_loc = data_dir + eid + loc_file_ext

    f = open(file_mt,'r')
    data_mt = json.load(f)
    f.close()

    f = open(file_loc,'r')
    data_loc = json.load(f)
    f.close()

    xvals.append(data_loc['longitude'])
    yvals.append(data_loc['latitude'])
    depths.append(data_loc['depth_in_m'] * 1E-3)

    mrr.append(data_mt['Mrr'])
    mtt.append(data_mt['Mtt'])
    mff.append(data_mt['Mpp'])
    mrt.append(data_mt['Mrt'])
    mrf.append(data_mt['Mrp'])
    mtf.append(data_mt['Mtp'])
    exp.append(12)

    # Compute moment magnitude for each moment tensor
    M = np.array([[data_mt['Mrr'],data_mt['Mrt'],data_mt['Mrp']],
                  [data_mt['Mrt'],data_mt['Mtt'],data_mt['Mtp']],
                  [data_mt['Mrp'],data_mt['Mtp'],data_mt['Mpp']]],dtype=float)

    M0 = (np.tensordot(M,M)/2.)**0.5
    Mw = 2./3.*(np.log10(M0) - 9.1)
    mags.append(Mw)

specs = {'mrr':mrr, 'mtt':mtt, 'mff':mff, 'mrt':mrt, 'mrf':mrf, 'mtf':mtf,'exponent':exp}

#debub
#print(specs)

# Plot beachballs on map

with pygmt.config(PS_MEDIA='letter', MAP_FRAME_TYPE='plain', PS_PAGE_ORIENTATION='portrait', PROJ_LENGTH_UNIT='inch', FORMAT_GEO_MAP='D'):
    
    fig = pygmt.Figure()

    fig.basemap(projection=J, region=R, frame=True)
    fig.coast(borders=['1/0.75p,0/0/0','2/0.5p,0/0/0'],resolution='i',area_thresh='0/0/1',shorelines='0.5p,0/0/0',water='150/255/255',land='gray',lakes='150/255/255',map_scale=sbarinfo,F='l+ggray90')

    pygmt.makecpt(cmap='seis',series=[0,200,1])
    fig.meca(spec=specs, scale=ballsize,longitude=xvals, latitude=yvals, depth=depths, plot_longitude=xvals, plot_latitude=yvals, C=True, L='0.5p,black')
    fig.colorbar(frame='+lDepth')

    fig.savefig('mt_map.png')
