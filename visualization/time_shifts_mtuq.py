import pygmt
import sys
import os
import json

import numpy as np


if len(sys.argv) != 2:
    raise Exception('proper useage: python time_shifts_pygmt.py eid')
else:
    print('event id: ',sys.argv[1])

eid = sys.argv[1]

data_dir = '../../mtuq/'

model_list = data_dir + str(eid) + 'DC_best_origin.json'
origin_list = data_dir + str(eid) + 'DC_origin.json'
mt_list = data_dir + str(eid) + 'DC_mt.json'
station_list = data_dir + str(eid) + 'DC_stations.json'

# Key Commands

# srad_km: max distance to stations
# Jbscale: size of basemap (increase number to decrease size)
# color scale for time shifts for surface waves (dtminS, dtmaxS) and P (dtminP, drmaxP)

srad_km = 300
Jbscale = 5000000
dtminS = -10
dtmaxS = 10
dtminP = -2
dtmaxP = 2

# read in list of all available waveforms (station, network, lat, lon)

try:
    f = open(station_list,'r')
    data = json.load(f)
    f.close()

    staname = []
    lat = []
    lon = []
    
    for key in data.keys():
        if len(key) <= 7:
            staname.append(key[3:-1])
        elif len(key) >= 8:
            staname.append(key)
        #staname.append(data[key]['station'])
        lat.append(round(data[key]['latitude'],4))
        lon.append(round(data[key]['longitude'],4))

    staname = np.array(staname)
    lon = np.array(lon)
    lat = np.array(lat)
        
except:
    raise Exception('check to see if %s exists or not' % (station_list))

# read results information

try:
    f = open(model_list)
    data = json.load(f)
    f.close()

    edepcap = data['depth_in_m']
    smodel = 'tactmod'
except:
    raise Exception('check to see if %s exists or not' % (model_list))

try:
    f = open(origin_list)
    data = json.load(f)
    f.close()

    elat = data['0']['latitude']
    elon = data['0']['longitude']
    edep = data['0']['depth_in_m']
except:
    raise Exception('check to see if %s exists or not' % (origin_list))

try:
    f = open(mt_list)
    data = json.load(f)
    f.close()

    new_keys = ['mrr', 'mtt', 'mff', 'mrt', 'mrf', 'mtf']
    count = 0
    spec = {}

    for key in data.keys():
        nkey = new_keys[count]
        count += 1
        spec[nkey] = data[key]

    spec['exponent'] = 16
        
except:
    raise Exception('check to see if %s exists or not' % (mt_list))

# default values

# basement
bmin = -7
bmax = 0

dmin = 0
dmax = 160
depinc = 40
deptick = 10

itopocolor = 4    # 1-7 (default 2 or 7; 5-6=Jeff) # repeated var assignment

topocorr = 0

Jbscalei = 50000000

iportrait = 1

res = 100

iB = 8
stag = 'Alaska'

plot_unused_stations = 1

ifkmod = 5

# correction for high elevation stations
edepcap = edepcap - topocorr

iRUL = 1

if iportrait == 1:
    orient = 'P'
    rotangle = 0

# plotting specifications
fsize0 = '24'
fsize1 = '18'
fsize2 = '12'
fsize3 = '10'
fontno = '1'    # 1 or 4
tick = '0.3c'
fpen = '2p'
tpen = '2p'

stainfo = '-Si12p -W1p/0/0/0'
stainfo0 = '%s -G255/255/255' % (stainfo)
station_info = '-Sc5p -W0.5p/255/255/255 -G255/0/0'
station_info2 = '-Si10p -W0.5p/0/0/0 -G255/0/0'

textinfo = '-G255 -S1.5p'
textinfo2 = '-G0/0/0 -S2p,255/255/255'
textinfo3 = '-G0/0/0 -S2p,255/255/0'
inset_box_info = '-W3.0p,0/0/255 -N'

# topography/bathymetry
Dlen1 = 2
Dx = 0
Dy = Dlen1/2
Dwid = 0.20
Dscale_topo = '-D%s/%s/%s/%s' % (Dx,Dy,Dlen1,Dwid) #why are these the same
Bscale_topo = '-B1000f500:\"Elevation (m)\":'

Dscale_grav = '-D%s/%s/%s/%s' % (Dx,Dy,Dlen1,Dwid) #why are these the same
Bscale_grav = '-B20f10:\"Isostatic Gravity Anomaly (mgal)\": -E10p'

Dscale_base = '-D%s/%s/%s/%s' % (Dx,Dy,Dlen1,Dwid)#why are these the same
Bscale_base = '-B1f0.5:\"Basement surface (km)\": -Eb10p'

# CMT depth scale
Dlen2 = 2
Dx = 0
Dy = Dlen2/2
depinc = 40
Dscale_dep = '-D%s/%s/%s/%s' % (Dx,Dy,Dlen1,Dwid)
cmtsize = 0.4

# plot title
J_title = '-JM7i'  # -JM7i
R_title = '-R0/1/0/1'
fsize_title1 = 16
fsize_title2 = 14

title1 = 'Event ' + str(eid)
title2 = 'model %s, depth %.0f km (catalog %.1f km)' % (smodel,edepcap/1000,edep/1000)

itopo = 1

#==================================================================================
# LOOP OVER DIFFERENT SCALAR QUANTITIES FROM CAP

# 23 columns of CAP: weights in columns 3,7,11,15,19
# NOTE: TIME SHIFTS FOR Z AND R ARE FIXED -- (BOTH PNL AND SURF)

caplabs = ['ZRdt','Tdt','Zcc','Rcc','Tcc','ZRPnldt','ZPnlcc','RPnlcc','PVratio','PRratio','SVratio','SRratio','STratio']
capvals = ['x+l"Rayleigh(Z,R) time shift, s"','x+l"Love time shift, s"','x+l"Rayleigh(Z) CC"','x+l"Rayleigh(R) CC"','x+l"Love CC"','x+l"P(Z,R) time shift, s"','x+l"P(Z) CC"','x+l"P(R) CC"','x+l"ln(Aobs/Asyn), P(Z)"','x+l"ln(Aobs/Asyn), P(R)"','x+l"ln(Aobs/Asyn), Rayleigh(Z)"','x+l"ln(Aobs/Asyn), Rayleigh(R)"','x+l"ln(Aobs/Asyn), Love"']

capcols = [19,33,18,25,32,5,4,11,6,13,20,27,34] # column of data to plot
wcols = [16,30,16,23,30,2,2,9,2,9,16,23,30] # columns with weight factors
icc = [0,0,1,1,1,0,1,1,2,2,2,2,2]           # whether CC (1), tshift (0), data-syn ratio(2)
iwave = [0,2,1,0,2,3,4,3,4,3,1,0,2] # whether _rayleigh.json (R,Z) (0,1), _love.json (T) (2), or _bw.json (R,Z) (3,4) - TIMESHIFTS SHOULD BE SAME FOR (R,Z) COMPONENTS

# color limits
ccmin = 30
ccmax = 100
ccinc = 10

# Surface wave min-max shift
dtminsurf = dtminS
dtmaxsurf = dtmaxS

# P wave min-max shift
dtminpnl = dtminP
dtmaxpnl = dtmaxP

dtinc = 1
P_amp_ratio_min = -3
P_amp_ratio_max = -P_amp_ratio_min
S_amp_ratio_min = -2
S_amp_ratio_max = -S_amp_ratio_min
Pinc = 1
P_amp_ratio_inc = 1
S_amp_ratio_inc = 1

cmins = [dtminsurf,dtminsurf,ccmin,ccmin,ccmin,dtminpnl,ccmin,ccmin,P_amp_ratio_min,P_amp_ratio_min,S_amp_ratio_min,S_amp_ratio_min,S_amp_ratio_min]
cmaxs = [dtmaxsurf,dtmaxsurf,ccmax,ccmax,ccmax,dtmaxpnl,ccmax,ccmax,P_amp_ratio_max,P_amp_ratio_max,S_amp_ratio_max,S_amp_ratio_max,S_amp_ratio_max]
cticks = [dtinc,dtinc,ccinc,ccinc,ccinc,dtinc,ccinc,ccinc,P_amp_ratio_inc,P_amp_ratio_inc,S_amp_ratio_inc,S_amp_ratio_inc,S_amp_ratio_inc]

#==================================================================================
# BASEMAP OPTIONS

dy = srad_km / 100   # degrees
fac = np.cos(elat*np.pi/180.0)
xmin = elon - dy / fac
xmax = elon + dy / fac
ymin = elat - dy
ymax = elat + dy
xran = xmax - xmin
yran = ymax - ymin
xL = xmin + 0.9 * xran
yL = ymin + 0.1 * yran
#sbarinfo = 'g%.1f/%.1f+w100km+p1.5p,0/0/0,solid+f255/255/255' % (xL,yL)
sbarinfo = 'g%.1f/%.1f+w100+ab+u' % (xL,yL)
icities = 0
ititle = 1
iscalecap = 1
ibasementc = 0
islab_extent = 0

xtick1 = 2
ytick1 = 1
xtick2 = 0.5
ytick2 = xtick2
emax = 3000
emin = -emax
topocolor = ['globe','topo','relief','gray']
itopocolor = 3
ifmt = 0
cmtsize = 0.4

xcen = (xmin + xmax) / 2
ycen = (ymin + ymax) / 2

J = 'b%.1f/%.1f/%.1f/%.1f/1:%i' % (xcen,ycen,ymin,ymax,Jbscale)
R = '%.1f/%.1f/%.1f/%.1f' % (xmin,xmax,ymin,ymax)

# scale for plots
Bopts = ['WESN','Wesn','wesN','wEsn','weSn','WesN','wEsN','wESn','WeSn','WEsn','weSN','wESN','WESn','WeSN','WEsN','wesn']
#B0 = 'a%3.3ff%3.3fd:\" \":/a%3.3ff%3.3fd:\" \"::.\" \":' % (xtick1,xtick2,ytick1,ytick2)
B01 = 'a%3.3ff%3.3fd' % (xtick1,xtick2)
B02 = 'a%3.3ff%3.3fd' % (ytick1,ytick2)

B = [Bopts[iB],B01,B02]

# Load relief grid (will download 1st time is run, then load from local directory after)
grid = pygmt.datasets.load_earth_relief(resolution='01m',region=[-173,-128,51,72],use_srtm=True)
grad = pygmt.grdgradient(grid,radiance='p') #shading


# KEY COMMAND: min and max indices for plotting the maps listed above
xpmin = 1; xpmax = len(caplabs);  # for full set of figures
#xpmin = 1; xpmax = xpmin;    # for testing
pmax = 1

print('\n PLOTTING MAPS FROM %i TO %i\n' % (xpmin,xpmax))
# centered on epicenter
print('\n CUSTOM REGION CENTERED ON EPICENTER \n')


with pygmt.config(PS_MEDIA='letter', MAP_FRAME_TYPE='plain', MAP_TICK_LENGTH=str(tick), FORMAT_GEO_MAP='D', PROJ_LENGTH_UNIT='inch', PS_PAGE_ORIENTATION='portrait'):
    for xx in range(xpmax):

        clab = caplabs[xx]
        fname = '%s_%s_%s' % (eid,xx,clab)
        psfile = fname + '.eps'
        pdffile = fname + '.pdf'

        fig = pygmt.Figure()
        # LOOP OVER TWO DIFFERENT MAP REGIONS (FOR EACH PAGE)
        for pp in range(pmax):
            
            fig.basemap(projection=J, region=R, frame=B)
            if itopo:
                dc = (emax - emin) / 20
                pygmt.makecpt(cmap=topocolor[itopocolor], series=[emin,emax,dc], continuous=True)
                fig.grdimage(grid,cmap=True,shading=grad,dpi=res)

            fig.coast(borders=['1/2p,0/0/0','2/1p,0/0/0'],resolution='i',area_thresh='0/0/1',shorelines='0.75p,0/0/0',water='150/255/255',lakes='150/255/255', map_scale=sbarinfo)

            # If station has waveforms but is not used, then plot an open triangle
            if plot_unused_stations:
                fig.plot(x=lon, y=lat, pen='2p,0/0/0', style='i12p')
                fig.plot(x=lon, y=lat, pen='1p,255/255/255', style='i12p')

            # PLOTTING MTUQ RESULTS

            # make colorpoint file
            cmin = cmins[xx]
            cmax = cmaxs[xx]
            ctick = cticks[xx]
            ctick2 = ctick/2

            colormap = 'seis'

            if icc[xx] == 0:    # Time-shift
                dc = (cmax - cmin) / 9
                dkey = 'time_shift'

            elif icc[xx] == 1:  # Cross correlation
                dc = (cmax - cmin) / 12
                dkey = 'normalized_cc_max'

            elif icc[xx] == 2:  # log(data/syn)
                dc = (cmax - cmin) / 9
                dkey = 'log_amplitude_ratio'


                
            if iwave[xx] == 0:
                fname = data_dir + str(eid) + '_rayleigh.json'
                f = open(fname,'r')
                data = json.load(f)
                f.close()

                ckey = 'R'

            elif iwave[xx] == 1:
                fname = data_dir + str(eid) + '_rayleigh.json'
                f = open(fname,'r')
                data = json.load(f)
                f.close()

                ckey = 'Z'

            elif iwave[xx] == 2:
                fname = data_dir + str(eid) + '_love.json'
                f = open(fname,'r')
                data = json.load(f)
                f.close()

                ckey = 'T'

            elif iwave[xx] == 3:
                fname = data_dir + str(eid) + '_bw.json'
                f = open(fname,'r')
                data = json.load(f)
                f.close()

                ckey = 'R'

            elif iwave[xx] == 4:
                fname = data_dir + str(eid) + '_rayleigh.json'
                f = open(fname,'r')
                data = json.load(f)
                f.close()

                ckey = 'Z'
                                         
            sta_name = []
            fplot = []

            for key in data.keys():
                if data[key]:
                    if len(key) <= 7:
                        sta_name.append(key[3:-1])
                    elif len(key) >= 8:
                        sta_name.append(key)

                    try:
                        if icc[xx] == 1:
                            fplot.append(data[key][ckey][dkey]*100)
                        else:
                            fplot.append(data[key][ckey][dkey])
                    except:
                        continue

            
            pygmt.makecpt(cmap=colormap, series=[cmin,cmax,dc])

            for ii in range(len(fplot)):
                sta_out = str(sta_name[ii])
                ind = np.nonzero(sta_out == staname)[0]
                

                x = np.insert(lon[ind],0,elon)
                y = np.insert(lat[ind],0,elat)
                
                    
                # Plot colored ray paths
                fig.plot(x=x, y=y, pen='2p,0/0/0')
                fig.plot(x=x, y=y, pen='1p', cmap=True, zvalue=str(fplot[ii]))

                # For waveforms used in the inversion, plot over with a colored symbol
                fig.plot(x=lon[ind], y=lat[ind], pen='1p,0/0/0', cmap=True, zvalue=fplot[ii], color='+z', style='i12p')
                fig.text(x=lon[ind], y=lat[ind], text=sta_out, font='6p,Helvetica-Bold,black', angle=0, justify='RB', fill='white')

            # Colorbar
            pos = 'x%.1f/%.1f+w3.5/%.1f+h+e' % (Dx,Dy,Dwid)

            af_frame = 'a%.1fg%.1f' % (ctick, ctick2)
            xshift = 'a0.0'
            yshift = 'a5.3'
        
            fig.colorbar(xshift=xshift, yshift=yshift, position=pos,cmap=True, frame=[af_frame, capvals[xx]])

            # Make cpt for moment tensor
            pygmt.makecpt(cmap=colormap, series=[dmin,dmax,depinc], continuous=True)

            # plot moment tensor using the Mij values
            fig.meca(spec=spec, scale='0.4c', longitude=elon, latitude=elat, depth=edep/1000, plot_longitude=elon, plot_latitude=elat, C=True)

            # Plot title and subtitle
            fig.text(x=0, y=0, text=title1, region='0/1/0/1', projection='M7i', xshift='a0.0', yshift='a7.10', angle=0, justify='LM', font='16p,Helvetica-Bold,black', no_clip=True)
            fig.text(x=0, y=0, text=title2, region='0/1/0/1', projection='M7i', xshift='a0.0', yshift='a6.75', angle=0, justify='LM', font='14p,Helvetica-Bold,black', no_clip=True)

  

        fig.savefig(fname=pdffile)
        print('done with ',pdffile)



print('plotting finished')
