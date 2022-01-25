import pygmt
import sys
import json

import numpy as np

import faulthandler
faulthandler.enable()


if len(sys.argv) != 2:
    raise Exception('proper useage: python depth_pygmt.py eid')
else:
    eid = sys.argv[1]
    print('event id: ', eid)

data_dir = '/home/ammcpherson/REPOSITORIES/mtuq/'

# TEMPORARY
smodel='tactmod'

show_vr = False # plot variance reduction on split axes
imodel = True # draw bars at known layer interfaces
ballsize = '2.5' # controls default beachball size (meca)
min0 = 1.0e+19 # impossibly large misfit value

# Load in catalog origin information
file = data_dir + str(eid) + 'DC_origin.json'
try:
    f = open(file, 'r')
    data = json.load(f)
    f.close()

    # Grab catalog depth
    cdep = data['0']['depth_in_m']
    cdep = float(cdep) / 1000
except:
    raise Exception('check to see if %s exists or not' % (file))



# Load in best fit MT at each depth
file = data_dir + str(eid) + '_depths_source.json'
try:
    f = open(file, 'r')
    data = json.load(f)
    f.close()

except:
    raise Exception('check to see if %s exists or not' % (file))


# Unpack into more usable arrays
depth = []
vr = []
mag = [] #implement later
misfit = []
specs = {}
new_keys = ['mrr', 'mtt', 'mff', 'mrt', 'mrf', 'mtf']
count = 0

for key in data.keys():
    count = 0
    spec = {}
    depth.append(float(key))
    d = data[key]
    for k in d.keys():
        if k == 'vr':
            vr.append(d[k])
        elif k == 'misfit':
            misfit.append(d[k])
        else:
            nkey = new_keys[count]
            spec[nkey] = d[k]
            count += 1
    spec['exponent'] = 1.0
    specs[key] = spec
    
imin = np.argmin(misfit)

# Compute moment magnitude for each moment tensor
for key in specs.keys():
    temp = specs[key]

    M = np.array([[temp['mrr'],temp['mrt'],temp['mrf']],
                    [temp['mrt'],temp['mtt'],temp['mtf']],
                    [temp['mrf'],temp['mtf'],temp['mff']]], dtype=float)

    M0 = (np.tensordot(M,M)/2.)**0.5
    Mw = 2./3.*(np.log10(M0) - 9.1)
    mag.append(Mw)

# Compute a relative measure of error (lerr = ln(misfit/misfit_min))
lerr = []
for ii in range(len(misfit)):
    dum = np.log(misfit[ii]/misfit[imin])
    lerr.append(dum)

emax = np.amax(lerr)

# Compute the parabolic equation using best 3 points only
x1 = depth[imin-1]; y1 = lerr[imin-1]
x2 = depth[imin]; y2 = lerr[imin]
x3 = depth[imin+1]; y3 = lerr[imin+1]
denom = (x1-x2)*(x1-x3)*(x2-x3)

a = ((x3*(y2-y1)) + (x2*(y1-y3)) + (x1*(y3-y2))) / denom
b = (x3*x3*(y1-y2) + x2*x2*(y3-y1) + x1*x1*(y2-y3)) / denom
c = (x2*x3*(x2-x3)*y1 + x3*x1*(x3-x1)*y2 + x1*x2*(x1-x2)*y3) / denom

# Compute minima of parabola (at dy/dx = 0)
min_depth = -b/(2*a)
Ydepth = a*min_depth*min_depth + b*min_depth + c

# Setup of plotting parabola and uncertainty
xinc = 0.1
tmp = 1000000 #temporary variable (start with ridiculously large misfit value to find the uncertainty)
err_cent = 0.1 #to compute uncertainty (err_cent*100 percent confidence interval)
l = np.arange(start=depth[0],stop=depth[-1],step=0.1)

# Plotting stuff
# Size of plot
zoom = 3
magpsx = 2.0*zoom
magpsy = 1.8*zoom

# Projection, title, other options
J = 'X%s/%s' % (magpsx, magpsy)

xmin = float(min(depth)); xmax = float(max(depth));
ymin = -10; ymax = 100;
ymin2 = -emax/10; ymax2 = emax;

yplot = 1.426*ymin2 # for plotting catalog depth and best fit depth
yplot = -0.015

R = '%.1f/%.1f/%.1f/%.1f' % (xmin,xmax,ymin,ymax)
R2 = '%.1f/%.1f/%.1f/%.1f' % (xmin,xmax,ymin2,ymax2)


xlab = 'Depth, km'
ylab1 = 'Misfit relative to minimum'
ylab2 = 'VR (gray)'

xoffset = 'X3.5c'
yoffset = 'Y3.5c'

# Plotting begins

with pygmt.config(MAP_FRAME_TYPE='plain', PROJ_LENGTH_UNIT='inch', FONT_ANNOT='16'):#, PS_CHAR_ENCODING='Standard+'):
    fig = pygmt.Figure()
    fig.basemap(projection=J, region=R2, frame=['St','xaf+lDepth(km)'])
    #print('made basemap')
    if show_vr:
        fig.basemap(region=R,frame=['E','yaf+l"VR (gray)"'])
        # Plot line segments behind VR
        for ii in range(len(depth)-1):
            x = [depth[ii], depth[ii+1]]
            y = [vr[ii], vr[ii+1]]
            fig.plot(x=x,y=y,pen='0.5p,lightgray')
        
            # Plot VR
            fig.plot(x=depth,y=vr,style='p0.25c',color='lightgray')
            #print('plotted vr')
            
        # Switch axis
        fig.basemap(region=R2, frame=['W','yaf+llog(misfit/misfit_min)'])
    else:
        fig.basemap(region=R2, frame=['yaf+llog(misfit/misfit_min)'])

    # Plot catalog depth as red inverted triangle
    fig.plot(x=cdep, y=yplot+.01, pen='1p,black', style='i0.5c', color='red', no_clip=True)

    # Plot best fit depth as white inverted triangle
    fig.plot(x=min_depth, y=yplot+.01, pen='1p,black', style='i0.5c', color='white', no_clip=True)

    # Plot line segments in between beach balls
    for ii in range(len(depth)-1):
        try:
            x = [depth[ii], depth[ii+1]]
            y = [lerr[ii], lerr[ii+1]]
            fig.plot(x=x, y=y, pen='2p')
        except:
            print('skipping...')


    #Compute + plot parabola
    ycord = []
    for jj in range(len(l)):
        dum = a*l[jj]*l[jj] + b*l[jj] + c
        ycord.append(dum)

        if abs(dum-err_cent) < tmp:
            tmp = abs(dum-err_cent)
            unc = abs(l[jj]-min_depth)
            xcord2 = l[jj]
            ycord2 = dum
        
    fig.plot(x=l, y=ycord, pen='2p,--')

    title = '%s | Model %s | Best depth %.1f \\261 %.1f km\n' % (eid,smodel,min_depth,unc)
    #title = '%s | Model %s | Best depth %.1f ± %.1f km\n' % (eid,smodel,min_depth,unc)
    
    # Plot 10% confidence interval
    fig.plot(x=[min_depth,xcord2], y=[err_cent,err_cent], pen='2p')
    fig.plot(x=[min_depth,min_depth], y=[Ydepth,err_cent], pen='1p,--')

    # Plot beachballs
    for kk in range(len(depth)):
        key = str(int(depth[kk]))
        spec = specs[key]
        fig.meca(spec=spec, scale=ballsize,longitude=depth[kk],latitude=lerr[kk],depth=1,plot_longitude=depth[kk], plot_latitude=lerr[kk], G='red')
        mtxt = '%.2f' % (mag[kk])
        fig.text(x=depth[kk],y=lerr[kk]+0.015,text=mtxt,font='6p,Helvetica,black')


    # Plot title
    xtitle = depth[0]
    ytitle = ymax2 + (ymax2-ymin2)*0.05
    fig.text(x=xtitle,y=ytitle, text=title, font='16p,Helvetica,black',justify='LM',no_clip=True)
    

    fig.savefig('%s_depth.pdf'% (eid))
    
    


