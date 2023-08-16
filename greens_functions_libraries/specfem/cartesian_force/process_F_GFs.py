import argparse
from obspy.core import read
from obspy import read
from obspy.geodetics.base import gps2dist_azimuth
from obspy.signal.rotate import rotate_ne_rt
import glob
import os
import numpy as np

# Get station locations 
def index_st(st_info,s):
    aux = s.split('/')[-1]
    st_name = aux.split('.')[1]
    for i in st_info:
        if st_name in i :
            pos = st_info.index(i)
            break
    stla = float(st_info[pos].split()[2])
    stlo = float(st_info[pos].split()[3])
    return(stla,stlo)

# Complete header information
def get_stations_info():
    open_STATIONS = open('STATIONS','r')
    read_STATIONS = open_STATIONS.readlines()
    open_STATIONS.close()
    
    return(read_STATIONS)

# Get hypocenter information    
def get_event_info():
    open_info = open('FORCESOLUTION','r')
    read_info = open_info.readlines()
    evla = float(read_info[3].split()[1])
    evlo = float(read_info[4].split()[1])
    evdp = float(read_info[5].split()[1])
    return(evla,evlo,evdp)

# Complete the SAC header
def complete_header(file,sta,net,evla,evlo,evdp,stla,stlo):
    
    st = read(file)
    st[0].stats.sac.station = sta
    st[0].stats.sac.evla = evla
    st[0].stats.sac.evlo = evlo
    st[0].stats.sac.evdp = float(evdp) * 1000
    st[0].stats.sac.stla = stla
    st[0].stats.sac.stlo = stlo
    st[0].stats.sac.knetwk = net
    st[0].stats.sac.khole = ''
    st[0].stats.network = net
    st[0].stats.location = '00'
    st[0].stats.lcalda ='TRUE'
    
    st.write(file,format='SAC')

# Renaming and rewriting the ASCII files to SAC
def rename_make_sac(list_txt,otime):
       
    st_info = get_stations_info()
    evla,evlo,evdp = get_event_info()
        
    for file in list_txt:
        
        stla,stlo = index_st(st_info,file)
        # Grab bits of old name to cobble together new name
        Fcomp = file.split('/')[1].split('_')[0]
        old_name = file.split('/')[2]
        net = old_name.split('.')[0]
        sta = old_name.split('.')[1]
        Wcomp = old_name.split('.')[2][2]
        #print(Fcomp)
        #print(Wcomp)
        
        # Have to do a little editing to be able to read in Obspy
        open_data=open(file,'r')
        data=open_data.readlines()
        open_data.close()
        
        samples = len(data)
        time_s = np.round(np.abs(float(data[0].split()[0])) + np.abs(float(data[-1].split()[0])))
        sps = np.round(float(samples/time_s))
        
        # Make new data series in SLIST format - this is an intermediate step to write ASCII files into SAC files
        data_series = []
        for j in range(len(data)):
            data_series.append(data[j].split()[1])
            
        # Make header for SLIST format
        # Will need to make the time variable more general
        #time = '2009-04-07T20:12:55.000000'
        #time = '2017-12-30T11:43:16.278000'
        time = otime
        name_stat = sta + '_' + Wcomp + '_00_' + net + '_R'
        header='TIMESERIES '+name_stat+','+' {} samples, {} sps, {}, SLIST, FLOAT,\n'.format(samples,sps,time)
            
        write_data=open(file+'.ascii','w')
        write_data.write(header)
        
        for j in range(len(data_series)):
            write_data.write(data_series[j]+'\n')
            
        write_data.close()
        
        # Read in file
        st = read(file+'.ascii')
        
        new_name = '%s.%s..%s.F%s.sac' % (net,sta,Wcomp,Fcomp.lower())
        
        write_name = './PROCESSED/%s' % (new_name)
        
        # Write out in SAC format
        st.write(write_name,format='SAC')
        
        # Remove the .ascii files
        remove_ascii ='rm %s.ascii' % (file)
        os.system(remove_ascii)
        
        
        complete_header(write_name,sta,net,evla,evlo,evdp,stla,stlo)


def rotate(m,scale_factor):
    read_name = 'PROCESSED/*..Z.%s.sac' % (m)
    list_st = glob.glob(read_name)
    #list_st = glob.glob('PROCESSED/*.Z.*.sac'.format(m))
    for t in list_st:
        stZ = read(t)
        #print(stZ[0].stats)
        net = stZ[0].stats.network
        source_lat = stZ[0].stats.sac.evla
        source_lon = stZ[0].stats.sac.evlo
        st_lat = stZ[0].stats.sac.stla
        st_lon = stZ[0].stats.sac.stlo
        baz = gps2dist_azimuth(source_lat,source_lon,st_lat,st_lon)
        
        st_name = t.split('.')[1]
        print('Rotating for {}, {}'.format(m,st_name))
        stN = read('PROCESSED/{}.{}..N.{}.sac'.format(net,st_name,m))
        stE = read('PROCESSED/{}.{}..E.{}.sac'.format(net,st_name,m))

        rotated = rotate_ne_rt(stN[0].data,stE[0].data,baz[1])
    
        stN[0].data = rotated[0]
        stN[0].stats.kcmpnm = 'BHR'
        stN[0].stats.channel = 'BHR'
        stN[0].stats.cmpaz = baz[2]

        stE[0].data = rotated[1]
        stE[0].stats.kcmpnm = 'BHT'
        stE[0].stats.channel = 'BHT'

        if baz[2] + 90 >= 360:
            stE[0].stats.sac.cmpaz = baz[2] + 90 - 360
        else:
            stE[0].stats.sac.cmpaz = baz[2] + 90

        r_name = 'PROCESSED/{}.{}..R.{}.sac'.format(net,st_name,m)
        t_name = 'PROCESSED/{}.{}..T.{}.sac'.format(net,st_name,m)
        
        # scaling factor (based on force factor 1e15)
        sc = scale_factor
        sc = float(sc)
        print('Scaling %s by %s' % (r_name.split('/')[1],sc))
        stN[0].data = stN[0].data / sc
        print('Scaling %s by %s' % (t_name.split('/')[1],sc))
        stE[0].data = stE[0].data / sc
        #if m == 'Fe':
        print('Scaling %s by %s' % (t.split('/')[1],sc))
        stZ[0].data = stZ[0].data / sc
        
        stN.write(r_name,format='SAC')
        stE.write(t_name,format='SAC')
        stZ.write(t,format='SAC')


# Actually running all the preprocessing steps
DIR = 'F_GFs/*/'
otime_2009 = '2009-04-07T20:12:55.000000'
scale_factor_2009 = 1e17
otime_2017 = '2017-12-30T11:43:16.278000'
scale_factor_2017 = 1e15

list_txt=glob.glob(DIR+'*semd')
rename_make_sac(list_txt,otime_2009)

M = ['Fe','Fn','Fz']
for m in M:
    rotate(m,scale_factor_2009)

