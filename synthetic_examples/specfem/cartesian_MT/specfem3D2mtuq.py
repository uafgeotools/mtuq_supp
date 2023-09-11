import argparse
from obspy.core import read
from obspy import read
from obspy.geodetics.base import gps2dist_azimuth
from obspy.signal.rotate import rotate_ne_rt
import glob
import os
import numpy as np
from pandas import DataFrame
import pyproj

class Station:
    def __init__(self,name,lat,lon,network):
        self.name = name
        self.lat = lat
        self.lon = lon
        self.network = network

def clean():
    print('rm -r PROCESSED')
    os.system('rm -r PROCESSED')
    
def get_event_info(path):
    print('\tGathering information about Specfem3D simulation\n')

    #Read Par_file to figure out if the simulation was run in lat-long or UTM coordinates
    open_par = open('{}/DATA/Par_file'.format(path),'r')
    read_par = open_par.readlines()
    open_par.close()

    for line in read_par:
        if line.find('UTM_PROJECTION_ZONE') != -1:
            utm_zone = int(line.split('=')[1])
        if line.find('SUPPRESS_UTM_PROJECTION') != -1:
            if line.find('.true.') != -1:
                utm_on = True
            else:
                utm_on = False

    if utm_on:
        print('\tYour coordinates are in utm in zone {}. It will be converted to Latitude-Longitude\n'.format(utm_zone))

    open_cmt = open('{}/DATA/CMTSOLUTION'.format(path),'r')
    header_cmt = open_cmt.readlines()[0].split()
    open_cmt.close()

    #print(header_cmt)
    year = int(header_cmt[1])

    month = int(header_cmt[2])
    if month < 10:
        month = '0{}'.format(month)

    day = int(header_cmt[3]) 
    if day < 10:
        day = '0{}'.format(day)

    hour = int(header_cmt[4])
    if hour < 10:
        hour= '0{}'.format(hour) 

    min = int(header_cmt[5])
    if min < 10:
        min= '0{}'.format(min) 

    sec = float(header_cmt[6])
    if sec < 10:
        sec= '0{}'.format(sec)

    ev_id = '{}{}{}{}{}{}'.format(year,month,day,hour,min,int(sec))

    #2017-12-01T02:32:44.0
    time = '{}-{}-{}T{}:{}:{}'.format(year,month,day,hour,min,sec)
    evla = float(header_cmt[7])
    evlo = float(header_cmt[8])
    evdp = float(header_cmt[9])

    if utm_on:
        evlo,evla = cart_to_geo(evlo,evla,utm_zone)

    return(evla,evlo,evdp,time,ev_id,utm_on,utm_zone)

def txt2sac(path,ev_id,time):
    print('\tConverting plain text seismograms into SAC files\n')
    
    #evla,evlo,evdp,time,ev_id,utm_on = get_event_info(path)
    mkdir_output = 'mkdir PROCESSED'
    os.system(mkdir_output)

    DIR = '{}/OUTPUT_FILES/*sem*'.format(path)

    list_txt=glob.glob(DIR)

    #Extracting the tail of every sem file for detecting the type of simulation
    type_sim = list_txt[0].split('/')[-1].split('.')[-1]

    if type_sim == 'semv':
        print('\t********** I detected velocity simulations **********')
    elif type_sim == 'semd':
        print('\t********** I detected displacement simulations **********')


    for i in range(len(list_txt)):
        print('Pre-processing {}'.format(list_txt[i]))
        open_data=open(list_txt[i],'r')
        data=open_data.readlines()
        open_data.close()

        samples = len(data)
        time_s = np.round(np.abs(float(data[0].split()[0])) + np.abs(float(data[-1].split()[0])))
        sps = np.round(float(samples/time_s))
   
        data_series=[]
        exp=[]
        for j in range(len(data)):
            data_series.append(data[j].split()[1])
        write_data=open(list_txt[i]+'.ascii','w')
        name_stat=list_txt[i].split('/')[-1]
        split_name=name_stat.split('.')
        name_st=split_name[1]
        component=split_name[2]
        network=split_name[0]
        #name_stat=name_st+'_'+component+'_00'+'_'+network+'_R'
        name_stat=network+'_'+name_st+'_00'+'_'+component+'_R'
        #header='TIMESERIES '+name_stat+','+' 42000 samples, 100 sps, 2017-12-01T02:32:44.000000, SLIST, FLOAT,\n'
        header='TIMESERIES '+name_stat+','+' {} samples, {} sps, {}, SLIST, FLOAT,\n'.format(samples,sps,time)
        write_data.write(header)
        for j in range(len(data_series)):
            write_data.write(data_series[j]+'\n')
        write_data.close()
        st = read(list_txt[i]+'.ascii')

        new_name = rename(list_txt[i],ev_id)
        st.write(new_name, format='SAC')
        print (st)

        remove_ascii ='rm {}.ascii'.format(list_txt[i]) 
        print(remove_ascii)
        os.system(remove_ascii)

        mv_sac = 'mv {} PROCESSED/'.format(new_name)
        print(mv_sac)
        os.system(mv_sac)
        print('\n')

def rename(old_name,ev_id):
    #Mtt/OUTPUT_FILES/IR.JHBN.HXZ.semv
    #Mtt/OUTPUT_FILES/20171201023244.IR.JHBN..HH.z
    aux=old_name.split('/')[-1].split('.')
    main_dir = old_name.split('/')[0]
    network = aux[0]
    st = aux[1]
    comp = aux[2][-1]
    if comp == 'X':
        comp = 'E'
    elif comp == 'Y':
        comp = 'N'
    new_name = '{}/OUTPUT_FILES/{}.{}.{}..HH.{}'.format(main_dir,ev_id,network,st,comp.lower())
    return(new_name)

def geo_to_cart(lat,long,z):

    p = pyproj.Proj(proj='utm', zone=z, ellps='WGS84')
    x,y = p(long,lat)
    return(x,y)

def cart_to_geo(x,y,z):
    #https://www.engineeringtoolbox.com/utm-latitude-longitude-d_1370.html
    p = pyproj.Proj(proj='utm', zone=z, ellps='WGS84')
    lon,lat = p(x,y,inverse=True)
    return(lon,lat)

def grab_stations(path,utm_on,utm_zone):
    print('\tGathering information STATIONS file\n')
    #RST1 IR 4130709.407463263 910822.7389732231 0 0
    open_st = open('{}/DATA/STATIONS'.format(path),'r')
    read_st = open_st.readlines()
    station_list = []
    for line in read_st:
        st_name = line.split()[0]
        network = line.split()[1]
        lat = line.split()[2]
        lon = line.split()[3]

        if utm_on == True:
            lon,lat = cart_to_geo(lon,lat,utm_zone)
  
        station_list.append(Station(st_name,lat,lon,network))

    return(station_list)
    
def complete_header(process_path,stations,ev_id,evla,evlo,evdp):
    print('\tAdding values to SAC header\n')

    stream = read('{}/{}.*'.format(process_path,ev_id))
    
    for trace in stream:
        st_name = trace.stats.station

        #ADD EVLA,EVLO,EVDP,STLA,STLO 
        trace.stats.sac['evla'] = evla
        trace.stats.sac['evlo'] = evlo
        trace.stats.sac['evdp'] = evdp
        trace.stats.sac['lovrok'] = 1
        trace.stats.sac['lcalda'] = 1
        trace.stats.sac.khole = ''
        trace.stats.location = ''

        for st in stations:
            if st.name == st_name:
                stla = st.lat
                stlo = st.lon
                break
        trace.stats.sac['stla'] = stla
        trace.stats.sac['stlo'] = stlo

        trace.stats.sac.o = 0.0

        if trace.stats.sac.kcmpnm[-1] == 'X':
            trace.stats.sac.kcmpnm = 'HHE'
            trace.stats.channel = 'HHE'
            trace.stats.sac['cmpinc'] = 90 
            trace.stats.sac['cmpaz'] = 90

        if trace.stats.sac.kcmpnm[-1] == 'Y':
            trace.stats.sac.kcmpnm = 'HHN'
            trace.stats.channel = 'HHN'
            trace.stats.sac['cmpinc'] = 90 
            trace.stats.sac['cmpaz'] = 0 
        
        if trace.stats.sac.kcmpnm[-1] == 'Z':
            trace.stats.sac.kcmpnm = 'HHZ'
            trace.stats.channel = 'HHZ'           
            trace.stats.sac['cmpinc'] = 0 
            trace.stats.sac['cmpaz'] = 0

        name = '{}/{}.{}.{}..HH.{}'.format(process_path,ev_id,trace.stats.network,st_name,trace.stats.channel[-1].lower())
        print(name)
        trace.write(name, format='SAC')

def rotate(process_path,ev_id):
    print('\tRotating the radial and transverse\n')
    #20171201023244.IR.JHBN..HH.z
	
    list_st = glob.glob('{}/*.z'.format(process_path))
    for t in list_st:

        stZ = read(t)
        st_name = stZ[0].stats.station
        net = stZ[0].stats.network
        source_lat = stZ[0].stats.sac.stla
        source_lon = stZ[0].stats.sac.stlo
        st_lat = stZ[0].stats.sac.evla
        st_lon = stZ[0].stats.sac.evlo
        baz = gps2dist_azimuth(source_lat,source_lon,st_lat,st_lon)
        
        
        print('Rotating for {},{}'.format(process_path,st_name))
        stN = read('{}/*{}*.n'.format(process_path,st_name))
        stE = read('{}/*{}*.e'.format(process_path,st_name))

        rotated = rotate_ne_rt(stN[0].data,stE[0].data,baz[1])
        
        stN[0].data = rotated[0]
        stN[0].stats.sac.kcmpnm = '{}R'.format(stN[0].stats.sac.kcmpnm[0:-1])
        stN[0].stats.channel = '{}R'.format(stN[0].stats.sac.kcmpnm[0:-1])
        stN[0].stats.sac.cmpaz = baz[2]

        stE[0].data = rotated[1]
        stE[0].stats.sac.kcmpnm = '{}T'.format(stE[0].stats.sac.kcmpnm[0:-1])
        stE[0].stats.channel = '{}T'.format(stE[0].stats.sac.kcmpnm[0:-1])

        if baz[2] + 90 >= 360:
            stE[0].stats.sac.cmpaz = baz[2] + 90 - 360
        else :
            stE[0].stats.sac.cmpaz = baz[2] + 90

        #20171201023244.IR.JHBN..HH.z
        r_name = '{}/{}.{}.{}..HH.{}'.format(process_path,ev_id,net,st_name,stN[0].stats.channel[-1].lower())
        t_name = '{}/{}.{}.{}..HH.{}'.format(process_path,ev_id,net,st_name,stE[0].stats.channel[-1].lower())

        stN.write(r_name,format='SAC')
        stE.write(t_name,format='SAC')

def padd_zeros(event,ev_id,time):
    print('\tAdding {}s of zeros to {} seismograms\n'.format(time,event))
    #id_event = event.split('/')[-1]
    data = glob.glob('{}/{}.*'.format(event,ev_id))
    extra_time = 60

    for d in data:
        st = read(d)
        ns = round(extra_time/st[0].stats.delta)
        extra_samples = np.zeros(ns)
        st[0].data = np.concatenate((extra_samples,st[0].data,extra_samples))
        st[0].stats.starttime = st[0].stats.starttime - extra_time
        #st[0].data = 100*st[0].data 
        #st[0].stats.endtime = st[0].stats.endtime + extra_time
        st.write(d,format='SAC')

def scale_amplitude(event,scale):
    print('\tMultiplying by {}, amplitude of {} seismograms\n'.format(scale,event))
    id_event = event.split('/')[-1]

    if len(id_event) == 0:
        id_event = event.split('/')[-2]

    data = glob.glob('{}/{}*'.format(event,id_event))

    for d in data:
        st = read(d)
        st[0].data = st[0].data*scale
        st.write(d,format='SAC')

def write_weight_all(processed_dir,components):

    print('\tWriting weights.dat files for events in {} '.format(processed_dir))
    list_events = glob.glob('{}/*'.format(processed_dir))
     
    data = glob.glob('{}/*.r'.format(processed_dir))

    open_w_all = open('{}/weights.dat'.format(processed_dir),'w')

    for t in range(len(data)):
        name = data[t].split('/')[-1][0:-1]
        st = read(data[t])
        distance = gps2dist_azimuth(st[0].stats.sac['evla'],st[0].stats.sac['evlo'],st[0].stats.sac['stla'],st[0].stats.sac['stlo'])
        if t < len(data) - 1:
            line=' {} {} {}   0.0   0.0      0      0      0\n'.format(name,np.round(distance[0]/1000),components)
        else:
            line=' {} {} {}   0.0   0.0      0      0      0'.format(name,np.round(distance[0]/1000),components)
        open_w_all.write(line)
    open_w_all.close()

if __name__=='__main__':

    path='20171201023244'
    #1.Gathering information about Specfem3D simulation
    evla,evlo,evdp,time,ev_id,utm_on,utm_zone = get_event_info(path)
    print(evla,evlo,evdp,time,ev_id,utm_on,utm_zone)

    #2. Converting plain text seismograms into SAC files\
    txt2sac(path,ev_id,time)
    
    #3.Gathering information STATIONS file
    #utm_on = True
    #utm_zone = 38
    stations = grab_stations(path,utm_on,utm_zone)
    for st in stations:
        print('st:{} lat:{} lon:{}'.format(st.name,st.lat,st.lon,st.network))

    #4. Adding values to SAC header
    process_path = 'PROCESSED'
    #ev_id = '20171201023244'
    complete_header(process_path,stations,ev_id,evla,evlo,evdp)

    #5. Rotating the radial and transverse
    rotate(process_path,ev_id)

    #6. Adding zeroes to the trace onset
    extra_time = 60
    padd_zeros(process_path,ev_id,time)

    #7. Multiply seismograms by an scale factor
    scale = 100
    scale_amplitude(process_path,scale)

    #8. Writing weight file
    components = '0 0 1 1 1'
    write_weight_all(process_path,components)
    

    
