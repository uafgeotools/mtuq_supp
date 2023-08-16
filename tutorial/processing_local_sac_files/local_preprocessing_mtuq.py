from datetime import datetime
from obspy import UTCDateTime
from obspy.geodetics.base import gps2dist_azimuth
from obspy.signal.rotate import rotate_ne_rt
from obspy.signal.rotate import rotate2zne
import datetime
from obspy import read
import glob
import argparse
import os
import numpy as np
import io
from obspy import Trace
from obspy.io.sac import sacpz
from pathlib import Path
from obspy import read, read_inventory
#https://www.learnbyexample.org/python-classes-and-objects/

#Remove Ext_Catalog

class Earthquake:
    def __init__(self,or_time,lat,lon,depth,m,directory,traces):
        self.or_time = or_time
        self.lat = lat
        self.lon = lon
        self.depth = depth
        self.m = m
        self.directory = directory
        self.traces = traces

    def add_traces(self,traces):
        self.traces = traces 

class Dir_earthquake:
    def __init__(self,or_time,directory):
        self.or_time = or_time
        self.directory = directory

class Ext_catalog:
    def __init__(self,or_time,depth,stations,magnitude):
        self.or_time = or_time
        self.depth = depth
        self.stations = stations
        self.magnitude = magnitude

def add_list_events(file):
    open_file = open(file,'r')
    read_file = open_file.readlines()
    ev_all_list = []
    for r in read_file:
        year = int(r.split()[0])
        month = int(r.split()[1])
        day = int(r.split()[2])
        hour = int(r.split()[3].split(':')[0])
        min = int(r.split()[3].split(':')[1])
        sec = int(r.split()[3].split(':')[2])
        or_time = UTCDateTime(year,month,day,hour,min,sec)
        lat = float(r.split()[4])
        lon = float(r.split()[5])
        m = float(r.split()[6])
        depth = float(r.split()[7])

        ev_all_list.append(Earthquake(or_time,lat,lon,depth,m,'',[]))
    return(ev_all_list)

def add_directory_events(main_dir):
    dir_list = glob.glob('{}/20*'.format(main_dir))
    list_selected = []
    for d in dir_list:
        aux = d.split('/')
        year = int(aux[-1][0:4])
        month = int(aux[-1][4:6])
        day = int(aux[-1][6:8])
        hour = int(aux[-1][8:10])
        min = int(aux[-1][10:12])
        sec = int(aux[-1][12:14])
        dir_name = d.split('/')[-1]
        or_time = UTCDateTime(year,month,day,hour,min,sec)
        list_selected.append(Dir_earthquake(or_time,dir_name))
    return(list_selected)

def merge_lists(ev_file_list,ev_dir_list,events_dir):
    def_event_list = []

    for sel in ev_dir_list:
        flag = 0
        for all in ev_file_list:
            diff_time = all.or_time - sel.or_time
            if np.abs(diff_time) <= 60:
                print('Match found: {} and {}'.format(all.or_time,sel.or_time))
                flag = 1
                traces=get_traces(events_dir,sel.directory)
                if len(traces) == 0:
                    print('No seismograms for {}'.format(all.or_time))
                else:
                    all.add_traces(traces)
                    all.directory = sel.directory
                    def_event_list.append(all)
        if flag == 0:
            print('I cannot find a match for {}'.format(sel.or_time))    
    return(def_event_list)

def get_traces(main_dir,directory):
    list_traces = glob.glob('{}/{}/*.sac'.format(main_dir,directory))
    if len(list_traces) == 0:
        data = []
    else:
        data = read('{}/{}/*.sac'.format(main_dir,directory))
    return(data)

def quality_control(def_event_list,main_dir):
    print('Quality control for all the events in {}'.format(main_dir))

    for ev in def_event_list:
        #This block checks if the seismograms begin at the origin time, then cuts the seismograms, and addd header information
        print('\tReviewing {}/{} '.format(main_dir,ev.directory))
        no_st = []

        filter_list_events = []
        flag_qc = 0
        traces_passed = []

        for i in range(len(ev.traces)):
            
            approved,ban_st,message = review_data_length(ev.traces[i],ev.or_time,main_dir,ev.directory)
            if ban_st == 'null' :
                pass
            else:
                no_st.append(ban_st)
            if approved and ev.traces[i].stats.station not in no_st:
                pzfile_flag,message = check_poles(ev.traces[i])

                if pzfile_flag:
                    traces_passed.append(ev.traces[i])
                else:
                    print('\t\tTrace {} rejected.'.format(ev.traces[i]))
                    print('{}\n'.format(message))
                    flag_qc = 1
            else:
                flag_qc = 1
                print('\t\tTrace {} rejected.'.format(ev.traces[i]))
                print('{}\n'.format(message))
        
        if flag_qc == 0:
            print('\tAll traces in {}/{} passed the quality control'.format(main_dir,ev.directory))

        ev.traces = traces_passed
        filter_list_events.append(ev)
    
    return(filter_list_events)

def review_data_length(trace,or_time,main_dir,directory):
    number_of_traces = glob.glob('{}/{}/*{}*.sac'.format(main_dir,directory,trace.stats.station))
    ban_st = 'null'
    message ='ok'

    if len(number_of_traces) > 3:
        approved = False
        message = '\t\tThere are more than 3 sac files for *{}*.sac. Please merge your traces.'.format(trace.stats.station)
        ban_st = trace.stats.station

    elif len(number_of_traces) < 3 and trace.stats.channel[-1] != 'z' :
        approved = False
        message = '\t\tYour horizontal components for *{}*.sac are incomplete '.format(trace.stats.station)
        ban_st = trace.stats.station
    
    elif trace.stats.starttime <= or_time and or_time <= trace.stats.endtime:

        #The end time has to be the same for all the traces: here we check this
        all_same_st = read('{}/{}/*{}*.sac'.format(main_dir,directory,trace.stats.station))
        array_end_time = []
        array_start_time = []

        for t in all_same_st:
            array_end_time.append(t.stats.endtime)
            array_start_time.append(t.stats.starttime)

        array_end_time = np.array(array_end_time)
        min_end_time = np.min(array_end_time)

        flag_start_time = 0

        for time in array_start_time:
            if min_end_time < time :
                flag_start_time = 1
                approved = False
                ban_st = trace.stats.station
                message = '\t\t At least one of your components for {} has and endtime shorter than at least one of the starttime of the other components.'
        
        if flag_start_time == 0:
            approved = True
    else :
        print('\t\tstarttime: {}\n\t\tendtime: {}\n\t\t or_time: {}'.format(trace.stats.starttime,trace.stats.endtime,or_time))
        message = '\t\tThe trace onset is after or its end time is earlier that the origin time.\n'
        approved = False 
        ban_st = trace.stats.station
    return(approved,ban_st,message)

def check_poles(trace):

    trace.stats.station

    if os.path.isfile('pzfiles/{}_{}.pz'.format(trace.stats.station,trace.stats.channel)):
        flag_poles = True
        message =''
    else:
        flag_poles = False
        message = '\t\tFile pzfiles/{}_{}.pz does not exists'.format(trace.stats.station,trace.stats.channel)

    return(flag_poles,message)

def preprocessing(joint_event_list,processed_dir,events_dir):

    for ev in joint_event_list:
        traces_mixed = []

        for i in range(len(ev.traces)):
        
            print('\t\tCutting trace {}.\n\t\t Making first sample matching the origin time: {}'.format(ev.traces[i],ev.or_time))
            trace_cut = cut_seismo(ev.traces[i],ev.or_time,events_dir,ev.directory)
            print('\t\tModifying the SAC header and making demean and detrend for:\n\t\t{}\n '.format(trace_cut))
            trace_head = complete_header(trace_cut,ev.lat,ev.lon,ev.depth)
            traces_mixed.append(trace_head)
        
        print('\t\tChecking the directory {}/{}'.format(processed_dir,ev.directory))
        id_event = make_dir_event(ev.directory,processed_dir,events_dir,ev.or_time)
        print('\t\tSaving seismograms in {}/{} if the stations coordinates are found in station_list.txt'.format(processed_dir,id_event))
        save_seismo(traces_mixed,processed_dir,id_event)
        print('\n')
       
def cut_seismo(trace,or_time,main_dir,directory):

    all_same_st = read('{}/{}/*{}*.sac'.format(main_dir,directory,trace.stats.station))
    array_end_time = []
    array_start_time = []

    for t in all_same_st:
        array_end_time.append(t.stats.endtime)
        array_start_time.append(t.stats.starttime)

    array_end_time = np.array(array_end_time)
    min_end_time = np.min(array_end_time)

    max_end_time = np.max(array_end_time)
    add_time = max_end_time - trace.stats.endtime
    ns = round(add_time/trace.stats.delta)

    if add_time > 0:
        extra_samples = np.zeros(ns)
        trace.data = np.concatenate((trace.data,extra_samples))

    #print(trace.stats.station)
    #print(min_end_time)
    print('\t\t{}: t1={}  t2={} \n'.format(t.stats.station,or_time,trace.stats.endtime))
    #trace.trim(starttime=or_time,endtime=min_end_time)
    trace.trim(starttime=or_time,endtime=trace.stats.endtime)
    str_otime = str(or_time)
    year = str_otime.split('-')[0]
    month = str_otime.split('-')[1]
    day = str_otime.split('-')[2].split('T')[0]
    hour = str_otime.split('T')[1].split(':')[0]
    min = str_otime.split('T')[1].split(':')[1]
    sec = str_otime.split('T')[1].split(':')[2].split('Z')[0]
    msec = '0.0'

    trace.stats.sac['nzyear'] = int(year)
    trace.stats.sac['nzjday'] = datestdtojd(year+'-'+month+'-'+day)
    trace.stats.sac['nzhour'] = int(hour)
    trace.stats.sac['nzmin'] = int(min)
    trace.stats.sac['nzsec'] = int(float(sec))
    trace.stats.sac['nzmsec'] = int(float(msec))

    return(trace)

def datestdtojd (stddate):
    #https://rafatieppo.github.io/post/2018_12_01_juliandate/
    fmt='%Y-%m-%d'
    sdtdate = datetime.datetime.strptime(stddate, fmt)
    sdtdate = sdtdate.timetuple()
    jdate = sdtdate.tm_yday
    return(jdate)

def complete_header(trace_cut,lat,lon,depth):
    
    trace_cut.detrend("linear")
    trace_cut.detrend("demean")
    trace_cut.stats.sac['depmin'] = np.min(trace_cut.data)
    trace_cut.stats.sac['depmax'] = np.max(trace_cut.data)
    trace_cut.stats.sac['depmen'] = np.mean(trace_cut.data)

    #ADD EVLA,EVLO,EVDP,STLA,STLO 
    trace_cut.stats.sac['evla'] = lat
    trace_cut.stats.sac['evlo'] = lon
    trace_cut.stats.sac['evdp'] = depth
    trace_cut.stats.sac['lovrok'] = 1
    trace_cut.stats.sac['lcalda'] = 1
    trace_cut.stats.sac.o = 0.0

    if hasattr(trace_cut.stats, 'khole'):
        pass
    else:
        trace_cut.stats.sac.khole = ''
        trace_cut.stats.location = ''
    
    return(trace_cut)

def make_dir_event(ev_directory,destiny_dir,main_dir,or_time):

    str_otime = str(or_time)
    year = str_otime.split('-')[0]
    month = str_otime.split('-')[1]
    day = str_otime.split('-')[2].split('T')[0]
    hour = str_otime.split('T')[1].split(':')[0]
    min = str_otime.split('T')[1].split(':')[1]
    sec = str_otime.split('T')[1].split(':')[2].split('Z')[0].split('.')[0]

    new_dir_name = '{}{}{}{}{}{}'.format(year,month,day,hour,min,sec)
    
    destiny_path = '{}/{}'.format(destiny_dir,new_dir_name)
    origin_path = '{}/{}'.format(main_dir,ev_directory)

    if os.path.exists(destiny_path):
        pass
    else:
        mkdir_event = 'mkdir -p {}'.format(destiny_path)
        print('\t\t{}'.format(mkdir_event))
        os.system(mkdir_event)

    return(new_dir_name)

def save_seismo(traces_mixed,destiny_dir,id_event):
    st_info = get_stations(destiny_dir,id_event)
    #20171201023244.IR.CHTH..BH.z
    for trace in traces_mixed:
        stla,stlo,network,flag= index_st(st_info,trace)
        #print(stla,stlo,trace.stats.station)
        if flag == 1:
            trace.stats.sac['stla'] = stla
            trace.stats.sac['stlo'] = stlo
            trace.stats.sac['knetwk'] = network
            trace.stats.network = network
            az = gps2dist_azimuth(trace.stats.sac['evla'],trace.stats.sac['evlo'],stla,stlo)
            trace.stats.sac['az'] = az[1]
            trace.stats.sac['dist'] = az[0]/1000
            new_name = '{}.{}.{}.{}.{}.{}'.format(id_event,network,trace.stats.station,trace.stats.location,trace.stats.channel[0]+trace.stats.channel[1],trace.stats.channel[-1].lower())
            save_file = '{}/{}/{}'.format(destiny_dir,id_event,new_name)
            print('\t\tCreating {}'.format(save_file))
            trace.write(save_file,format='SAC')
        else:
            print('\t\tCannot find station coordinates for {}/{}/{}'.format(destiny_dir,id_event,trace.stats.station))

def get_stations(destiny_dir,ev_directory):
    file_name = 'station_list.txt'.format(destiny_dir,ev_directory)
    open_file = open(file_name,'r')
    read_STATIONS = open_file.readlines()
    open_file.close()
    return(read_STATIONS)

def index_st(st_info,trace):
    flag = 0
    st_name = trace.stats.station
    for i in st_info:
        if st_name in i :
            pos = st_info.index(i)
            flag = 1
            stla = float(st_info[pos].split()[1])
            stlo = float(st_info[pos].split()[2])
            network = st_info[pos].split()[4]
            break
    if flag == 0:
        stla = 'n'
        stlo = 'n'

    return(stla,stlo,network,flag)

def rotate_true_north(destiny_dir):
    print ('\t\tYou are going to rotate the seismograms according to the true north deviation.')
    print ('\t\tBe sure that your seismograms were not rotated already or make a backup.')  
    list_events = glob.glob('{}/*'.format(destiny_dir))

    #https://github.com/obspy/obspy/blob/a8cf88bfc28b7d06b88427ca626f67418922aa52/obspy/signal/rotate.py#L188-251
    #os.chdir('{}/{}'.format(destiny_dir,id_event))

    for ev in list_events:
        id_event = ev.split('/')[-1]
        list_st = glob.glob('{}/{}/*.z'.format(destiny_dir,id_event))
        non_rot = []
        for t in list_st:
            st_name =  t.split('/')[-1].split('.')[2]
            st_trace = glob.glob('{}/{}/*{}*'.format(destiny_dir,id_event,st_name))
        
            if len(st_trace) == 3 :
                print ('\t\tChecking for {}'.format(st_name))
                for s in st_trace:
                    if s.split('.')[-1] == 'n':
                        s1 = read(s)
                        tn = s
                    elif s.split('.')[-1] == 'e':
                        s2 = read(s)
                        te = s
                    elif s.split('.')[-1] == 'z':
                        s3 = read(s)
                        tz = s
                    else :
                        print ('\t\t We have an issue with {}. I cannot find a SAC file ending with n,e, or z. Check the name of your SAC file '.format(s))
                pz_file = glob.glob('pzfiles/{}_*N.pz'.format(st_name))
                open_pz = open(pz_file[0],'r')
                read_pz = open_pz.readlines()
                open_pz.close()
            
                for j in range(len(read_pz)):
                    if read_pz[j][0] == '*':
                        if read_pz[j].split()[1] == 'AZIMUTH':
                            print('\t\t {}'.format(read_pz[j].split()))
                            #print(read_pz[j])
                            #az = -1*int(read_pz[j].split()[3])
                            az = int(read_pz[j].split()[3])
                            break       
                if az == 0:
                    print ('\t\t {} station orientation is correct. No need of rotation'.format(st_name))
                else:
                    print ('\t\t {} station is {} off from the true north. A rotation will be applied'.format(st_name,az))
                    rotated = rotate2zne(s3[0].data,0,-90,s1[0].data,az,0,s2[0].data,az+90,0)  
                    s3[0].data = rotated[0]
                    s1[0].data = rotated[1]
                    s2[0].data = rotated[2]

                s3[0].stats.sac['cmpinc'] = 0 
                s3[0].stats.sac['cmpaz'] = 0 

                s1[0].stats.sac['cmpinc'] = 90 
                s1[0].stats.sac['cmpaz'] = 0 

                s2[0].stats.sac['cmpinc'] = 90 
                s2[0].stats.sac['cmpaz'] = 90

                s3.write(tz,format='SAC')
                s1.write(tn,format='SAC')
                s2.write(te,format='SAC')
            else:
                print('\t\tYour components are incomplete for {}. Rotation is not possible'.format(st_name))
                non_rot.append(st_name)

def get_pz_stations():
    list_pz = glob.glob('pzfiles/*.pz')
    pz_st = []

    for pz in list_pz:
        station_name = pz.split('/')[-1].split('_')[0]
        list_all = glob.glob('pzfiles/*{}*.pz'.format(station_name))
        if len(list_all) == 3:
            pz_st.append(station_name)
    return(pz_st)

def paz2resp(processed_dir):
    print('\n\t\tMaking RESP files from processed seismograms in {} \n'.format(processed_dir))
    list_ev = glob.glob('{}/*'.format(processed_dir))

    pz_list = []
    for ev in list_ev:
        list_seis = glob.glob(ev+'/*')
        for seis in list_seis:
            st = read(seis)
            st_name = st[0].stats.station
            channel = st[0].stats.channel
            location = st[0].stats.location
            network = st[0].stats.network
            pz_file = 'pzfiles/{}_{}.pz'.format(st_name,channel)
            if pz_file in pz_list:
                continue
            elif os.path.isfile(pz_file):
                pz_list.append(pz_file)
                poles,zeros,constant = grab_paz_data(pz_file)
                #print(poles)
                #print(zeros)
                #print(constant)
                #NHDN_BHZ.pz
                writeresp(poles,zeros,constant,st_name,channel,location,network)
            else:
                print('\n\t\tCannot find {} \n'.format(pz_file))
        
def grab_paz_data(paz_file):
    #Borrowed from https://github.com/obspy/obspy/blob/master/obspy/io/sac/sacpz.py

    poles = []
    zeros = []

    if isinstance(paz_file, (str, Path)):
        paz_file = open(paz_file, 'r')
        is_filename = True
    else:
        is_filename = False

    try:
        while True:
            line = paz_file.readline()
            if not line:
                break
            # lines starting with * are comments
            if line.startswith('*'):
                continue
            if line.find('ZEROS') != -1:
                a = line.split()
                noz = int(a[1])
                for _k in range(noz):
                    line = paz_file.readline()
                    a = line.split()
                    if line.find('POLES') != -1 or \
                       line.find('CONSTANT') != -1 or \
                       line.startswith('*') or not line:
                        while len(zeros) < noz:
                            zeros.append('0.000000E+00  0.000000E+00')
                        break
                    else:
                        zeros.append('{}  {}'.format(a[0], a[1]))

            if line.find('POLES') != -1:
                a = line.split()
                nop = int(a[1])
                for _k in range(nop):
                    line = paz_file.readline()
                    a = line.split()
                    if line.find('CONSTANT') != -1 or \
                       line.find('ZEROS') != -1 or \
                       line.startswith('*') or not line:
                        while len(poles) < nop:
                            poles.append('0.000000E+00  0.000000E+00')
                        break
                    else:
                        poles.append('{}  {}'.format(a[0], a[1]))
            if line.find('CONSTANT') != -1:
                a = line.split()
                # in the observatory this is the seismometer gain [muVolt/nm/s]
                # the A0_normalization_factor is hardcoded to 1.0
                constant = a[1]
    finally:
        if is_filename:
            paz_file.close()

    return(poles,zeros,constant)

def writeresp(poles,zeros,constant,st_name,channel,location,network):
    name_resp = 'RESP.{}.{}.{}.{}'.format(network,st_name,location,channel)
    print('\n\t\tMaking RESP_FILES/{} \n'.format(name_resp))

    if os.path.isdir('RESP_FILES'):
        pass
    else:
        mkdir_resp = 'mkdir RESP_FILES'
        print(mkdir_resp)
        os.system(mkdir_resp)

    name_resp = 'RESP.{}.{}.{}.{}'.format(network,st_name,location,channel)
    open_resp = open('RESP_FILES/{}'.format(name_resp),'w')
    open_resp.write('#               << IRIS SEED Reader, Release 4.6 >>\n#\n')
    open_resp.write('#               ======== CHANNEL RESPONSE DATA ========\n')
    open_resp.write('B050F03     Station:     {}\n'.format(st_name))
    open_resp.write('B050F16     Network:     {}\n'.format(network))
    open_resp.write('B052F03     Location:    {}\n'.format(location))
    open_resp.write('B052F04     Channel:     {}\n'.format(channel))
    open_resp.write('B052F22     Start date:  1900,1,00:00:00\n')
    open_resp.write('B052F23     End date:    No Ending Time\n')
    open_resp.write('#               =======================================\n')
    open_resp.write('#               +               +--------------------------------------------+               +\n')
    open_resp.write('#               +               |   Response (Poles & Zeros),  {} ch {}   |               +\n'.format(st_name,channel))
    open_resp.write('#               +               +--------------------------------------------+               +\n#\n')
    open_resp.write('B053F03     Transfer function type:                A [Laplace Transform (Rad/sec)])\n')
    open_resp.write('B053F04     Stage sequence number:                 1\n')
    open_resp.write('B053F05     Response in units lookup:              M - Displacement in Meters\n')
    open_resp.write('B053F06     Response out units lookup:             V - Volts\n')
    open_resp.write('B053F07     A0 normalization factor:               {}\n'.format(constant))
    open_resp.write('B053F08     Normalization frequency:               1.0\n')
    open_resp.write('B053F09     Number of zeroes:                      {}\n'.format(len(zeros)))
    open_resp.write('B053F14     Number of poles:                       {}\n'.format(len(poles)))
    open_resp.write('#               Complex zeroes:\n')
    open_resp.write('#                 i  real          imag          real_error    imag_error\n')
    cont_zeros = 0
    for zero in zeros:
        open_resp.write('B053F10-13    {} {}  {}  0.000000E+00  0.000000E+00\n'.format(cont_zeros,zero.split()[0],zero.split()[1]))
        cont_zeros = cont_zeros + 1
    open_resp.write('#               Complex poles:\n')
    open_resp.write('#                 i  real          imag          real_error    imag_error\n')
    cont_poles = 0
    for pole in poles:
        open_resp.write('B053F15-18    {} {}  {}  0.000000E+00  0.000000E+00\n'.format(cont_poles,pole.split()[0],pole.split()[1]))
        cont_poles = cont_poles + 1
    open_resp.write('#\n#               +                  +---------------------------------------+               +\n')
    open_resp.write('#               +                  |       Channel Gain,  {} ch {}      |               +\n'.format(st_name,channel))
    open_resp.write('#               +                  +---------------------------------------+               +\n')
    open_resp.write('#\nB058F03     Stage sequence number:                 1\n')
    open_resp.write('B058F04     Gain:                                  1.0\n')
    open_resp.write('B058F05     Frequency of gain:                     1.0 HZ\n')
    open_resp.write('B058F06     Number of calibrations:                0\n')
    open_resp.write('#\n#               +               +-------------------------------------------+               +\n')
    open_resp.write('#               +               |   Response (Coefficients),  {} ch {}   |               +\n'.format(st_name,channel))
    open_resp.write('#               +               +-------------------------------------------+               +\n')
    open_resp.write('#\nB054F03     Transfer function type:                D\n')
    open_resp.write('B054F04     Stage sequence number:                 2\n')
    open_resp.write('B054F05     Response in units lookup:              V - Volts\n')
    open_resp.write('B054F06     Response out units lookup:             COUNTS - Digital Counts\n')
    open_resp.write('B054F07     Number of numerators:                  0\n')
    open_resp.write('B054F10     Number of denominators:                0\n')
    open_resp.write('#\n##               +                      +------------------------------+               +\n')
    open_resp.write('#               +                      |   Decimation,  {} ch {}   |               +\n'.format(st_name,channel))
    open_resp.write('#               +                      +------------------------------+               +\n')
    open_resp.write('#\nB057F03     Stage sequence number:                 2\n')
    open_resp.write('B057F04     Input sample rate:                     5.120000E+03\n')
    open_resp.write('B057F05     Decimation factor:                     1\n')
    open_resp.write('B057F06     Decimation offset:                     0\n')
    open_resp.write('B057F07     Estimated delay (seconds):             0.000000E+00\n')
    open_resp.write('B057F08     Correction applied (seconds):          0.000000E+00\n')
    open_resp.write('#\n#               +                  +---------------------------------------+               +\n')
    open_resp.write('#               +                  |       Channel Gain,  {} ch {}      |               +\n'.format(st_name,channel))
    open_resp.write('#               +                  +---------------------------------------+               +\n')
    open_resp.write('#\nB058F03     Stage sequence number:                 2\n')
    open_resp.write('B058F04     Gain:                                  1.0\n')
    open_resp.write('B058F05     Frequency of gain:                     0.000000E+00 HZ\n')
    open_resp.write('B058F06     Number of calibrations:                0\n')
    open_resp.write('#\n#               +                  +---------------------------------------+               +\n')
    open_resp.write('#               +                  |   Channel Sensitivity,  {} ch {}   |\n'.format(st_name,channel))
    open_resp.write('#               +                  +---------------------------------------+               +\n')
    open_resp.write('#\nB058F03     Stage sequence number:                 0\n')
    open_resp.write('B058F04     Sensitivity:                           1.0\n')
    open_resp.write('B058F05     Frequency of sensitivity:              1.0 HZ\n')
    open_resp.write('B058F06     Number of calibrations:                0\n')
    open_resp.write('#')

def remove_instrumental_response_dir(processed_dir):
    print('\n\t\tRemoving instrumental response for events in {}'.format(processed_dir))
    list_events = glob.glob('{}/*'.format(processed_dir))

    for ev in list_events:
        id_event = ev.split('/')[-1]
        list_stream = glob.glob('{}/{}*'.format(ev,id_event))
        for seismogram in list_stream:
            print('\n\t\tRemoving instrumental response for {}'.format(seismogram))
            st = read(seismogram)
            tr = st[0]
            resp_path = 'RESP_FILES/RESP.{}.{}.{}.{}'.format(tr.stats.network,tr.stats.station,tr.stats.location,tr.stats.channel)
            inv = read_inventory('{}'.format(resp_path))
            tr.detrend("linear")
            tr.detrend("demean")
            pre_filt = [0.001, 0.005, 45, 50]

            tr.remove_response(inventory=inv,pre_filt=pre_filt,output="VEL",water_level=None,plot=False,taper=True,taper_fraction=0.05)
            tr.stats.sac.idep = 7
            tr.write(seismogram,format='SAC')

def rotate_radial_transverse_dir(processed_dir):
    print('\n\t\tRotate to radial and transverse events in {}'.format(processed_dir))
    list_events = glob.glob('{}/*'.format(processed_dir))

    for ev in list_events:
        id_event = ev.split('/')[-1]
        list_stream = glob.glob('{}/{}*.n'.format(ev,id_event))
        
        for seismogram in list_stream:
            print('\n\t\tRotating {},e'.format(seismogram))
            stN = read(seismogram)
            traceN = stN[0]
            aux = seismogram.split('/')[-1]
            chan_no_comp = aux.split('.')[4]
            name_e = '{}.{}.{}.{}.{}.e'.format(id_event,traceN.stats.network,traceN.stats.station,traceN.stats.location,traceN.stats.channel[0]+traceN.stats.channel[1])
            path_e = '{}/{}'.format(ev,name_e)

            if os.path.isfile(path_e):

                stE = read('{}'.format(path_e))

                st_name = stE[0].stats.station 
                net =  stE[0].stats.network
                location = stE[0].stats.location
                source_lat = stE[0].stats.sac.stla
                source_lon = stE[0].stats.sac.stlo
                st_lat = stE[0].stats.sac.evla
                st_lon = stE[0].stats.sac.evlo
                baz = gps2dist_azimuth(source_lat,source_lon,st_lat,st_lon)

                rotated = rotate_ne_rt(stN[0].data,stE[0].data,baz[1])

                stN[0].data = rotated[0]
                stN[0].stats.sac.kcmpnm = '{}R'.format(chan_no_comp)
                stN[0].stats.channel = '{}R'.format(chan_no_comp)
                stN[0].stats.sac.cmpaz = baz[2]

                stE[0].data = rotated[1]
                stE[0].stats.sac.kcmpnm = '{}T'.format(chan_no_comp)
                stE[0].stats.channel = '{}T'.format(chan_no_comp)

                if baz[2] + 90 >= 360:
                    stE[0].stats.sac.cmpaz = baz[2] + 90 - 360
                else :
                    stE[0].stats.sac.cmpaz = baz[2] + 90

                r_name = '{}/{}.{}.{}.{}.{}.r'.format(ev,aux.split('.')[0],net,st_name,location,chan_no_comp)
                t_name = '{}/{}.{}.{}.{}.{}.t'.format(ev,aux.split('.')[0],net,st_name,location,chan_no_comp)

                stN.write(r_name,format='SAC')
                stE.write(t_name,format='SAC')
            else:
                print('I cannot find {}'.format(path_e))

def padd_zeros(event,time):
    print('\n\t\tAdding {}s of zeros to {} seismograms'.format(time,event))
    id_event = event.split('/')[-1]
    data = glob.glob('{}/{}*'.format(event,id_event))
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
    print('\n\t\tMultiplying by {}, amplitude of {} seismograms'.format(scale,event))
    id_event = event.split('/')[-1]

    if len(id_event) == 0:
        id_event = event.split('/')[-2]

    data = glob.glob('{}/{}*'.format(event,id_event))

    for d in data:
        st = read(d)
        st[0].data = st[0].data*scale
        st.write(d,format='SAC')

def write_weight_all(processed_dir,components):

    print('\n\t\t Writing weights.dat files for events in {} '.format(processed_dir))
    list_events = glob.glob('{}/*'.format(processed_dir))
     
    for ev in list_events:

        print('\n\t\t Writing weights.dat for {} '.format(ev))
        data = glob.glob('{}/*.r'.format(ev))

        if len(data) > 0:
            open_w_all = open('{}/weights.dat'.format(ev),'w')

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
   
######
#mpirun -np 8 python FK_GFs_GridSearch.DoubleCouple_SW_options.py -event 20171201023244 -evla 30.734 -evlo 57.390 -evdp 5000 -mw 5.9 -time 2017-12-01T02:32:44.00000Z -np 30 -fb 10-100

if __name__=='__main__':   

    #1. Read all events info text file creating a list where each element is an object of the class Earthquake. The traces attribute will be empty at this point.
    file='events_list.txt'
    ev_all_list=add_list_events(file)

    #2. Read events info from RAW_DATA directory creating a list where each element is an object of the class Dir_earthquake
    events_dir = 'DATA/'
    list_selected = add_directory_events(events_dir)

    #3. We are going to compare the attribute or_time in the lists ev_all_list and list_selected for double checking the existence of SAC files for all events listed in in hojedk_events_IRSCloc_v0922.txt.
    #Once there is a match, then the SAC files are read and added to the traces attribute from each object of the list ev_all_list. 
    #If there is an event in hojedk_events_IRSCloc_v0922.txt that does not have data in ../RAW_DATA/ there will be a print message advertising the situation.
    #The list def_event_list is the final list for beginning the pre-processing. 
    joint_event_list=merge_lists(ev_all_list,list_selected,events_dir)

    #4. Quality control
    filter_list_events = quality_control(joint_event_list,events_dir)

    #5. Preprocessing
    #.Now let's to proceed with the pre-processing part 1: Change the name of the original files add certain values to the header SAC files and rotate .
    #process_seismograms_pt1(joint_event_list,destiny_dir,events_dir)
    processed_dir ='PROCESSED_DATA/'
    preprocessing(filter_list_events,processed_dir,events_dir)

    #6. Rotate to the true north
    #print('\t\tWorking on the new seismograms in {}/{} Rotating them to the true north in case to exist any deviation'.format(processed_dir,id_event))
    rotate_true_north(processed_dir)

    #7. Make RESP files for the processed seismograms
    paz2resp(processed_dir)

    #8. Remove Instrumental Response
    remove_instrumental_response_dir(processed_dir)

    #9. Rotate to RADIAL TRANSVERSE:
    rotate_radial_transverse_dir(processed_dir)

    #10. Add zero header marker: I did it on fix header

    #11.  Padd zeros
    event = '{}/20171201023244'.format(processed_dir)
    time = 60
    padd_zeros(event,time)

    #12. Scale Amplitude
    processed_dir = 'PROCESSED_DATA/'
    event = '{}/20171201023244/'.format(processed_dir)
    scale = 100
    scale_amplitude(event,scale)

    #12. Write Weight file
    processed_dir = 'PROCESSED_DATA/'
    components = '0 0 1 0 1'
    write_weight_all(processed_dir,components)


    #13. Run MTUQ



  