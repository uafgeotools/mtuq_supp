#!/usr/bin/env python

import os
import numpy as np

from mtuq import read, open_db, download_greens_tensors
from mtuq.event import Origin
import argparse
from mtuq.graphics import plot_data_greens1, plot_data_greens2, plot_beachball, plot_misfit_lune, plot_misfit_dc
from mtuq.grid import DoubleCoupleGridRegular, DeviatoricGridSemiregular, FullMomentTensorGridSemiregular
from mtuq.grid_search import grid_search
from mtuq.misfit import Misfit
from mtuq.process_data import ProcessData
from mtuq.util import fullpath, merge_dicts, save_json
from mtuq.util.cap import parse_station_codes, Trapezoid

class Header:
    def __init__(self,station,component,time_shift,cc):
        self.station = station
        self.component = component
        self.time_shift = time_shift
        self.cc = cc

def _getattr(trace, name, *args):
    if len(args)==1:
        if not hasattr(trace, 'attrs'):
            return args[0]
        else:
            return getattr(trace.attrs, name, args[0])
    elif len(args)==0:
        return getattr(trace.attrs, name)
    else:
        raise TypeError("Wrong number of arguments")
    
def get_headerinfo(data,greens,misfit,stations,origin,source):

    synthetics = misfit.collect_synthetics(data, greens.select(origin), source)

    header_info = []

    for _i in range(len(stations)):
        stream_dat = data[_i]
        stream_syn = synthetics[_i] 

        for dat in stream_dat:
            component = dat.stats.channel[-1].upper()
            try:
                syn = stream_syn.select(component=component)[0]
            except:
                warn('Missing component, skipping...')
                continue

            time_shift = 0.
            time_shift += _getattr(syn, 'time_shift', np.nan)
            time_shift += _getattr(dat, 'static_time_shift', 0)

            s = syn.data
            d = dat.data
            # display maximum cross-correlation coefficient
            Ns = np.dot(s,s)**0.5
            Nd = np.dot(d,d)**0.5

            if Ns*Nd > 0.:
                max_cc = np.correlate(s, d, 'valid').max()
                max_cc /= (Ns*Nd)
            else:
                max_cc = np.nan
                
            header_info.append(Header(stations[_i]['station'],component,np.round(time_shift,2),np.round(max_cc,2)))
            #print('{},{}: {} {}'.format(stations[_i]['station'],component,np.round(time_shift,2),np.round(max_cc,2)))

    return(header_info)

def wrap_up(ts_list,cc_list,station):
    total_ts = np.round(np.nansum(np.abs(ts_list))/(ts_list.size - np.count_nonzero(np.isnan(ts_list))),2)
    total_cc = np.round(np.nansum(cc_list)/(cc_list.size - np.count_nonzero(np.isnan(cc_list))),2)
            
    line = '{} {} {} {} {} {} {} {} {}'.format(station,ts_list[0],ts_list[1],ts_list[2],total_ts,cc_list[0],cc_list[1],cc_list[2],total_cc)
    return(line)

def write_headers(header_info,event_id):

    open_header_file=open('{}DC_header_info.txt'.format(event_id),'w')

    open_header_file.write('STATIONS  ts_Z ts_R ts_T abs_av_shift cc_Z cc_R cc_T av_cc\n')

    init_stat = header_info[0].station
    ts_list = np.array([np.nan,np.nan,np.nan])
    cc_list = np.array([np.nan,np.nan,np.nan])

    for i in range(len(header_info)):

        if header_info[i].station != init_stat:

            line = wrap_up(ts_list,cc_list,init_stat)
            open_header_file.write(line+'\n')
            #print(line)

            ts_list = np.array([np.nan,np.nan,np.nan])
            cc_list = np.array([np.nan,np.nan,np.nan])

            init_stat = header_info[i].station

        if header_info[i].station == init_stat:
            if header_info[i].component == 'Z':
                ts_list[0] = header_info[i].time_shift
                cc_list[0] = header_info[i].cc

            if header_info[i].component == 'R':
                ts_list[1] = header_info[i].time_shift
                cc_list[1] = header_info[i].cc

            if header_info[i].component == 'T':
                ts_list[2] = header_info[i].time_shift
                cc_list[2] = header_info[i].cc

            if i == len(header_info)-1:
                line = wrap_up(ts_list,cc_list,header_info[i].station)
                open_header_file.write(line)
                #print(line)

def save_results(filename,save_dir):
    print('\tSaving results in {}'.format(save_dir))
    if os.path.isdir(save_dir):
        pass
    else:
        mkdir = 'mkdir -p {}'.format(save_dir)
        os.system(mkdir)
    
    cp_file = 'mv {}* {}/'.format(filename,save_dir)
    os.system(cp_file)

#The common input parameters are
#Event ID: 20171201023244
#Green's function directory: greens/ir/ 
#Event latitude: 
#Event longitude:
#Event depth:
#Event magnitude range [mw_min,mw_max,mw_step]:
#Origin time:
#Point per axis:
#Band-pass filter:
#Band-pass filters [BW,SW] (for both BW and SW option):
#Window-lenght:
#Window-lengths [BW,SW] (for both BW and SW option):
#Maximum time shift:

## Double Couple moment tensor ###

def launch_MTUQ(gs_type,vel_model,waves,event,gf_dir,evla,evlo,evdp,mw,ot,ppa,filter,wl,ts,data_unit):

    '''
    Subroutine for launching different MTUQ moment tensor estimations.\n\n
    The expected parameters are:\n

    1. Type of grid-search (string): DC (double couple), DEV (deviatoric) and FMT(full moment tensor).\n
    2. Type of velocity model to use (string): 1D or 3D.\n
    3. Type of seismograms to use in the moment tensor estimation (integer): 
        1 (body waves), 2 (surface waves), 3 (both).\n
    4. Event_id (string): data directory. It follows the format yyyyMMddhhmmss. 
        y: year, M: month, d: day, m: month, s: seconds. E.g., 20171201023244.
        In addition to the SAC files, the data directory should contain the MTUQ parameters file (weights.dat).\n
    5. Green Functions (GFs) directory: path for finding the GFs.
        5.1 In case of using a 1D velocity model, the path must follow the format green/model/model_sourcedepth. 
            E.g., green/ir/ir_6.
        5.2 In case if using a 3D velocity model, the path must follow the format 3D/event_id/sourcedepth. 
            E.g., 3D/20171201023244/6.
        For the 3D GFs path it is only neccessary to speficify the main dir (e.g., 3D) since event_id and depth
        are found in the other options.\n
    6. Event latitude (float).\n
    7. Event longitude (float).\n
    8. Event depth (integer).\n
    9. Magnitude range (float array): range of magnitudes to explore [Mw_min,Mw_max,Mw_step].
       E.g., The array [5.8,6.2,0.1] will test the mangitudes: 5.8,5.9,6.0,6.1,6.2.\n
    10. Earthuake origin Time (string): it must follow the format of the ObsPy UTCDateTime object
       (year-month-dayThour:minute:second.milisecondsZ). E.g, 2017-12-01T02:32:44.000000Z.\n
    11. Number of point per axis (integer): controls how coarse or fine the grid search will be. \n
    12. Waveforms bandpass filter (string array): period corners for bandpass filtering the waveforms.
        12.1 In case of using body or surface waves (option 3 set to 1 or 2) a single pair of periods is required. 
            E.g., ['15-33'] for a bandpass filter between 15 and 33s.
        12.2 In case of using body and surface waves (option 3 set to 3) two pair of periods is required. 
            E.g., ['3-15','15-33']. The first pair is for body waves and the second for surface ones.\n
    13. Window length (integer array): lenght in seconds to the waveforms to use for estimating the moment tensor.
        13.1 In case of using body or surface waves a single length is required. 
            E.g., [150] for using a 150s seismogram.
        13.2 In  case of using body and surface waves two lengths are required. 
            E.g., [25,150]. The first value is for body waves and the second for surface ones.\n
    14. Maximum time-shift allowed (integer array).\n
        14.1 In case of using body or surface waves a single time-shift is required. 
            E.g., [15] for allowing a maximum shift of 25s.\n
        14.2 In case of using body and surface waves two time-shift are required. E.g., [5,15]. 
            The first value is for body waves and the second for surface ones.\n 
    15. Type of data unit (string): type of units of the observed seismograms. E.g., 'velocity', or 'displacement'.\n

    A full suite of parameters for running a double couple gridsearch, using body waves, 
    GFs calculated with a 1D velocity model should look like:\n
   
    gs_type = 'DC'
    vel_model = '1D'
    waves = 1 
    event = '20171201023244'
    gf_dir ='greens/ir'
    evla = 30.734
    evlo = 57.39
    evdp = 6000
    mw = [5.8,6.0,0.1]
    ot = '2017-12-01T02:32:44.000000Z'
    ppa = 30
    filter = ['15-33']
    wl = [150]
    ts = [15]
    data_unit = 'velocity'\n

    and the MTUQ grid-search is launched in this way:

    MTUQ_GridSearch_master.launch_MTUQ(gs_type,vel_model,waves,event,gf_dir,evla,evlo,evdp,mw,ot,ppa,filter,wl,ts,data_unit)
    '''

    #check_parameters(gs_type,vel_model,waves,event,gf_dir,evla,evlo,evdp,mw,ot,ppa,filter,wl,ts)

    mdir = os.getcwd()

    path_data=    fullpath('{}/{}/*.[zrt]'.format(mdir,event))
    path_weights= fullpath('{}/{}/weights.dat'.format(mdir,event))
    event_id=     '{}'.format(event)

   
    if vel_model == '1D':
        model= '{}'.format(gf_dir.split('/')[1])
        db = open_db('{}'.format(gf_dir),format='FK')
    else:
        model='ak135'
        depth = int(evdp/1000)
        db = open_db('{}/{}/{}'.format(gf_dir,event,depth),format="SPECFEM3D")

    if waves == 3:
        
        freqs_bw = filter[0].split('-')
        freqs_sw = filter[1].split('-')

        if vel_model == '1D':
            process_bw = ProcessData(
                filter_type='Bandpass',
                freq_min= 1/float(freqs_bw[1]),
                freq_max= 1/float(freqs_bw[0]),
                pick_type='FK_metadata',
                FK_database='{}/{}'.format(mdir,gf_dir),
                window_type='body_wave',
                window_length=wl[0],
                capuaf_file=path_weights,
                )
    
            process_sw = ProcessData(
                filter_type='Bandpass',
                freq_min=1/float(freqs_sw[1]),
                freq_max=1/float(freqs_sw[0]),
                pick_type='FK_metadata',
                FK_database='{}/{}'.format(mdir,gf_dir),
                window_type='surface_wave',
                window_length=wl[1],
                capuaf_file=path_weights,
                )
            
        elif vel_model == '3D':

            process_bw = ProcessData(
                filter_type='Bandpass',
                freq_min= 1/float(freqs_bw[1]),
                freq_max= 1/float(freqs_bw[0]),
                pick_type='taup',
                taup_model=model,
                #FK_database='{}/greens/ir'.format(mdir),
                window_type='body_wave',
                window_length=wl[0],
                capuaf_file=path_weights,
                )
    
            process_sw = ProcessData(
                filter_type='Bandpass',
                freq_min=1/float(freqs_sw[1]),
                freq_max=1/float(freqs_sw[0]),
                pick_type='taup',
                taup_model=model,
                #FK_database='{}/greens/ir'.format(param.mdir),
                window_type='surface_wave',
                window_length=wl[1],
                capuaf_file=path_weights,
                )

        misfit_bw = Misfit(
            norm='L2',
            time_shift_min=-1*ts[0],
            time_shift_max=ts[0],
            time_shift_groups=['ZR'],
            )
            
        misfit_sw = Misfit(
            norm='L2',
            time_shift_min=-1*ts[1],
            time_shift_max=ts[1],
            time_shift_groups=['ZR','T'],
            )


    elif waves == 1:
        freqs = filter[0].split('-')
        if vel_model == '1D':
            process_bw = ProcessData(
                filter_type='Bandpass',
                freq_min=1/float(freqs[1]),
                freq_max=1/float(freqs[0]),
                pick_type='FK_metadata',
                FK_database='{}/{}'.format(mdir,gf_dir),
                window_type='body_wave',
                window_length=wl[0],
                capuaf_file=path_weights,
                )
            
        elif vel_model == '3D':
            process_bw = ProcessData(
                filter_type='Bandpass',
                freq_min=1/float(freqs[1]),
                freq_max=1/float(freqs[0]),
                pick_type='taup',
                taup_model=model,
                window_type='body_wave',
                window_length=wl[0],
                capuaf_file=path_weights,
                )
        
        misfit_bw = Misfit(
            norm='L2',
            time_shift_min=-1*ts[0],
            time_shift_max=ts[0],
            time_shift_groups=['ZR'],
            )
        
    else :
        freqs = filter[0].split('-')
        if vel_model == '1D':
            process_sw = ProcessData(
                filter_type='Bandpass',
                freq_min=1/float(freqs[1]),
                freq_max=1/float(freqs[0]),
                pick_type='FK_metadata',
                FK_database='{}/{}'.format(mdir,gf_dir),
                window_type='surface_wave',
                window_length=wl[0],
                capuaf_file=path_weights,
                )
        elif vel_model == '3D':
            process_sw = ProcessData(
                filter_type='Bandpass',
                freq_min=1/float(freqs[1]),
                freq_max=1/float(freqs[0]),
                pick_type='taup',
                taup_model=model,
                window_type='surface_wave',
                window_length=wl[0],
                capuaf_file=path_weights,
                )
        misfit_sw = Misfit(
            norm='L2',
            time_shift_min=-1*ts[0],
            time_shift_max=ts[0],
            time_shift_groups=['ZR','T'],
            )
    
    station_id_list = parse_station_codes(path_weights)
    if mw[0] == mw[1]:
        magnitudes = np.array(mw[0])
    else:
        magnitudes=np.arange(mw[0],mw[1],mw[2])
    
    wavelet = Trapezoid(
        magnitude=np.median(magnitudes))
    
    if gs_type == 'DC':
        grid = DoubleCoupleGridRegular(
            npts_per_axis=ppa,
            magnitudes= magnitudes.tolist())
    elif gs_type == 'DEV':
        grid = DeviatoricGridSemiregular(
            npts_per_axis=ppa,
            magnitudes= magnitudes.tolist())
    else:
        grid = FullMomentTensorGridSemiregular(
            npts_per_axis=ppa,
            magnitudes= magnitudes.tolist()
            )
    
    origin = Origin({
        'time': '{}'.format(ot),
        'latitude': evla,
        'longitude': evlo,
        'depth_in_m': evdp,
        })
    
    from mpi4py import MPI
    comm = MPI.COMM_WORLD

    if comm.rank==0:
        print('Reading data...\n')
        data = read(path_data, format='sac', 
            event_id=event_id,
            station_id_list=station_id_list,
            tags=['units:{}'.format(data_unit), 'type:{}'.format(data_unit)])
        
        data.sort_by_distance()
        stations = data.get_stations()
        
        print('Reading Greens functions...\n')
        greens = db.get_greens_tensors(stations,origin)
        greens.convolve(wavelet)
    
        
        if  waves == 1:
            print('Processing data... \n')
            data_bw = data.map(process_bw)
            print('Processing Greens functions...\n')
            greens_bw = greens.map(process_bw)

        elif waves == 2:
            print('Processing data... \n')
            data_sw = data.map(process_sw)
            print('Processing Greens functions...\n')
            greens_sw = greens.map(process_sw)
        else:
            print('Processing data... \n')
            data_bw = data.map(process_bw)
            data_sw = data.map(process_sw)

            print('Processing Greens functions...\n')
            greens_bw = greens.map(process_bw)
            greens_sw = greens.map(process_sw)

    else:
        stations = None
        if waves == 1:
            data_bw = None
            greens_bw = None
        elif waves == 2:
            data_sw = None
            greens_sw = None
        else:
            data_bw = None
            data_sw = None
            greens_bw = None
            greens_sw = None

    if comm.rank==0:
        if waves == 1:
            print('Evaluating body wave misfit...\n')
            results_bw = grid_search(
            data_bw, greens_bw, misfit_bw, origin, grid)

            results = results_bw

        elif waves == 2:
            print('Evaluating surface wave misfit...\n')
            results_sw = grid_search(
            data_sw, greens_sw, misfit_sw, origin, grid)

            results = results_sw
        
        else:
            print('Evaluating body wave misfit...\n')
            results_bw = grid_search(
            data_bw, greens_bw, misfit_bw, origin, grid)

            print('Evaluating surface wave misfit...\n')
            results_sw = grid_search(
            data_sw, greens_sw, misfit_sw, origin, grid)

            results = results_bw + results_sw

    if comm.rank==0:

        idx = results.source_idxmin()
        best_source = grid.get(idx)
        lune_dict = grid.get_dict(idx)
        mt_dict = grid.get(idx).as_dict()

        print('Generating figures...\n')

        if waves == 1:
            plot_data_greens1(event_id+'{}_waveforms.png'.format(gs_type),
                data_bw, greens_bw, process_bw, 
                misfit_bw, stations, origin, best_source, lune_dict)
                
            header_info = get_headerinfo(data_bw,greens_bw,misfit_bw,stations,origin,best_source)
            write_headers(header_info,event_id)          
            
        elif waves == 2:
            plot_data_greens1(event_id+'{}_waveforms.png'.format(gs_type),
                data_sw, greens_sw, process_sw, 
                misfit_sw, stations, origin, best_source, lune_dict)
                
            header_info = get_headerinfo(data_sw,greens_sw,misfit_sw,stations,origin,best_source)
            write_headers(header_info,event_id)
            
        else:
            plot_data_greens2(event_id+'{}_waveforms.png'.format(gs_type),
                data_bw, data_sw, greens_bw, greens_sw, process_bw, process_sw, 
                misfit_bw, misfit_sw, stations, origin, best_source, lune_dict)
                
            #header_info = get_headerinfo(data_sw,greens_sw,misfit_sw,stations,origin,best_source)
            #write_headers(header_info,event_id)
        
        if gs_type == 'DC':
            plot_misfit_dc(event_id+'DC_misfit.png', results)

        else:
            plot_misfit_lune(event_id+'{}_misfit.png'.format(gs_type), results)
            plot_misfit_lune(event_id+'{}_misfit_mt.png'.format(gs_type), results, show_mt=True)
            plot_misfit_lune(event_id+'{}_misfit_tradeoff.png'.format(gs_type), results, show_tradeoffs=True)

        plot_beachball(event_id+'{}_beachball.png'.format(gs_type),
                    best_source, stations, origin)   
         
        print('Saving results...\n')

        merged_dict = merge_dicts(lune_dict, mt_dict, origin,
            {'M0': best_source.moment(), 'Mw': best_source.magnitude()})
        
        # save best-fitting source
        save_json(event_id+'{}_solution.json'.format(gs_type), merged_dict)

        # save misfit surface
        results.save(event_id+'{}_misfit.nc'.format(gs_type))

        if waves == 1:
            type_waves = 'BW'
            type_filter = filter[0]
        elif waves == 2:
            type_waves = 'SW'
            type_filter = filter[0]
        else:
            type_waves = 'BW_SW'
            type_filter = '{}_{}'.format(filter[0],filter[1])

        save_dir = 'SOLUTIONS/{}/{}/{}/{}_ppa/{}/{}/'.format(vel_model,event,gs_type,ppa,type_waves,type_filter)

        save_results('{}{}'.format(event,gs_type),save_dir)

        print('\nFinished\n')

def check_parameters(gs_type,vel_model,waves,event,gf_dir,evla,evlo,evdp,mw,ot,ppa,filter,wl,ts):
        
    print('I do not understand grid-search type {}'.format(gs_type))
    #Kill run
            
    print('For using body and surface waves you need to provide two frequency bands in filter')
    #Kill run
    
    print('I do not understand the type of waves to use. The options are:\n\t1 Body waves \n\t2 Surface waves \n\t3 Both.')
    #Kill run

    print('I do not understand the type of Grid Search')
    #Kill run

if __name__=='__main__':
    pass