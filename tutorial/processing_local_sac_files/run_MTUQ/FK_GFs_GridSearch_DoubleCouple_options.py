#!/usr/bin/env python

import os
import numpy as np
import argparse

from mtuq import read, open_db, download_greens_tensors
from mtuq.event import Origin
from mtuq.graphics import plot_data_greens2, plot_beachball, plot_misfit_dc
from mtuq.grid import DoubleCoupleGridRegular
from mtuq.grid_search import grid_search
from mtuq.misfit import Misfit
from mtuq.process_data import ProcessData
from mtuq.util import fullpath, merge_dicts, save_json
from mtuq.util.cap import parse_station_codes, Trapezoid
#import pygmt

def launch_gs(event,evla,evlo,evdp,mw,time,nppa,fb_bw,fb_sw,w_bw,w_sw):
    #
    # Carries out grid search over 64,000 double couple moment tensors
    #
    # USAGE
    #   mpirun -n <NPROC> python GridSearch.DoubleCouple.py
    #
    # For a simpler example, see SerialGridSearch.DoubleCouple.py, 
    # which runs the same inversion in serial
    #

    #
    # We will investigate the source process of an Mw~4 earthquake using data
    # from a regional seismic array
    #

    """ Launch an MTUQ Double-Double grid-search the input paramenters are:
        launch_gs(event,evla,evlo,evdp,mw,time,nppa,fb_bw,fb_sw)
        -event: event directory (must be in main dir): 20140825161903
        -evla: event latitude.
        -evlo: event longitude.
        -evdp: event depth in meters.
        -mw: event magnitude.
        -time: origin time in ObsPy UTC format. Example: 2014-08-25T16:19:03.00000Z).
        -nppa: number of point per axis for defining how coarse or fine will be the grid-search.
        -fb_bw: frequency band for filtering body waves (in seconds). Example: 3-15. 
        -fb_sw: frequency band for filtering surface waves (in seconds). Example: 15-33.
        -w_bw: window lenght for body waves. 
        -w_sw: window lenght for surface waves. 

        An example for calling this method will be:
        launch_gs('20171201023244',30.734,57.39,6000.0,6.0,'2017-12-01T02:32:44.000000Z',30,'3-15','15-33',25,150)
    """   
    print(event)
    mdir = os.getcwd()

    path_data=    fullpath('{}/{}/*.[zrt]'.format(mdir,event))
    path_weights= fullpath('{}/{}/weights.dat'.format(mdir,event))
    event_id=     '{}'.format(event)
    model=        'ir'
    db = open_db('{}/greens/ir'.format(mdir),format='FK')

    #
    # Body and surface wave measurements will be made separately
    #

    freqs_bw = fb_bw.split('-')
    freqs_sw = fb_sw.split('-')

    process_bw = ProcessData(
        filter_type='Bandpass',
        freq_min= 1/float(freqs_bw[1]),
        freq_max= 1/float(freqs_bw[0]),
        pick_type='FK_metadata',
        FK_database='{}/greens/ir'.format(mdir),
        window_type='body_wave',
        window_length=w_bw,
        capuaf_file=path_weights,
        )
    
    process_sw = ProcessData(
        filter_type='Bandpass',
        freq_min=1/float(freqs_sw[1]),
        freq_max=1/float(freqs_sw[0]),
        pick_type='FK_metadata',
        FK_database='{}/greens/ir'.format(mdir),
        window_type='surface_wave',
        window_length=w_sw,
        capuaf_file=path_weights,
        )

    #
    # For our objective function, we will use a sum of body and surface wave
    # contributions
    #

    misfit_bw = Misfit(
        norm='L2',
        time_shift_min=-5.,
        time_shift_max=+5.,
        time_shift_groups=['ZR'],
        )
    
    misfit_sw = Misfit(
        norm='L2',
        time_shift_min=-15.,
        time_shift_max=+15.,
        time_shift_groups=['ZR','T'],
        )

    #
    # User-supplied weights control how much each station contributes to the
    # objective function
    #
    station_id_list = parse_station_codes(path_weights)
    #
    # Next, we specify the moment tensor grid and source-time function
    #
    print(mw)
    magnitudes=np.arange(mw-0.2,mw+0.2,0.1)
    
    grid = DoubleCoupleGridRegular(
        npts_per_axis=nppa,
        magnitudes=mw)
        #magnitudes= magnitudes.tolist())

    wavelet = Trapezoid(
        magnitude=mw)

    #
    # Origin time and location will be fixed. For an example in which they 
    # vary, see examples/GridSearch.DoubleCouple+Magnitude+Depth.py
    #
    # See also Dataset.get_origins(), which attempts to create Origin objects
    # from waveform metadata
    #

    origin = Origin({
        'time': '{}'.format(time),
        'latitude': evla,
        'longitude': evlo,
        'depth_in_m': evdp,
        })


    from mpi4py import MPI
    comm = MPI.COMM_WORLD

    #
    # The main I/O work starts now
    #

    if comm.rank==0:
        print('Reading data...\n')
        data = read(path_data, format='sac', 
            event_id=event_id,
            station_id_list=station_id_list,
            tags=['units:cm', 'type:velocity']) 


        data.sort_by_distance()
        stations = data.get_stations()

        print('Processing data...\n')
        data_bw = data.map(process_bw)
        data_sw = data.map(process_sw)

        #for i in data_bw:
        #    print(i[1].stats.station)
        #    if i[1].stats.station == 'K250':
        #        i[1].write('kalfz_took_from_mtuq.sac',format='SAC')
        #print(data_bw)

        print('Reading Greens functions...\n')
        greens = db.get_greens_tensors(stations,origin)

        print('Processing Greens functions...\n')
        greens.convolve(wavelet)
        greens_bw = greens.map(process_bw)
        greens_sw = greens.map(process_sw)


    else:
        stations = None
        data_bw = None
        data_sw = None
        greens_bw = None
        greens_sw = None


    stations = comm.bcast(stations, root=0)
    data_bw = comm.bcast(data_bw, root=0)
    data_sw = comm.bcast(data_sw, root=0)
    greens_bw = comm.bcast(greens_bw, root=0)
    greens_sw = comm.bcast(greens_sw, root=0)

    #
    # The main computational work starts now
    #

    if comm.rank==0:
        print('Evaluating body wave misfit...\n')
    
    results_bw = grid_search(
        data_bw, greens_bw, misfit_bw, origin, grid)
    
    if comm.rank==0:
        print('Evaluating surface wave misfit...\n')

    results_sw = grid_search(
        data_sw, greens_sw, misfit_sw, origin, grid)

    if comm.rank==0:

        results = results_bw + results_sw

        # array index corresponding to minimum misfit
        #idx = results.idxmin('source') #Old version
        idx = results.source_idxmin() #New version

        best_source = grid.get(idx)
        lune_dict = grid.get_dict(idx)
        mt_dict = grid.get(idx).as_dict()

        #
        # Generate figures and save results
        #

        print('Generating figures...\n')

        plot_data_greens2(event_id+'DC_waveforms.png',
            data_bw, data_sw, greens_bw, greens_sw, process_bw, process_sw, 
            misfit_bw, misfit_sw, stations, origin, best_source, lune_dict)

        plot_beachball(event_id+'DC_beachball.png',
            best_source, stations, origin)

        plot_misfit_dc(event_id+'DC_misfit.png', results)

        print('Saving results...\n')

        merged_dict = merge_dicts(lune_dict, mt_dict, origin,
            {'M0': best_source.moment(), 'Mw': best_source.magnitude()})

        # save best-fitting source
        save_json(event_id+'DC_solution.json', merged_dict)

        # save misfit surface
        results.save(event_id+'DC_misfit.nc')

        print('\nFinished\n')
                
if __name__=='__main__':
    launch_gs('20171201023244',30.734,57.39,6000.0,6.0,'2017-12-01T02:32:44.000000Z',30,'3-15','15-33',25,150)

    
