#!/usr/bin/env python

import os
import numpy as np
import argparse

from mtuq import read, open_db, download_greens_tensors
from mtuq.event import Origin
from mtuq.graphics import plot_data_greens1, plot_beachball, plot_misfit_dc
from mtuq.grid import DoubleCoupleGridRegular
from mtuq.grid_search import grid_search
from mtuq.misfit import Misfit
from mtuq.process_data import ProcessData
from mtuq.util import fullpath, merge_dicts, save_json
from mtuq.util.cap import parse_station_codes, Trapezoid

def parse_args():

    parser = argparse.ArgumentParser(
    description="Input event info run MTUQ",
    formatter_class=argparse.RawTextHelpFormatter,
                                     )

    #parser.add_argument("-mdir",type=str,help="Main dir: -mdir /Users/felix/Documents/INVESTIGACION/2_FEB_JUL_2022/IRIS_WORKSHOP/MTUQ_INVERSIONS/FK_VS_1D_CUBIC_MESH/gs_mtuq")
    parser.add_argument("-event",type=str,help="event (event directory must be in main dir): -event 20140823183304000 ")
    parser.add_argument("-evla",type=str,help="Event latitude: -evla 64.68 ")
    parser.add_argument("-evlo",type=str,help="Event longitude: -evla -98.2 ")
    parser.add_argument("-evdp",type=float,help="Event depth in m: -evdp 10000.0 ")
    parser.add_argument("-mw",type=float,help="Event magnitude: -mw 4.7")
    parser.add_argument("-time", type=str,help="Earthquake origin time: -time 2014-08-25T16:19:03.00000Z")
    parser.add_argument("-np", type=int,help="Number of points per axis: -np 10")
    parser.add_argument("-fb",type=str,help="Frequency band for filtering data in seconds: -fb 10-100")


    return parser.parse_args()
                     

if __name__=='__main__':
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
        
    param = parse_args()
    print(param.event)
    mdir = os.getcwd()

    path_data=    fullpath('{}/{}/*.[zrt]'.format(mdir,param.event))
    path_weights= fullpath('{}/{}/weights.dat'.format(mdir,param.event))
    event_id=     '{}'.format(param.event)
    #model=        'ir'
    model=        'ak135'
    #db = open_db('{}/greens/ir'.format(param.mdir),format='FK')
    depth = int(param.evdp/1000)
    db = open_db('GFs/{}'.format(event_id),format="SPECFEM3D")


    #
    # Body and surface wave measurements will be made separately
    #

    freqs = param.fb.split('-')
    process_sw = ProcessData(
        filter_type='Bandpass',
        freq_min=1/float(freqs[1]),
        freq_max=1/float(freqs[0]),
        pick_type='taup',
        taup_model=model,
        #FK_database='{}/greens/ir'.format(param.mdir),
        window_type='surface_wave',
        window_length=250.,
        capuaf_file=path_weights,
        )

    #
    # For our objective function, we will use a sum of body and surface wave
    # contributions
    #

    misfit_sw = Misfit(
        norm='L2',
        time_shift_min=0.,
        time_shift_max=0.,
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
    #magnitudes=np.arange(param.mw-0.2,param.mw+0.2,0.1)
    magnitudes=param.mw
    
    grid = DoubleCoupleGridRegular(
        npts_per_axis=param.np,
        magnitudes= param.mw)
        #magnitudes= magnitudes.tolist())

    wavelet = Trapezoid(
        magnitude=param.mw)

    #
    # Origin time and location will be fixed. For an example in which they 
    # vary, see examples/GridSearch.DoubleCouple+Magnitude+Depth.py
    #
    # See also Dataset.get_origins(), which attempts to create Origin objects
    # from waveform metadata
    #

    origin = Origin({
        'time': '{}'.format(param.time),
        'latitude': param.evla,
        'longitude': param.evlo,
        'depth_in_m': param.evdp,
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
        data_sw = data.map(process_sw)

        print('Reading Greens functions...\n')
        greens = db.get_greens_tensors(stations,origin)
        #greens_tensors = db.get_greens_tensors(stations, origin)

        print('Processing Greens functions...\n')
        greens.convolve(wavelet)
        greens_sw = greens.map(process_sw)
        


    else:
        stations = None
        data_sw = None
        greens_sw = None


    stations = comm.bcast(stations, root=0)
    data_sw = comm.bcast(data_sw, root=0)
    greens_sw = comm.bcast(greens_sw, root=0)


    #
    # The main computational work starts now
    #

    if comm.rank==0:
        print('Evaluating surface wave misfit...\n')

    results_sw = grid_search(
        data_sw, greens_sw, misfit_sw, origin, grid)


    if comm.rank==0:

        results = results_sw

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

        for i in greens_sw:
            #print(i[1].stats.station)
            if i[1].stats.station == 'TPRV':
                i[1].write('TPRV_GF_took_from_mtuq.sac',format='SAC')
        #print(data_sw)

        #for i in data_sw:
            #print(i[1].stats.station)
        #    if i[1].stats.station == 'TPRV':
        #        i[1].write('TPRV_took_from_mtuq.sac',format='SAC')
        #print(data_bw)

        print(process_sw)

        plot_data_greens1(event_id+'DC_waveforms.png',
            data_sw, greens_sw, process_sw, 
            misfit_sw, stations, origin, best_source, lune_dict)


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
