#!/usr/bin/env python

import os
import numpy as np

from mtuq import read, open_db, download_greens_tensors
from mtuq.event import Origin
from mtuq.graphics import plot_data_greens2, plot_misfit_force, plot_likelihood_force, plot_magnitude_tradeoffs_force
from mtuq.grid import ForceGridRegular, Grid
from mtuq.grid_search import grid_search
from mtuq.misfit import Misfit
from mtuq.process_data import ProcessData
from mtuq.util import fullpath, merge_dicts, save_json
from mtuq.util.cap import parse_station_codes, Trapezoid
from mtuq.misfit.waveform import estimate_sigma
from mtuq.grid.force import to_force



if __name__=='__main__':
    #
    # Example force search script
    #
    # USAGE
    #   mpirun -n <NPROC> python ForceGridSearch.py
    #
    # For a simpler example, see SerialGridSearch.DoubleCouple.py, 
    # which runs the same inversion in serial
    #


    #
    # We will investigate the source process of an Mw~4 earthquake using data
    # from a regional seismic array
    #

    path_data=    fullpath('data/examples/20090407201255351/*.[zrt]')
    path_weights= fullpath('data/examples/20090407201255351/weights.dat')
    event_id=     '20090407201255351'
    model=        'ak135'

    SF_path = '/Users/amanda/REPOSITORIES/mtuq_supp/greens_functions_libraries/specfem/cartesian_force/PROCESSED'
    #SF_path = './SPECFEM_GF'
    db = open_db(SF_path,format='SPECFEM3D',include_mt=False, include_force=True)

    #
    # Body and surface wave measurements will be made separately
    #

    process_bw = ProcessData(
        filter_type='Bandpass',
        freq_min= 0.1,
        freq_max= 0.333,
        pick_type='taup',
        taup_model=model,
        window_type='body_wave',
        window_length=15.,
        capuaf_file=path_weights,
        )

    process_sw = ProcessData(
        filter_type='Bandpass',
        freq_min=0.025,
        freq_max=0.0625,
        pick_type='taup',
        taup_model=model,
        window_type='surface_wave',
        window_length=150.,
        capuaf_file=path_weights,
        )


    #
    # For our objective function, we will use a sum of body and surface wave
    # contributions
    #

    misfit_bw = Misfit(
        norm='L2',
        time_shift_min=-2.,
        time_shift_max=+2.,
        time_shift_groups=['ZR'],
        )

    misfit_sw = Misfit(
        norm='L2',
        time_shift_min=-10.,
        time_shift_max=+10.,
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

    grid = ForceGridRegular(
        npts_per_axis=25,
        magnitudes_in_N=10.**np.arange(11.,12.,0.005))
    #grid = Grid(dims=('F0','phi','h'),coords=(9.89e11,338,-0.7986),callback=to_force)

    wavelet = Trapezoid(
        magnitude=4.5)


    #
    # Origin time and location will be fixed. For an example in which they 
    # vary, see examples/GridSearch.DoubleCouple+Magnitude+Depth.py
    #
    # See also Dataset.get_origins(), which attempts to create Origin objects
    # from waveform metadata
    #

    origin = Origin({
        'time': '2009-04-07T20:12:55.000000Z',
        'latitude': 61.454200744628906,
        'longitude': -149.7427978515625,
        'depth_in_m': 33033.599853515625,
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


        print('Reading Greens functions...\n')
        # If you want to see what the force inversion looks like using axiSEM instead of SPECFEM GFs, you can use the following line of code. Be sure to comment out Line 159.
        #greens = download_greens_tensors(stations, origin, model,
        #    include_mt=False, include_force=True)
        greens = db.get_greens_tensors(stations, origin)

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

        # `grid` index corresponding to minimum misfit
        idx = results.source_idxmin()

        best_source = grid.get(idx)
        force_dict = grid.get_dict(idx)
        print('Up-South-East Force Vector:')
        print(grid.get(idx).as_vector())
        #mt_dict = best_mt.as_dict()


        #
        # Generate figures and save results
        #

        print('Generating figures...\n')

        plot_data_greens2(event_id+'_GF_DC_waveforms.png',
            data_bw, data_sw, greens_bw, greens_sw, process_bw, process_sw, 
            misfit_bw, misfit_sw, stations, origin, best_source, force_dict)


        plot_misfit_force(event_id+'_GF_misfit_force.png',
                          results, title='L2 misfit')


        plot_magnitude_tradeoffs_force(event_id+'_force_tradeoffs.png',
                                       results, title='Magnitude tradeoffs')

        sigma = estimate_sigma(data_sw, greens_sw,
            best_source, misfit_sw.norm, ['Z','R'],
            misfit_sw.time_shift_min, misfit_sw.time_shift_max)

        plot_likelihood_force(event_id+'_GF_likelihood_force.png', 
            results, sigma**2, title='Maximum likelihoods')


        #print('Saving results...\n')

        # collect information about best-fitting source
        #merged_dict = merge_dicts(
        #    mt_dict,
        #    lune_dict,
        #    {'M0': best_mt.moment()},
        #    {'Mw': best_mt.magnitude()},
        #    origin,
        #    )

        # save best-fitting source
        #save_json(event_id+'DC_solution.json', merged_dict)


        # save misfit surface
        #results.save(event_id+'DC_misfit.nc')


        print('\nFinished\n')

