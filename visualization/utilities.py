import os
import sys

import numpy as np

# Utility functions from an old version of MTUQ used to gather station information, timeshifts, and cross correlations from the best fit moment tensor.

# Necessary to gather information for the misfit vs depth plotting tool, depth_pygmt.py, and used during your main run script for MTUQ.

def _get_depths(origins):
    depths = []
    for origin in origins:
        depths += [float(origin.depth_in_m)]
    return np.array(depths)

def _get_sources(sources, indices):
    return [sources.get(index) for index in indices]

def _min_dataarray(ds):
    values, indices = [], []
    for _i in range(ds.shape[-1]):
        sliced = ds[:,:,:,:,:,:,_i]
        values += [sliced.values.min()]
        indices += [int(sliced.values.argmin())]
    return np.array(values), indices

# Example useage:
#depths = _get_depths(origins)
#values, indices = _min_dataarray(results)
#best_sources = _get_sources(grid,indices)
#source_depth = {}
#for ii in range(len(best_sources)):
#    temp = best_sources[ii].as_dict()
#    temp['misfit'] = values[ii]
#    vr_temp = 100*(1-values[ii]/data_norm)
#    temp['vr'] = vr_temp
#    km = int(depths[ii]/1000)
#    source_depth[km] = temp
#    
#save_json(eid+'_depths_source.json',source_depth)
