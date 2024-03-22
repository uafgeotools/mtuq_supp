import numpy as np
from mtuq import open_ds
from mtuq_changes.depth_changes import plot_misfit_depth
from obspy.core.event import Origin
from itertools import chain
import warnings
warnings.filterwarnings('ignore')
result_All = open_ds('./input_data/output_result_All.hf5', format='HDF5')
print("readind data points")
df = result_All.reset_index()

depths = np.array(
     # depth in meters
    [25000., 30000., 35000., 40000.,                    
     45000., 50000., 55000., 60000., 65000.])
event_id= '20090407201255351'
catalog_origin = Origin({
    'time': '2009-04-07T20:12:55.000000Z',
    'latitude': 61.454200744628906,
    'longitude': -149.7427978515625,
    'depth_in_m': 33033.599853515625,
    })
npts = 64000
min_misfit = []
mt_parameters = []
unc = []
for i in range(len(depths)):
    min_misfit += [min(df[0][npts*i:npts*(i+1)])]
    df[0:npts] = df[npts*i:npts*(i+1)]
    result = df[0:npts]
    dl = result.loc[result[0].argmin()]
    mt_parameters += [dl.values[2:8]]
mt_parameters_all = np.array(mt_parameters)
    

origins = []
for depth in depths:
    origins += [catalog_origin.copy()]
    setattr(origins[-1], 'depth_in_m', depth)
import pandas as pd

# Specify the path to DAT file
file_path = "./input_data/unc_values.dat"

# Read the DAT file using pandas
data = pd.read_csv(file_path, delimiter="\t")

unc = list(chain.from_iterable(data.values.tolist()))
print("generating figures")
plot_misfit_depth(event_id + 'DC+Z_misfit_depth_tradeoffs.png', result_All, unc, mt_parameters_all, min_misfit, origins, show_tradeoffs=True, show_magnitudes=True, title=event_id)
print("finish")
