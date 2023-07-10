from mtuq import open_ds
from mtuq_changes import omega_changes
import warnings
warnings.filterwarnings('ignore')

print("reading data points")
df = open_ds('./input_data/20090407201255351DC_misfit.hf5', format='HDF5')
print("generating figures")
omega_changes.misfit_vs_omega('misfit_vs_omega.png', df)
omega_changes.plot_pdf('pdf.png', df, var= 50, nbins=40)
omega_changes.plot_cdf('cdf.png', df, var=50, nbins=40)
omega_changes.plot_rho_vs_V('rho_vs_v.png', df, var=50, nbins=40)
print("Finish")
