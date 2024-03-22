from matplotlib import pyplot
import numpy as np
from sklearn.metrics import auc
from scipy.interpolate import make_interp_spline


def _plot_omega_matplotlib_test(filename, omega, values,
    title=None, xlabel='Angular distance', ylabel=None, figsize=(6., 6.), fontsize=16.):

    pyplot.figure(figsize=figsize)
    pyplot.plot(omega, values, 'k.')
    pyplot.xlim([0., 180.])

    if title:
        pyplot.title(title, fontsize=fontsize)

    if xlabel:
         pyplot.xlabel(xlabel, fontsize=fontsize)

    if ylabel:
         pyplot.ylabel(ylabel, fontsize=fontsize)

    pyplot.savefig(filename)
    
def _plot_omega_matplotlib(filename, omega, values1, values2,
    title=None, xlabel='Angular distance', ylabel=None, figsize=(15., 6.), fontsize=16.):

    pyplot.figure(figsize=figsize)
    X_Y_Spline1 = make_interp_spline(omega, values1)
    X_Y_Spline2 = make_interp_spline(omega, values2)
    X_ = np.linspace(omega.min(), omega.max(), 1000)
    Y_1 = X_Y_Spline1(X_)
    Y_2 = X_Y_Spline2(X_)
    pyplot.plot(X_, Y_1, 'g-')
    pyplot.plot(X_, Y_2, 'b-')
    pyplot.xlim([0., 180.])

    if title:
        pyplot.title(title, fontsize=fontsize)

    if xlabel:
         pyplot.xlabel(xlabel, fontsize=fontsize)

    if ylabel:
         pyplot.ylabel(ylabel, fontsize=fontsize)

    pyplot.savefig(filename)
    
    
def _plot_rho_av(filename, x, y, xlabel='V', ylabel=(r"$\rho(V)$"), figsize=(6., 6.), fontsize=10.):
    pyplot.figure(figsize=figsize)
    p = [0,1]
    q = [0,1]
    pyplot.plot(x, y, 'r-')
    pyplot.plot(p, q, 'k--')
    pyplot.xlim([0., 1.])
    pyplot.ylim([0., 1.])
    pyplot.fill_between(
        x= x, 
        y1= y, 
        color= "k",
        alpha= 0.2)
    k = auc(x,y)
    #unc += [k]
    pyplot.text(0.7, 0.05, r"$\rho_{av} =$"+ str(k)[0:4])
    if xlabel:
         pyplot.xlabel(xlabel, fontsize=fontsize)

    if ylabel:
         pyplot.ylabel(ylabel, fontsize=fontsize)

    pyplot.savefig(filename)
