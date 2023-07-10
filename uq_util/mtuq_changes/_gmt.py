
import numpy as np
import subprocess

from mtuq.graphics._gmt import exists_gmt, gmt_not_found_warning, _parse_filetype, _safename
from mtuq.util import fullpath
from mtuq.graphics.uq._gmt import _parse_title, _savetxt, _parse_lune_array2



def _plot_depth_gmt(filename, depths, values, unc,
        magnitudes=None, lune_array=None,
        title='', xlabel='', ylabel='', fontsize=16.):

    # parse filenames
    filename, filetype = _parse_filetype(filename)

    ascii_file_1 = _safename('tmp_'+filename+'_ascii1.txt')
    ascii_file_2 = _safename('tmp_'+filename+'_ascii2.txt')
    ascii_file_3 = _safename('tmp_'+filename+'_ascii3.txt')


    # parase title and labels
    title, subtitle = _parse_title(title)

    xlabel = "'%s'" % xlabel
    ylabel = "'%s'" % ylabel


    data = np.column_stack((depths, values))
    minval, maxval, exp = _parse_limits(data[:,-1])

    # write values to be plotted as ASCII table
    _savetxt(ascii_file_1, data)

    if lune_array is not None:
        data2 = _parse_lune_array2(data[:,0], data[:,1], lune_array)
        data2[::,2] = unc
        _savetxt(ascii_file_2, data2)

    if magnitudes is not None:
        data3 = np.column_stack((data[:,0], data[:,1], magnitudes))
        _savetxt(ascii_file_3, data3, fmt='%e %e %.2f')

    # call bash script
    if exists_gmt():
        subprocess.call("%s %s %s %s %s %s %f %f %d %s %s %s %s" %
           (fullpath('../mtuq_supp/uq_util/mtuq_changes/plot_depth'),
            filename,
            filetype,
            ascii_file_1,
            ascii_file_2,
            ascii_file_3,
            minval,
            maxval,
            exp,
            title,
            subtitle,
            xlabel,
            ylabel,
            ),
            shell=True)
    else:
        gmt_not_found_warning(
            values_ascii)

def _parse_limits(values):

    masked = np.ma.array(values, mask=np.isnan(values))

    minval = masked.min()
    maxval = masked.max()
    exp = np.floor(np.log10(np.max(np.abs(masked))))

    if -1 <= exp <= 2:
        return minval, maxval, 0

    else:
        minval /= 10**exp
        maxval /= 10**exp
        masked /= 10**exp
        return minval, maxval, exp

