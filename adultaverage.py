# -*- coding: utf-8 -*-
"""
Examples of plots and calculations using the tmm package.
"""

from __future__ import division, print_function, absolute_import

from tmm_core import (coh_tmm, unpolarized_RT, ellips,
                       position_resolved, find_in_structure_with_inf)

from numpy import pi, linspace, inf, array
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
import numpy

try:
    import colorpy.illuminants
    import colorpy.colormodels
    from . import color
    colors_were_imported = True
except ImportError:
    # without colorpy, you can't run sample5(), but everything else is fine.
    colors_were_imported = False


# "5 * degree" is 5 degrees expressed in radians
# "1.2 / degree" is 1.2 radians expressed in degrees
degree = pi/180

#load experimental data so they can be plotted together with simulation

expdata = []
expdata = numpy.genfromtxt("adultaverage.txt",delimiter="\t")
expx = expdata[:,0]
expy = expdata[:,1]

def sample1():
    """
    Here's a thin non-absorbing layer, on top of a thick absorbing layer, with
    air on both sides. Plotting reflected intensity versus wavenumber, at two
    different incident angles.
    """

        
    # list of wavenumbers to plot in nm^-1
    
    ks=linspace(0.0025,.00125,num=800)
    
    # list of layer thicknesses in nm

    d0 = inf
    d1 = 92
    d2 = 90
    d3 = inf
    
    d_list = [d0, d1, d2, d1, d2, d1, d2, d1, d2, d1, d2, d1, d2, d1, d3]


    material_nk_data2 = array([[0.0025, 1.575+0.17j], # @400nm
                              [0.002, 1.5552+0.15j], # @500nm
                              [0.001667, 1.5444+0.14j], # @600nm
                              [0.001429, 1.53796+0.13j], # @700nm
                              [0.00125, 1.53375+0.12j]]) # @800nm
    material_nk_fn2 = interp1d(material_nk_data2[:,0].real,
                              material_nk_data2[:,1], kind='quadratic')

    material_nk_data1 = array([[0.0025, 1.705+0.06j], # @400nm
                              [0.002, 1.6852+0.04j], # @500nm
                              [0.001667, 1.6744+0.03j], # @600nm
                              [0.001429, 1.66796+0.02j], # @700nm
                              [0.00125, 1.66375+0.014j]]) # @800nm
 
    material_nk_fn1 = interp1d(material_nk_data1[:,0].real,
                              material_nk_data1[:,1], kind='quadratic')   
    
    # initialize lists of y-values to plot
    Rnorm=[] 
    R37=[]
    for k in ks:
		# For normal incidence, s and p polarizations are identical.
		# I arbitrarily decided to use 's'.
        n_list = [1, material_nk_fn1(k), material_nk_fn2(k), material_nk_fn1(k), material_nk_fn2(k), material_nk_fn1(k), material_nk_fn2(k), material_nk_fn1(k), material_nk_fn2(k), material_nk_fn1(k), material_nk_fn2(k), material_nk_fn1(k), material_nk_fn2(k), material_nk_fn1(k), 1.615+0.085j]
        Rnorm.append(coh_tmm('s',n_list, d_list, 0, 1/k)['R'])
        #R37.append(unpolarized_RT(n_list, d_list, 37*degree, 1/k)['R'])
    kcm = ks * 1e7 #ks in cm^-1 rather than nm^-1
    lamb = 1/ks
    plt.figure()
    plt.plot(lamb,Rnorm,'blue',expx,expy,'purple')
    plt.xlabel('Wavelength (/nm)')
    plt.ylabel('Fraction reflected')
    plt.title('Simulated reflection of unpolarized light at 0$^\circ$ incidence (blue), '
                'Experimental data (purple)')
    plt.axis([400, 800, 0, 0.25])
'''
# Save simulation data in a text file made out of columns delimited by a space
    import numpy
    datafile_path = "/users/Olimpia/Desktop/datafile.txt"
    datafile_id = open(datafile_path, 'w+')
    data = numpy.array([lamb, Rnorm])
    data = data.T
    numpy.savetxt(datafile_id, data, fmt = ['%.4f', '%.4f'])
    datafile_id.close()
'''

sample1()
plt.show()


