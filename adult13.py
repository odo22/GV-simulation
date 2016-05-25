from __future__ import division, print_function, absolute_import

from tmm_core import (inc_tmm, unpolarized_RT, ellips,
                       position_resolved, find_in_structure_with_inf)

from numpy import pi, linspace, inf, array
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
import numpy
import math

try:
    import colorpy.illuminants
    import colorpy.colormodels
    from . import color
    colors_were_imported = True
except ImportError:
    colors_were_imported = False

#WARNING: the following statement is redundant
ks=linspace(0.0025,.00125,num=800)

# "5 * degree" is 5 degrees expressed in radians
# "1.2 / degree" is 1.2 radians expressed in degrees
degree = pi/180

#Load experimental spectra so they can be plotted together with simulation

expdata = []
expdata = numpy.genfromtxt("adult13.txt",delimiter="\t")
expx = expdata[:,0]
expy = expdata[:,1]

#Function to read off the layers' positions and gray values from ImageJ profile.
#Subtracts relative position values in the list to get absolute thickeness (in nm).

def loadGrayvalues(FIN):
    t_list = []
    gray_list = []
    f = open(FIN,'r')
    f.next()    #Skip the first line saying "X Y"
    f.next()    #Skip the first data point in 0.0
    for line in f:
        line = line.split()
        if len(t_list) == 0:
                thickness = float(line[0])
        else:
                for i in range (len(t_list),1):
                    thickness = t_list[i] - t_list[i-1] 
        t_list = t_list + [thickness]
        gray_list = gray_list + [float(line[1])]
    t_list = array(t_list)
    gray_list = array(gray_list)
    return t_list,gray_list

t_list,gray_list = loadGrayvalues("Plot Values 1.xls")

#Conversion from gray values to refractive index
maximum = 60 #numpy.amax(gray_list)
minimum = numpy.amin(gray_list)

for k in ks:
    
    N1R600 = 1.55 +0j #Real part of n1 at 600nm 1.45
    N1I600 = 0 + 0.001j #Imainary part of n1 at 600nm
    BR1 = 8800 # B Coefficient in Cauchy's equation (real)
    AR1 = N1R600 - BR1/(600**2) #Calculates Cauchy's A from given values of B and n at 600nm
    BI1 = 0 + 5188j #B coefficient in Cauchy's equation (imaginary)
    AI1 = N1I600 - BI1/(600**2) #Calculates Cauchy's A from given values of B and n at 600nm

    N2R600 = 1.72+0j #Real part of n1 at 600nm 1.65
    N2I600 = 0 + 0.06j #Imaginary part of n1 at 600nm
    BR2 = 23700 # B Coefficient in Cauchy's equation (real)
    AR2 = 1.648 #Calculates Cauchy's A from given values of B and n at 600nm
    BI2 = 270 #B coefficient in exponential equation(imaginary part)
    AI2 = 0.56#A coefficient in exponential equation (im part)

    n1 = AR1 + BR1*(k**2) + AI1 + BI1*(k**2)
    nmelanin = AR2 + BR2*(k**2) + (AI2*math.exp(-1/(270*k)))*1j 
    qpigm = 0.70
    n2 = n1*(1-qpigm) + nmelanin*qpigm
    
    conversion = (n2-n1)/(maximum-minimum)
    n_list_converted = n2 - conversion*(gray_list-minimum)
   
    #print (n1, k)

    '''
    wavel=1/k
    n1i = n1.imag
    n2i = n2.imag
    plt.figure(3)
    plt.plot(wavel, n1i, marker='o', ms = 10, alpha=1, color='b', label='Chitin layer')
    plt.plot(wavel, n2i, marker='o', ms = 10, alpha=1, color='k', label='Melanin layer')
    plt.xlabel('Wavelenght (/nm)')
    plt.ylabel('Refractive index')
    plt.title('Dispersion relations, imaginary part. Blue = Chitin layer. Black = Melanin layer')
    '''
    
#Adds bottom infinite layer at the end of the n list (average of ns)
    n_lists = numpy.append(n_list_converted, [n2])

#Adds air layer at the top of the n list
    air = 1.00029
    n_list = numpy.insert(n_lists, 0, air)

#Adds top and bottom inf layers to d list
t_lists = numpy.append(t_list, inf)
d_lists = numpy.insert(t_lists, 0, inf)
d_list = d_lists*0.95

# C_list: same lenght as t_list. First and last value must be 'i'.
a = []
a = ['c']*len(t_list)
b = numpy.append(a, 'i')
c_list = numpy.insert(b, 0, 'i')

def sample1():
        
    # list of wavenumbers to plot in nm^-1
    
    ks=linspace(0.0025,.00125,num=800)
    
    # initialize lists of y-values to plot
    Rnorm=[] 
    R37=[]
    for k in ks:
		# For normal incidence, s and p polarizations are identical.

		# I arbitrarily decided to use 's'.
              Rnorm.append(inc_tmm('s',n_list, d_list, c_list, 0, 1/k)['R'])
              #R37.append(unpolarized_RT(n_list, d_list, 37*degree, 1/k)['R'])
    kcm = ks * 1e7 #ks in cm^-1 rather than nm^-1
    lamb = 1/ks
    plt.figure(0)
    plt.plot(lamb,Rnorm,'blue',expx,expy,'purple')
    plt.xlabel('Wavelength (/nm)')
    plt.ylabel('Fraction reflected')
    plt.title('Simulated reflection of unpolarized light at 0$^\circ$ incidence (blue), '
                'Experimental data (purple)')
    plt.axis([400, 800, 0, 0.25])

# Save simulation data in a text file made out of columns delimited by a space
    import numpy
    datafile_path = "/users/Olimpia/Desktop/datafile1.txt"
    datafile_id = open(datafile_path, 'w+')
    data = numpy.array([Rnorm])
    data = data.T
    numpy.savetxt(datafile_id, data, fmt = ['%.4f'])
    datafile_id.close()

sample1()
#plt.show()


def profile():    
    position = []
    position = numpy.genfromtxt("Plot Values 1.xls",delimiter="\t")
    position = position[:,0]
    pos = numpy.delete(position, 0)
    poss = numpy.delete(pos,1)
    n_list_imaginary = n_list_converted.imag
    
    k= 0.002 #500nm
    plt.figure(1)
    plt.plot(poss, n_list_converted.imag, 'red')
    plt.xlabel('Distance (/nm)')
    plt.ylabel('Im{refractive index}')
    #plt.title('Simulated reflection of unpolarized light at 0$^\circ$ incidence (blue), '
                #'Experimental data (purple)')
    plt.axis([0, 2000, -0.10, 0.10])
profile()

plt.show()
