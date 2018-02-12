##  USAGE:
##  python gaussians.py gaussians.cfg

import sys
import numpy as np


#########################################################################
##########################     SUBROUTINEs     ##########################
#########################################################################


def ParseConfigFile(cfg_file):
	global gauss_log_file,stick_out_file,gauss_out_file,fwhm,v_min,v_step,scale_factor,lineshape
	f = open(cfg_file)
	for line in f:
		# first remove comments
		if '#' in line:
			line, comment = line.split('#',1)
		if '=' in line:
			option, value = line.split('=',1)
			option = option.strip()
			value = value.strip()
			print "Option:", option, " Value:", value
			# check value
			if option.lower()=='gauss_log_file':
				gauss_log_file = value
			elif option.lower()=='stick_out_file':
				stick_out_file = value
			elif option.lower()=='gauss_out_file':
				gauss_out_file = value
			elif option.lower()=='fwhm':
				fwhm = float(value)
			elif option.lower()=='wavenumber_min':
				v_min = float(value)
			elif option.lower()=='wavenumber_step':
				v_step = float(value)
			elif option.lower()=='vibrational_scale_factor':
				scale_factor = float(value)
			elif option.lower()=='lineshape':
				lineshape = str(value)
			else :
				print "Option:", option, " is not recognized"
	
	f.close()

def ParseLogFile(log_file):

    f = open(log_file)
    freq_flag = "false"
    g09_flag = "false"
    g16_flag = "false"
    frequencies = []
    intensities = []
    for line in f:
        # first remove comments
        if '#' in line:
            line, comment = line.split('#',1)
        if 'Harmonic frequencies (cm**-1), IR intensities (KM/Mole)' in line:
            print "Frequency calculation recognized"
            freq_flag = "true"
        if 'Gaussian 09:' in line:
            print 'Gaussian 09 calculation recognized'
            g09_flag = 'true'
        if 'Gaussian 16:' in line:
            print 'Gaussian 16 calculation recognized'
            g16_flag = 'true'
        # Gaussian 16 flags
        if freq_flag == "true" and g16_flag == 'true' and "Frequencies --" in line:
            blah, info = line.split('--',1)
            data = info.split()
            for i in range(len(data)):
                frequencies.append(data[i])
        if freq_flag == "true" and g16_flag == 'true' and "IR Inten    --" in line:
            blah, info = line.split('--',1)
            data = info.split()
            for i in range(len(data)):
                intensities.append(data[i])
        # Gaussian 09 flags
        if freq_flag == "true" and g09_flag == 'true' and "Frequencies ---" in line:
            blah, info = line.split('---',1)
            data = info.split()
            for i in range(len(data)):
                frequencies.append(data[i])
        if freq_flag == "true" and g09_flag == 'true' and "IR Intensities ---" in line:
            blah, info = line.split('---',1)
            data = info.split()
            for i in range(len(data)):
                intensities.append(data[i])
    #print frequencies
    #print intensities
    print len(frequencies), " frequencies identified"
    stick = np.empty((len(frequencies),2),dtype=float)
    for i in range(len(frequencies)):
        stick[i,0] = float(frequencies[i])
        stick[i,1] = float(intensities[i])
    #print stick
    f.close()
    return stick


##########################################################################
##########################     MAIN PROGRAM     ##########################
##########################################################################

ParseConfigFile(sys.argv[1])

stick = ParseLogFile(gauss_log_file)
#stick = np.loadtxt(gauss_log_file)## load stick spectrum file
for i in range(stick.shape[0]):# Scale frequencies
    stick[i,0] *= scale_factor

v_end = int(stick[-1][0])
v = np.arange(v_min,v_end+1000,v_step)## Wavenumber axis
modes = np.zeros(len(v))## Intensities stick spectrum axis
peaks = np.zeros((len(v),1))## Intensities gaussian convolution axis
fwhm_steps = int(fwhm / v_step)## Number of steps on one side of gaussian



## Bin the sticks into the correct place on the 'v' axis based on the first column of the sticks array.

n_steps = stick.shape[0] # [Rows] Number of modes from freq calculation
fwhm_constants = 4.0*np.log(2.0)/(fwhm**2) ## the group of constants that make the FWHM easy to set rather than the standard deviation

if lineshape == "G":
    print 'Gaussinan lineshape selected'
elif lineshape == "L":
    print 'Lorentzian lineshape selected'
else:
    print 'Please enter valid \"lineshape\" value.'

for i in range(n_steps):
    index = int((stick[i][0] - v_min)/v_step)## index stick modes at nearest integer number steps from v_min.
    modes[index] += stick[i][1]## stick spectrum array
    #print 'Frequency %s to index %d...' %(stick[i][0],index)
    peaks[index] += modes[index]## convoluted spectrum array, convolution below happens exclusive of the stick frequency -- thus the addition of the stick intensity to the peaks array.

    for j in range(1,10*fwhm_steps + 1):
        if lineshape == "G":
            gaussian = modes[index] * np.exp(-fwhm_constants * (j*v_step) * (j*v_step))
            peaks[index + j][0] += gaussian
            peaks[index - j][0] += gaussian
        elif lineshape == "L":
            lorentz = modes[index] * 1.0 / ( ((j*v_step)*(j*v_step) / (fwhm * 2.0)) + 1.0 )
            peaks[index + j][0] += lorentz
            peaks[index - j][0] += lorentz


stick_spectrum = np.column_stack((v,modes))
np.savetxt(stick_out_file,stick_spectrum)

gauss_spectrum = np.column_stack((v,peaks))
np.savetxt(gauss_out_file,gauss_spectrum)
