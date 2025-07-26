# IMPORT STATEMENTS:
import numpy as np
from matplotlib.pyplot import *

##################################################################################################################################################

def getObs(snr=10.0, ldist=10.0, **extras):
    """Build a dictionary of observational data. The data should 
    consist of photometry for a single galaxy.

    :param snr:
        The S/N to assign to the photometry.

    :param ldist:
        The luminosity distance to assume for translating absolute 
        magnitudes into apparent magnitudes.

    :returns obs:
        A dictionary of observational data to use in the fit.
    """
    from prospect.utils.obsutils import fix_obs
    import sedpy

    obs = {}

    filternames = ['jwst_f090w', 'jwst_f115w', 'jwst_f150w', 'jwst_f200w', 'jwst_f277w',
                   'jwst_f356w', 'jwst_f410m', 'jwst_f444w', 'jwst_f770w', 'jwst_f1800w']
    obs['filters'] = sedpy.observate.load_filters(filternames)

    ################################################# FIX THE UNITS OF THIS DATA
    M_AB = np.array([0.0148, 0.028, 0.0269, 0.322, 1.090, 1.334, 1.597, 1.725, 1.54, 1.96]) # These are currently in units of microJanskeys (fix)
    ################################################# FIX THE UNITS OF THIS DATA

    obs['maggies'] = 10**(-0.4*M_AB)
    obs['maggies_unc'] = (1./snr) * obs['maggies']

    obs['phot_wave'] = np.array([f.wave_effective for f in obs['filters']])

    obs['wavelength'] = None
    obs['spectrum'] = None
    obs['unc'] = None
    obs['mask'] = None

    obs = fix_obs(obs)

    return obs

##################################################################################################################################################

def getSps(zcontinuous=1, **extras):
    """
    :param zcontinuous:
        A value of 1 ensures that we use interpolation between SSPs to
        have a continuous metallicity parameter (`logzsol`)
        See python-FSPS documentation for details
    """

    # imf_upper_limit = 300
    # then imf_type = 0 (power law) (go through each of first 5 imf_type
    from prospect.sources import CSPSpecBasis
    sps = CSPSpecBasis(zcontinuous=zcontinuous)
    wspec = sps.wavelengths
    return sps, wspec

##################################################################################################################################################

def get_fixed_params():
    # IMFs = theta
    params = {'zred': 3.548,
              'mass': 1e8,
              'logzsol': 0.0,
              'dust2': 0.0,
              'sfh': 0,       # SSP
              'imf_type': 2,
              'dust_type': 0,
              'tage': 1.0,
             }
    
    return params

##################################################################################################################################################
# Defining functions for plotting

def plot_binned_imf(weights):

    imf_bins = np.array([[0.08, 0.5], [0.5, 2.0], [2.0, 4.0], [4.0, 5.0], [5.0, 120.0]])

    figure(figsize=(8,5))
    
    bin_labels = [f'{str(imf[0])}, {str(imf[1])}' for imf in imf_bins]
    imf_labels = [f'imf_{_}' for _ in range(1,len(imf_bins)+1)]
    
    N = 2
    
    n = 0
    x_arrs = np.empty([len(imf_bins), N])
    y_arrs = np.empty([len(imf_bins), N])
    for start, end, weight in zip(imf_bins[:,0],imf_bins[:,1], imf_weights):
        x_arrs[n] = np.linspace(start, end, N)
        y_arrs[n] = np.ones(N) * weight
    
        fill_between(x_arrs[n], y_arrs[n], 0, alpha=0.4, label=f'{imf_labels[n]}={weight:.2f}')
        plot(x_arrs[n], y_arrs[n])
        
        n+=1
    
    xscale('log')
    ylim([0,1])
    
    xlabel(r'Mass [$M_{\odot}$]')
    ylabel('dN/dm')
    
    legend()

    return






##################################################################################################################################################
# Defining functions for conversions



##################################################################################################################################################




##################################################################################################################################################




##################################################################################################################################################




##################################################################################################################################################
