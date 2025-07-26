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
# Defining all of the conversion functions to be used for priors

# def sf_start_from_fstart(**sfh):
#     pset = getNamespace(**sfh)
#     return (pset.tage * pset.fstart)

def sf_trunc_from_ftrunc(tage,ftrunc):
    if ftrunc == 0.:
        return 0.0
    end_trunc = tage * ftrunc
    return tage - end_trunc

def tau_from_logtau(logtau):
    return 10**logtau
    
##################################################################################################################################################

# def getNamespace(**params):
#     from argparse import Namespace
#     pset = Namespace(**params)
#     return pset

##################################################################################################################################################

def getBreakIndices(wspec):
    
    blue, red = 3620, 4000 # Indices used in De Graaff, 2025
    
    all_indices = np.zeros(len(wspec), dtype=bool)
    blue_lower, blue_upper, red_lower, red_upper = False, False, False, False
    
    for i, wave in enumerate(wspec):
        if wave > blue and blue_lower==False:
            all_indices[i], blue_lower = True, True
        elif wave > blue+100 and blue_upper==False:
            all_indices[i-1], blue_upper = True, True
        elif wave > red and red_lower==False:
            all_indices[i], red_lower = True, True
        elif wave > red+100 and red_upper==False:
            all_indices[i-1], red_upper = True, True
    
    true_indices = [index for index, value in enumerate(all_indices) if value]
    
    return true_indices

##################################################################################################################################################

# Make sure that the parameter(s) for logzsol (and tage) are maximized as according to the original parameter fits from earlier work
def get_fixed_params():
    # tage, logtau, fburst, fage_burst, const, ftrunc = theta
    params = {'zred': 3.548,
              'mass': 1e8,
              'logzsol': 0.0,
              'dust2': 0.0,
              'sfh': 1,
              'imf_type': 2,
              'dust_type': 0,
              'tage': 1.0,
              'logtau':-1.0,
              'tau': 0.1,
              'fburst': 0.0,
              'fage_burst': 0.0,
              'tburst': 0.0,
              'const': 0.0,
              'fstart':0.0,
              'ftrunc':0.0,
              'sf_start': 0.0,
              'sf_trunc': 0.0}
    return params

##################################################################################################################################################

def compute_break_strength(spec):
    """This function will compute the balmer break strengh for prospecter spectra
    by averaging the flux from 4000-4100 and from 3620-3720 Angstroms (rest-frame)
    then dividing the redder average flux by the bluer (De Graaff, 2025).
    """
    blue_l, blue_u, red_l, red_u = 382, 492, 804, 914
    return spec[red_l:red_u].mean() / spec[blue_l:blue_u].mean()

##################################################################################################################################################

def init_prospect_generation():
    params = get_fixed_params()
    sps, wspec = getSps()
    obs = getObs()

    return params, obs, sps, wspec

##################################################################################################################################################

def emcee_model(theta, params, obs, sps):
    from prospect.models.transforms import tburst_from_fage

    tage, logtau, fburst, fage_burst, const, ftrunc = theta
    
    params['tage'] = tage
    params['logtau'] = logtau
    params['tau'] = tau_from_logtau(logtau)
    
    # params['add_burst'] = True
    params['fburst'] = fburst
    params['fage_burst'] = fage_burst
    params['tburst'] = tburst_from_fage(tage=tage, fage_burst=fage_burst)
    
    params['const'] = const
    
    # params['add_trunc'] = True
    params['ftrunc'] = ftrunc
    params['sf_trunc'] = sf_trunc_from_ftrunc(tage=tage,ftrunc=ftrunc)

    spec = getSpec(obs, sps, **params)

    Bpred = compute_break_strength(spec)

    return Bpred

##################################################################################################################################################

def getSpec(obs, sps, add_burst=True, add_trunc=True, **params):
    
    from prospect.models import SpecModel
    from prospect.models.templates import TemplateLibrary

    model_params = TemplateLibrary['parametric_sfh']

    if add_burst:
        model_params.update(TemplateLibrary['burst_sfh'])
        model_params['const'] = {'N': 1, 'isfree': False, 'init': 0.0, 'units': 'Solar Masses per year'}
        if add_trunc:
            model_params['sf_start'] = {'N': 1, 'isfree': False, 'init': 0.0, 'units': 'Gyrs'}
            model_params['sf_trunc'] = {'N': 1, 'isfree': False, 'init': 0.0, 'units': 'Gyrs'}

    for key in list(model_params.keys()):
        model_params[key]['init'] = params[key]
        model_params[key]['isfree'] = False

    model = SpecModel(model_params)

    spec, _, _, = model.predict(model.theta, obs=obs, sps=sps)

    return spec

##################################################################################################################################################

def log_prior(theta):
    
    tage, logtau, fburst, fage_burst, const, ftrunc = theta

    # This is where all the priors will be set for each parameter in theta
    tage_prior = 0.1 <= tage < 1.775 # Max age is the age of the universe at this redshift
    logtau_prior = -1 <= logtau <= 2
    fburst_prior = 0.0 <= fburst <= 1.0
    fage_burst_prior = 0.0 <= fage_burst <= 1.0
    const_prior = 0.0 <= const <= 1.0
    ftrunc_prior = 0.0 <= ftrunc < 1.0

    f_mass_prior = (const + fburst) <= 1.0

    if tage_prior and logtau_prior and fburst_prior and fage_burst_prior and const_prior and ftrunc_prior and f_mass_prior:
        return 0.0
    return -np.inf #, [tage_prior, logtau_prior, fburst_prior, fage_burst_prior, const_prior, ftrunc_prior]

##################################################################################################################################################

def log_likelihood(theta, B, Berr, params, obs, sps):

    # B, Berr = 6.9, 1.5
    
    Bpred = emcee_model(theta, params, obs, sps)
    return -0.5 * (B - Bpred)**2 / Berr**2

##################################################################################################################################################

def log_probability(theta,B,Berr,params,obs,sps):
    lp = log_prior(theta)

    if not np.isfinite(lp):
        return -np.inf

    return lp + log_likelihood(theta,B,Berr,params,obs,sps)

##################################################################################################################################################
# Defining all the functions for useful plots

def plot_spec(wspec, spec):
    figure(figsize=(8,4))
    plot(wspec[182:4713], spec[182:4713])
    xlabel(r'Rest-frame wavelength [$\AA$]')
    ylabel('Flux [maggies]')
    return

def plot_sfh(params):
    from prospect.plotting.sfh import parametric_sfr
    fig, ax = subplots()
    times = np.linspace(0,params['tage'],10000)
    sfrs = parametric_sfr(times=times,**params)
    ax.plot(times,sfrs,label=f'tage={params['tage']:.2f} Gyrs')

    ax.invert_xaxis()
    ax.set_xlabel('lookback time')
    ax.set_ylabel('sfr')
    return

def plot_imf(bins):

    return

##################################################################################################################################################



##################################################################################################################################################



##################################################################################################################################################
