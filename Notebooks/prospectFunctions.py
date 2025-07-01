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

# def getModel(mass=None, zred=None, logzsol=None, tage=None, dust2=None, ldist=10.0, **extras):
#     """Build a prospect.models.SedModel object

#     :param mass: (optional, default:None)
#         If given, produce spectra for this mass. Otherwise the mass will
#         be 1e8 solar masses.

#     :param zred: (optional, default: None)
#         If given, produce spectra and observed frame photometry appropriate
#         for this redshift. Otherwise the redshift will be zero.

#     :param logzsol: (optional, default: None)
#         If given, fix the model metallicity (:math: `log(Z/Z_sun)`) to the given value.
#         Otherwise the metallicity will be set to -0.5.
        
#     :param tage: (optional, default: None)
#         If given, produce spectra and model photometry appropriate for
#         this galactic age. Otherwise the age will be set to 13. Gyrs.

#     :param dust2: (optional, default: None)
#         If given, produce spectra that are appropriate for provided dust
#         attenuation. Otherwise attenuation will be set to 0.6.

#     :param ldist: (optional, default: 10)
#         The luminosity distance (in Mpc) for the model. Spectra and observed
#         frame (apparent) photometry will be appropriate for this luminosity distance.

#     :returns model:
#         An instance of prospect.models.SedModel
#     """
#     from prospect.models.sedmodel import SedModel
#     from prospect.models.templates import TemplateLibrary

#     model_params = TemplateLibrary['parametric_sfh']

#     model_params['lumdist'] = {'N':1, 'isfree':False, 'init':ldist, 'units':'Mpc'}

#     # Change `isfree` so that all parameters that will be kept track of are identified 
#     # in the `model` object as `free_params`
#     model_params['zred']['isfree'] = True
#     model_params['tau']['isfree'] = False

#     if zred is None:
#         model_params['zred']['init'] = 0.0
#     else:
#         model_params['zred']['init'] = zred

#     if mass is not None:
#         model_params['mass']['init'] = mass

#     if logzsol is not None:
#         model_params['logzsol']['init'] = logzsol

#     if tage is None:
#         model_params['tage']['init'] = 13.
#     else:
#         model_params['tage']['init'] = tage

#     if dust2 is not None:
#         model_params['dust2']['init'] = dust2

#     model = SedModel(model_params)

#     return model

def getModel(mass=None, zred=None, logzsol=None, tage=None, dust2=None, **extras):
    """Build a prospect.models.SpecModel object

    :param mass: (optional, default:None)
        If given, produce spectra for this mass. Otherwise the mass will
        be 1e10 solar masses.

    :param zred: (optional, default: None)
        If given, produce spectra and observed frame photometry appropriate
        for this redshift. Otherwise the redshift will be zero.

    :param logzsol: (optional, default: None)
        If given, fix the model metallicity (:math: `log(Z/Z_sun)`) to the given value.
        Otherwise the metallicity will be set to -0.5.
        
    :param tage: (optional, default: None)
        If given, produce spectra and model photometry appropriate for
        this galactic age. Otherwise the age will be set to 13. Gyrs.

    :param dust2: (optional, default: None)
        If given, produce spectra that are appropriate for provided dust
        attenuation. Otherwise attenuation will be set to 0.6.

    :returns model:
        An instance of prospect.models.SedModel
    """
    from prospect.models import SpecModel
    from prospect.models.templates import TemplateLibrary

    model_params = TemplateLibrary['ssp']

    # Change `isfree` so that all parameters that will be kept track of are identified 
    # in the `model` object as `free_params`
    model_params['zred']['isfree'] = True

    if zred is None:
        model_params['zred']['init'] = 0.0
    else:
        model_params['zred']['init'] = zred

    if mass is not None:
        model_params['mass']['init'] = mass

    if logzsol is not None:
        model_params['logzsol']['init'] = logzsol

    if tage is None:
        model_params['tage']['init'] = 13.
    else:
        model_params['tage']['init'] = tage

    if dust2 is not None:
        model_params['dust2']['init'] = dust2

    model = SpecModel(model_params)

    return model

##################################################################################################################################################

def getSps(zcontinuous=1, **extras):
    """
    :param zcontinuous:
        A value of 1 ensures that we use interpolation between SSPs to
        have a continuous metallicity parameter (`logzsol`)
        See python-FSPS documentation for details
    """
    from prospect.sources import CSPSpecBasis
    sps = CSPSpecBasis(zcontinuous=zcontinuous)
    return sps

##################################################################################################################################################

def getTheta(model_obj, obs_obj, sps_obj):
    theta = model_obj.theta.copy()
    init_spec, init_phot, init_mfrac = model_obj.sed(theta, obs=obs_obj, sps=sps_obj)
    return init_spec, init_phot, init_mfrac

##################################################################################################################################################

def getWave(obs=None, sps=None, zred=0.0, **extras):

    wphot = obs['phot_wave']

    a = 1.0 + zred

    wspec = sps.wavelengths
    wspec *= a

    return wspec, wphot

##################################################################################################################################################

def plotBalmerBreak(zred=None, **extras):

    a = 1.0 + zred
    
    balmer = np.ones(2)*3646*a
    
    y = np.linspace(0,1e10,2)
    
    lower_1 = np.ones(2)*3620*a
    upper_1 = np.ones(2)*3720*a
    
    lower_2 = np.ones(2)*4000*a
    upper_2 = np.ones(2)*4100*a

    alpha = 0.6
    lw = .7
    ls = '--'
    
    # plot(balmer,y,color='black',alpha=alpha+.2)
    
    plot(lower_1,y,color='purple',ls=ls,lw=lw,alpha=alpha)
    plot(upper_1,y,color='purple',ls=ls,lw=lw,alpha=alpha)
    
    plot(lower_2,y,color='red',ls=ls,lw=lw,alpha=alpha)
    plot(upper_2,y,color='red',ls=ls,lw=lw,alpha=alpha)

##################################################################################################################################################

def plotBalmerBreakD4000(zred=None, **extras):

    a = 1.0 + zred
    
    balmer = np.ones(2)*3646*a
    
    y = np.linspace(0,1e10,2)
    
    lower_1 = np.ones(2)*3850*a
    upper_1 = np.ones(2)*3950*a
    
    lower_2 = np.ones(2)*4000*a
    upper_2 = np.ones(2)*4100*a

    alpha = 0.6
    lw = .7
    ls = '--'
    
    # plot(balmer,y,color='black',alpha=alpha+.2)
    
    plot(lower_1,y,color='blue',ls=ls,lw=lw,alpha=alpha)
    plot(upper_1,y,color='blue',ls=ls,lw=lw,alpha=alpha)
    
    plot(lower_2,y,color='red',ls=ls,lw=lw,alpha=alpha)
    plot(upper_2,y,color='red',ls=ls,lw=lw,alpha=alpha)

##################################################################################################################################################

def plotBalmerBreakRestFrame(zred=None, **extras):

    # a = 1.0 + zred
    
    balmer = np.ones(2)*3646
    
    y = np.linspace(0,1e10,2)
    
    lower_1 = np.ones(2)*3620
    upper_1 = np.ones(2)*3720

    lower_2 = np.ones(2)*3850
    upper_2 = np.ones(2)*3950
    
    lower_3 = np.ones(2)*4000
    upper_3 = np.ones(2)*4100

    alpha = 0.6
    lw = 1
    ls = '--'
    
    plot(lower_1,y,color='purple',ls=ls,lw=lw,alpha=alpha)
    plot(upper_1,y,color='purple',ls=ls,lw=lw,alpha=alpha)
    
    plot(lower_2,y,color='blue',ls=ls,lw=lw,alpha=alpha)
    plot(upper_2,y,color='blue',ls=ls,lw=lw,alpha=alpha)

    plot(lower_3,y,color='red',ls=ls,lw=lw,alpha=alpha)
    plot(upper_3,y,color='red',ls=ls,lw=lw,alpha=alpha)

##################################################################################################################################################

    
# def getBalmerStrength(spec):
#     """Compute the balmer strength of a spectra as calculated in De Graaff et. al. 2025.
#     :param spec:
#         The spectra of an object outputted after using the model.sed() function on 
#         a prospect.models.sedmodel.SedModel() (model) object. The average is calculated
#         between the ranges of 3620-3720 and 4000-4100 Angstroms.
#     :returns balmer_strength:
#     """
#     return spec[111:113].mean() - spec[119:121].mean()

##################################################################################################################################################
def getBreakBoundsAnna(wspec, zred=None, **extras):

    a = 1.0 + zred

    for i,s in enumerate(wspec>3620*a):
        if s:
            # print(s, i, wspec[i], wspec[i]/a)
            blue_lower = i
            break
    
    for i,s in enumerate(wspec<3720*a):
        if not s:
            # print(s, i-1, wspec[i-1], wspec[i-1]/a)
            blue_upper = i-1
            break
    
    for i,s in enumerate(wspec>4000*a):
        if s:
            # print(s, i, wspec[i], wspec[i]/a)
            red_lower = i
            break
    
    for i,s in enumerate(wspec<4100*a):
        if not s:
            # print(s, i-1, wspec[i-1], wspec[i-1]/a)
            red_upper = i-1
            break

    return {'blue':[blue_lower, blue_upper], 'red':[red_lower, red_upper]}

##################################################################################################################################################

def getBreakBoundsD4000(wspec, zred=None, **extras):

    a = 1.0 + zred

    for i,s in enumerate(wspec>3850*a):
        if s:
            # print(s, i, wspec[i], wspec[i]/a)
            blue_lower = i
            break
    
    for i,s in enumerate(wspec<3950*a):
        if not s:
            # print(s, i-1, wspec[i-1], wspec[i-1]/a)
            blue_upper = i-1
            break
    
    for i,s in enumerate(wspec>4000*a):
        if s:
            # print(s, i, wspec[i], wspec[i]/a)
            red_lower = i
            break
    
    for i,s in enumerate(wspec<4100*a):
        if not s:
            # print(s, i-1, wspec[i-1], wspec[i-1]/a)
            red_upper = i-1
            break

    return {'blue':[blue_lower, blue_upper], 'red':[red_lower, red_upper]}

##################################################################################################################################################

def getBreakBounds(wspec, start, zred=None, **extras):
    a = 1.0 + zred

    for i,s in enumerate(wspec>start*a):
        if s:
            # print(i, wspec[i]/a)
            blue_lower = i
            break
    
    for i,s in enumerate(wspec<(start+100)*a):
        if not s:
            # print(i-1, wspec[i-1]/a)
            blue_upper = i
            break
    
    for i,s in enumerate(wspec>4000*a):
        if s:
            # print(i, wspec[i]/a)
            red_lower = i
            break
    
    for i,s in enumerate(wspec<4100*a):
        if not s:
            # print(i-1, wspec[i-1]/a)
            red_upper = i
            break
            
    return {'blue':[blue_lower, blue_upper], 'red':[red_lower, red_upper]}

##################################################################################################################################################
    
def getParams():
    grid_ranges = {}
    grid_ranges['logzsol'] = np.linspace(-1,.5,10)
    grid_ranges['dust2'] = np.zeros(1)
    grid_ranges['tage'] = np.logspace(-2,1,10)
    
    run_params = {}
    run_params['zred'] = 3.548
    run_params['mass'] = 1e8
    run_params['add_duste'] = False
    run_params['zcontinuous'] = 1
    return grid_ranges, run_params
    
##################################################################################################################################################

