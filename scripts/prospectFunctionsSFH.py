# IMPORT STATEMENTS:
import numpy as np

##################################################################################################################################################

# def getModel(mass=None, zred=None, logzsol=None, tage=None, dust2=None, imf_type=None, 
#              sfh_type=None, sf_start=None, sf_trunc=None, fburst=None, fage_burst=None, **extras):
#     """Build a prospect.models.SpecModel object

#     :param mass: (optional, default:None)
#         If given, produce spectra for this mass. Otherwise the mass will
#         be 1e10 solar masses.

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

#     :returns model:
#         An instance of prospect.models.SedModel
#     """
#     from prospect.models import SpecModel
#     from prospect.models.templates import TemplateLibrary

#     model_params = TemplateLibrary['parametric_sfh']
#     model_params.update(TemplateLibrary['burst_sfh'])

#     # Change `isfree` so that all parameters that will be kept track of are identified 
#     # in the `model` object as `free_params`
#     # model_params['zred']['isfree'] = True

#     if sf_start is not None:
#         model_params['sf_start'] = {'N': 1, 'isfree': False, 'init': sf_start, 'units': 'Gyrs'}
#     if sf_trunc is not None:
#         model_params['sf_trunc'] = {'N': 1, 'isfree': False, 'init': sf_trunc, 'units': 'Gyrs'}

#     if zred is not None:
#         model_params['zred']['init'] = zred

#     if mass is not None:
#         model_params['mass']['init'] = mass

#     if logzsol is not None:
#         model_params['logzsol']['init'] = logzsol

#     if tage is not None:
#         model_params['tage']['init'] = tage

#     if dust2 is not None:
#         model_params['dust2']['init'] = dust2

#     if imf_type is not None:
#         model_params['imf_type']['init'] = imf_type

#     if sfh_type is not None:
#         model_params['sfh']['init'] = sfh_type

#     if add_burst:
#         model_params.update(TemplateLibrary['burst_sfh'])

#     else:
#         for key in extras.keys()
    

#     model = SpecModel(model_params)

#     return model

##################################################################################################################################################

def model_function(add_burst=False, add_trunc=False, **extras):

    from prospect.models import SpecModel
    from prospect.models.templates import TemplateLibrary

    model_params = TemplateLibrary['parametric_sfh']

    if add_burst:
        model_params.update(TemplateLibrary['burst_sfh'])
        model_params['const'] = {'N': 1, 'isfree': False, 'init': 0.0, 'units': 'Solar Masses per year'}
        if add_trunc:
            model_params['sf_start'] = {'N': 1, 'isfree': False, 'init': 0.0, 'units': 'Gyrs'}
            model_params['sf_trunc'] = {'N': 1, 'isfree': False, 'init': 0.0, 'units': 'Gyrs'}

    for key in list(extras.keys()):
        model_params[key]['init'] = extras[key]

    model = SpecModel(model_params)
    
    return model

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
    return sps
    
##################################################################################################################################################









##################################################################################################################################################









##################################################################################################################################################









##################################################################################################################################################









##################################################################################################################################################









##################################################################################################################################################
