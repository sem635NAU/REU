from matplotlib.pyplot import *

def getModel(mass=None, zred=None, tage=None, metallicity=None, dust2=None, ldist=10.0, add_duste=False, **extras):
    from prospect.models.sedmodel import SedModel
    from prospect.models.templates import TemplateLibrary
    from prospect.models import priors
    
    model_params = TemplateLibrary['parametric_sfh']
    
    # If necessary, this is the step where you would need to change the initial values for the parameters
    model_params['lumdist'] = {'N':1, 'isfree':False, 'init':ldist, 'units':'Mpc'}
    # model_params['dust2']['init'] = 0.05 # 0,3 or 4 (n=6)
    model_params['logzsol']['init'] = -0.5 #-1,0.5 (n=4)
    # model_params['tage']['init'] = 13.
    
    # Make it so these parameters can be varied
    # Variable masses
    if mass == None:
        model_params['mass']['init'] = 1e8    #################################################################################################
        model_params['mass']['prior'] = priors.LogUniform(mini=1e6, maxi=1e10)  # MAKE SURE THAT THE PRIORS CHANGE NOTHING ABOUT THE SED MODELS
        model_params['mass']['disp_floor'] = 1e6    ###########################################################################################
    else:
        model_params['mass']['init'] = mass
        model_params['mass']['prior'] = priors.LogUniform(mini=10**(np.log10(mass)-2), maxi=10**(np.log10(mass)+2))
        model_params['mass']['disp_floor'] = 10**(np.log10(mass)-2)
    
    # Variable redshifts
    if zred == None:
        model_params['zred']['init'] = 0.0
    else:
        model_params['zred']['isfree'] = False
        model_params['zred']['init'] = zred
        
    # Variable age of the galaxy
    if tage == None:
        model_params['tage']['init'] = 13.
    else:
        model_params['tage']['init'] = tage

    # Variable metallicities
    if metallicity == None:
        model_params['logzsol']['init'] = -0.5
    else:
        model_params['logzsol']['isfree'] = False
        model_params['logzsol']['init'] = metallicity

    # Variable dust parameters
    if dust2 == None:
        model_params['dust2']['init'] = 0.05
    else:
        model_params['dust2']['init'] = dust2
        
    # Changing the priors for some of the free parameters
    model_params["dust2"]["prior"] = priors.TopHat(mini=0.0, maxi=2.0)
    model_params['tau']['prior'] = priors.LogUniform(mini=1e-1, maxi=1e2)
    
    # Providing minimum scale for the cloud of walkers for emcee
    model_params['tau']['disp_floor'] = 1.0
    model_params['tage']['disp_floor'] = 1.0
    
    # If optional (**extras) parameters are provided...
    # if fixed_metallicity is not None:
    #     model_params['logzsol']['isfree'] = False
    #     # make it a fixed parameter
    #     model_params['logzsol']['init'] = fixed_metallicity
        
    # if object_redshift is not None:
    #     model_params['zred']['isfree'] = False
    #     # make it a fixed parameter
    #     model_params['zred']['init'] = object_redshift
    
    if add_duste:
        model_params.update(TemplateLibrary['dust_emission'])
        
    model = SedModel(model_params)
    
    return model

# -------------------------------------------------------------------------------------------------------------------

# Building the SPS function
def getSps(zcontinuous=1, **extras):
    from prospect.sources import CSPSpecBasis
    sps = CSPSpecBasis(zcontinuous=zcontinuous)
    return sps

# -------------------------------------------------------------------------------------------------------------------

# Building the observations function
def getObs(snr=10, ldist=10.0, **extras_):
    
    from prospect.utils.obsutils import fix_obs
    import sedpy
    import numpy as np
    
    # obs is a dictionary of observational data to use for the fit
    obs = {}
    
    galex = ['galex_'+a for a in ['FUV', 'NUV']]
    sdss = [f'sdss_{c}0' for c in ['u', 'g', 'r', 'i', 'z']]
    spitzer = ['spitzer_irac_ch'+b for b in ['1', '2', '3', '4']]
    filternames = galex + sdss + spitzer
    
    obs['filters'] = sedpy.observate.load_filters(filternames)
    
    M_AB = np.array(obj)

    dm = 25 + 5.0 * np.log10(ldist)
    mags = M_AB + dm
    
    obs['maggies'] = 10**(-0.4*mags)
    obs['maggies_unc'] = (1./snr) * obs['maggies']
    
    obs['phot_mask'] = np.array(['spitzer' not in f.name for f in obs['filters']])
    obs['phot_wave'] = np.array([f.wave_effective for f in obs['filters']])

    obs['wavelength'] = None
    obs['spectrum'] = None
    obs['unc'] = None
    obs['mask'] = None
    
    obs = fix_obs(obs)
    
    return obs

# -------------------------------------------------------------------------------------------------------------------

# Generate the theta for the model and the SED
def getTheta(model_obj, obs_obj, sps_obj):
    theta = model_obj.theta.copy()
    init_spec, init_phot, init_mfrac = model_obj.sed(theta, obs=obs_obj, sps=sps_obj)
    return init_spec, init_phot, init_mfrac

# -------------------------------------------------------------------------------------------------------------------

# To plot the model from the model, obs, sps, and SED
def getModelPlot(model,obs,sps,init_spec,init_phot):
    
    figure(figsize=(16,8))
    
    wphot = obs['phot_wave']
    
    a = 1.0 + model.params.get('zred',0.0)
    
    if obs['wavelength'] is None:
        wspec = sps.wavelengths
        wspec *= a
    else:
        wspec = obs['wavelength']
    
    loglog(wspec, init_spec)
    errorbar(wphot, init_phot, marker='s', color='gray', alpha=0.8,
             ls='', lw=3, markeredgewidth=3)
    
    xmin, xmax = np.min(wphot)*0.8, np.max(wphot)/0.8
    temp = np.interp(np.linspace(xmin,xmax,10000), wspec, init_spec)
    ymin, ymax = temp.min()*0.8, temp.max()/0.4

    title(f"mass={model.params.get('mass',0.0)[0]:.1E}, zred={model.params.get('zred',0.0)[0]}, tage={model.params.get('tage',0.0)[0]}", fontsize=15)
    xlabel('Wavelength [A]', fontsize=15)
    ylabel('Flux Density [maggies]', fontsize=15)
    xlim([xmin, xmax])
    ylim([ymin, ymax])
    xscale('log')
    yscale('log')
#     legend(fontsize=15)
    tight_layout()
    show()
    