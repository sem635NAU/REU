import prospect
import os
import fsps
import dynesty
import sedpy
import h5py, astropy
import numpy as np
import pandas as pd
import astroquery
import astropy.units as u
from tqdm import tqdm
from prospect.models.templates import TemplateLibrary, describe
from astropy.cosmology import Planck13
from prospect.models import SpecModel
from prospect.sources import CSPSpecBasis, FastStepBasis
from sedpy import observate 
import seaborn as sns
import matplotlib.pyplot as plt

def get_lsfr_ratio(lmass, z, lsfr):
    """
    We're modeling using a two-value star formation history, both of which are constant. 
    Prospector takes as input the log of the ratio of those two star formation rates. 
    We have the current star formation rate; we next have to find the 'old' star formation rate
        such that the total mass is what we want it to be. 
    This function does that algebra for us. 
    """
    # determine lsfr1 / lsfr2 = lsfr1 / lsfr

    # calculate sfr1
    sfr = 10 ** lsfr
    mass = 10 ** lmass

    t_sfr1 = Planck13.lookback_time(1000).value - Planck13.lookback_time(z).value - t_current_sfr
    sfr1 = (mass - (sfr * t_current_sfr * 1e9)) / (t_sfr1 * 1e9)
    assert sfr1 >= 0, 'sfr1 < 0'

    return np.log10(sfr / sfr1)  

def get_lr(line1, e_line1, line2, e_line2):
    """
    Compute line ratios and uncertainties to avoid doing error propagation every time.
    """
    llr = np.log10(line1 / line2)
    e_lr = abs(llr) / np.log(10) * np.sqrt((e_line1 / line1) ** 2 + (e_line2 / line2) ** 2)
    return llr, e_lr

# fit logU - sSFR relation from https://ui.adsabs.harvard.edu/abs/2018MN/lread_csv('/usr/data/oso2/zlewis/rubies/sel_func/logu_ssfr.csv', names=['ssfr', 'logu'])
m_logu, b_logu = np.polyfit(logu_csv.ssfr, logu_csv.logu, deg=1)

def get_logu(sSFR):
    """
    get logU (excitation parameter) from sSFR (specific star formation rate)
    https://ui.adsabs.harvard.edu/abs/2018MNRAS.477.5568K/abstract
    """
    return b_logu + m_logu * sSFR


# ==================MAKE GRID OF PROSPECTOR GALAXIES======================

t_current_sfr = 0.05 # Gyr; amount of time spent at 'current' star formation rate

n_elem = 30 # number of elements for M*, SFR, z, 12+log(O/H). this is the "resolution"
n_param = 4 # number of parameters we'll vary
lmass_arr = np.linspace(7, 11.5, n_elem)
lsfr_arr = np.linspace(-0.5, 3, n_elem)
z_arr = np.linspace(3, 7, n_elem)
metal_arr = np.linspace(6.5, 9, n_elem)

pros_mat = np.array(np.meshgrid(lmass_arr, lsfr_arr, z_arr, metal_arr)).T
pros_arr = pros_mat.reshape(int(np.prod(pros_mat.shape) / n_param), n_param)

print('Prospector matrix shape:', pros_arr.shape)

# generate empty arrays to fill
f150w_pros = np.zeros(pros_arr.shape[0])
f444w_pros = np.zeros(pros_arr.shape[0])
mass_pros, lsfr_pros, z_pros, metal_pros = pros_arr[:, 0], pros_arr[:, 1], pros_arr[:, 2], pros_arr[:, 3]
has = np.zeros(pros_arr.shape[0])
hbs = np.zeros(pros_arr.shape[0])
oiiis = np.zeros(pros_arr.shape[0])
siis0 = np.zeros(pros_arr.shape[0])
siis1 = np.zeros(pros_arr.shape[0])
weights_sfms = np.zeros(pros_arr.shape[0])
weights_mf = np.zeros(pros_arr.shape[0])
mass_fracs = np.zeros(pros_arr.shape[0]) 

# DEFINE BASICS
sps = FastStepBasis(zcontinuous=1)
filters = observate.load_filters(['jwst_f150w', 'jwst_f444w']) # we'll want F150W and F444W
# see link below to add additional photometry
# https://github.com/bd-j/sedpy/blob/main/sedpy/data/filters/README.md
obs = dict(wave_effective=[1494.4, 4378.7], spectrum=None, unc=None,
           maggies=np.ones(len(filters))*1e-10, maggies_unc=np.ones(len(filters))*1e-10, filters=filters)

# START LOOP
print('Generating prospector models')
for pp in tqdm(range(pros_arr.shape[0])):
    if 10 ** mass_pros[pp] >= (10 ** lsfr_pros[pp]) * t_current_sfr * 1e9: 
    # check that mass > mass formed during current star formation, otherwise nonsensical!
        model_params = TemplateLibrary["continuity_sfh"]
        model_params.update(TemplateLibrary["nebular"])

        # initialize parameters
        model_params['mass']['init'] = 10 ** mass_pros[pp]
        model_params['zred']['init'] = z_pros[pp]
        model_params['gas_logz']['init'] = metal_pros[pp] - 8.69 # units of log(Z / Z_sol)
        model_params['logzsol']['init'] = metal_pros[pp] - 8.69 # units of log(Z / Z_sol)
        
        # model_params['gas_logu']['init'] = -2

        t_now = Planck13.lookback_time(z_pros[pp]).value
        model_params['agebins']['init'] = [[0,                      np.log10((t_now-0.1)*1e9)], 
                                        [np.log10((t_now-0.1)*1e9), np.log10( t_now     *1e9)]] # units of log(years)
        model_params['logsfr_ratios']['init'] = get_lsfr_ratio(mass_pros[pp], z_pros[pp], lsfr_pros[pp])
        model = SpecModel(model_params)

        try:
            mass_fracs[pp] = sps.ssp.stellar_mass / sps.ssp.formed_mass
            # prospector is weird. we give it the total mass, but only get the stellar mass out after the fact.
        except: mass_fracs[pp] = np.nan

        ssfr = 10 ** (lsfr_pros[pp]) / ((10 ** mass_pros[pp]) * mass_fracs[pp])
        model_params['gas_logu']['init'] = get_logu(np.log10(ssfr))
        
        model = SpecModel(model_params)

        # current_parameters = ",".join([f"{p}={v}" for p, v in zip(model.free_params, model.theta)])
        spec, phot, mfrac = model.predict(model.theta, obs=obs, sps=sps)

        emlines = (model._eline_lum * u.L_sun / u.M_sun * (10 ** mass_pros[pp] * u.M_sun))
        emlines = emlines / (4 * np.pi * Planck13.luminosity_distance(z_pros[pp]) ** 2)
        emlines = emlines.to(u.erg / u.s / u.cm ** 2).value
        # https://github.com/cconroy20/fsps/blob/master/data/emlines_info.dat
        has[pp], hbs[pp], oiiis[pp], siis0[pp], siis1[pp] = emlines[74], emlines[59], emlines[62], emlines[77], emlines[78]

        f150w_pros[pp], f444w_pros[pp] = phot[0], phot[1] 

    else: 
        has[pp], hbs[pp], oiiis[pp], siis0[pp], siis1[pp] = np.nan, np.nan, np.nan, np.nan, np.nan
        f150w_pros[pp], f444w_pros[pp] = np.nan, np.nan 

# save into dataframe
pros_df = pd.DataFrame()
pros_df['mass'] = mass_pros
pros_df['lsfr'] = lsfr_pros
pros_df['z'] = z_pros
pros_df['metal'] = metal_pros
pros_df['ha'] = has
pros_df['hb'] = hbs
pros_df['oiii'] = oiiis
pros_df['sii0'] = siis0
pros_df['sii1'] = siis1
pros_df['f444w'] = -2.5 * np.log10(f444w_pros)
pros_df['f150w'] = -2.5 * np.log10(f150w_pros)
pros_df['mass_frac'] = mass_fracs
pros_df['mstar'] = np.log10(10 ** mass_pros * mass_fracs)
pros_df[(pros_df.mstar.values < 0)] = np.nan
pros_df.to_csv(f'/usr/data/oso2/zlewis/rubies/sel_func/pros_df_{n_elem}elem.csv') # change!
print('dataframe saved')