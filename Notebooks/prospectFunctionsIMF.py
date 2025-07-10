


def getModel(mass=None, zred=None, logzsol=None, tage=None, dust2=None, imf_type=2, sfh_type=4, **extras):
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
        

    model_params['imf_type']['init'] = imf_type

    model = SpecModel(model_params)

    return model


##################################################################################################################################################









##################################################################################################################################################








