import numpy as np
import sys

sys.path.append('/home/samantha/Documents/REU/scripts/') # Use the absolute path to the directory
from prospectFunctionsSFH import *

print('Initializing parameters...\n')
run_params, obs, sps, wspec = init_prospect_generation()

B = 6.9 # Lower
Berr = 1.2

tage = 0.3
tau = 0.1
# fburst = 
# fage_burst = 
theta = [tage, tau]#, fburst, fage_burst, ftrunc, const] #fburst, fage instead
# Which accounts for tage, tau, fburst, tburst, sf_start, sf_trunc, const

print('Calculating probability...\n')
# test_prob = log_probability(theta,B,Berr,run_params,obs,sps)
data = (B,Berr,run_params,obs,sps)


print(test_prob)
print('Done')

# To run, open terminal and type, `python SFHModelsGeneration.py`

