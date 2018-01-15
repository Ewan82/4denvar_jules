import numpy as np
import netCDF4 as nc


def extract_twin_data(mod_truth='output/model_truth/mod_truth.daily.nc'):
    mod_dat = nc.Dataset(mod_truth, 'r')
    npp_ob_position = np.array([179,  79, 208, 106,  33,  42, 319, 175, 152,  21, 326, 251, 124,
       131, 327, 263, 192, 177,  99, 182, 287, 268, 148, 146, 193, 209,
       260, 154, 348, 201, 119,  36, 355,  14,  32, 104, 336, 226,  94,
       354, 165,  13, 318, 159, 361, 254, 340,   5, 364, 189])
    sm_obs_position = np.array([ 17, 296, 274, 194, 299, 360, 267, 252, 127, 214, 305, 237, 334,
        89,  78, 254, 249,  62, 228, 126,   6, 310, 134,  26, 117, 338,
       122, 287,  75,  63, 355, 157, 242, 290, 283, 322, 230,  64,  65,
        21, 171, 358, 149,  36, 251,  92,  13,  88,  18,  31])
    npp_obs = mod_dat.variables['npp'][npp_ob_position, 5, 0, 0]
    sm_obs = mod_dat.variables['smcl'][sm_obs_position, 0, 0, 0]
    return npp_obs, npp_ob_position, sm_obs, sm_obs_position


def extract_day_data(obs='metolius_day_assim.csv'):
    assim_dat = np.loadtxt(obs)
    gpp = assim_dat[:, 8]
    gpp_qc = assim_dat[:, 7]
    gpp_err = assim_dat[:, 9]
    gpp_ob_position = np.where(gpp_qc==1.0)[0]
    gpp_obs = gpp[gpp_ob_position]
    gpp_err = gpp_err[gpp_ob_position]
    return gpp_obs, gpp_ob_position, gpp_err