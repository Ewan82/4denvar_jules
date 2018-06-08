import numpy as np
import netCDF4 as nc
import matplotlib.mlab as mlab


def extract_twin_data(mod_truth='output/test/williams_truth.daily.nc'):
    mod_dat = nc.Dataset(mod_truth, 'r')
    gpp_ob_position = np.array([ 119, 123, 125, 126, 127, 128, 129, 130, 131,
       134, 135, 139, 141, 143, 146, 148, 155, 157, 158, 160, 161, 162,
       163, 164, 165, 166, 167, 168, 169, 170, 171, 172, 173, 175, 176,
       177, 179, 180, 182, 183, 184, 186, 187, 190, 191, 194, 196, 197,
       198, 199, 200, 201, 202, 203, 204, 206, 207, 208, 209, 210, 211,
       213, 214, 215, 218, 219, 220, 222, 223, 225, 228, 229, 230, 231,
       233, 234, 239, 240, 241, 242, 243, 244, 245, 247, 248, 249, 251,
       252, 259, 260])

        #np.array([179,  79, 208, 106,  33,  42, 319, 175, 152,  21, 326, 251, 124,
       #131, 327, 263, 192, 177,  99, 182, 287, 268, 148, 146, 193, 209,
       #260, 154, 348, 201, 119,  36, 355,  14,  32, 104, 336, 226,  94,
       #354, 165,  13, 318, 159, 361, 254, 340,   5, 363, 189, 230, 206, 181, 183, 186, 190])
    lai_obs_position = np.array([ 119, 123, 125, 126, 127, 128, 129, 130, 131,
       134, 135, 139, 141, 143, 146, 148, 155, 157, 158, 160, 161, 162,
       163, 164, 165, 166, 167, 168, 169, 170, 171, 172, 173, 175, 176,
       177, 179, 180, 182, 183, 184, 186, 187, 190, 191, 194, 196, 197,
       198, 199, 200, 201, 202, 203, 204, 206, 207, 208, 209, 210, 211,
       213, 214, 215, 218, 219, 220, 222, 223, 225, 228, 229, 230, 231,
       233, 234, 239, 240, 241, 242, 243, 244, 245, 247, 248, 249, 251,
       252, 259, 260])

        #np.array([ 17, 296, 274, 194, 299, 360, 267, 252, 127, 214, 305, 237, 334,
        #89,  78, 254, 249,  62, 228, 126,   6, 310, 134,  26, 117, 338,
       #122, 287,  75,  63, 355, 157, 242, 290, 283, 322, 230,  64,  65,
        #21, 171, 358, 149,  36, 251,  92,  13,  88,  18,  31, 150, 175, 180, 183, 187, 191, 192, 193])
    gpp_obs_tr = 1000 * 60 * 60 * 24 * mod_dat.variables['gpp'][gpp_ob_position, 7, 0, 0]
    lai_obs_tr = mod_dat.variables['lai'][lai_obs_position, 7, 0, 0]
    gpp_err = np.mean(gpp_obs_tr[gpp_obs_tr > 0]) * 0.1  # 10% error in mod obs
    lai_err = np.mean(lai_obs_tr[lai_obs_tr > 0]) * 0.05
    gpp_obs = np.array([gpp + np.random.normal(0.0, gpp_err) for gpp in gpp_obs_tr])
    lai_obs = np.array([lai + np.random.normal(0.0, lai_err) for lai in lai_obs_tr])
    return gpp_obs, gpp_ob_position, lai_obs, lai_obs_position


def extract_day_data(obs='mead_dat_assim_2010.csv'):
    assim_dat = np.loadtxt(obs)
    gpp = assim_dat[:, 8]
    gpp_qc = assim_dat[:, 7]
    gpp_err = assim_dat[:, 9]
    gpp_ob_position = np.where(gpp_qc==1.0)[0]
    gpp_obs = gpp[gpp_ob_position]
    gpp_err = gpp_err[gpp_ob_position]
    return gpp_obs, gpp_ob_position, gpp_err


def extract_gpp_day_data(obs='mead_dat_assim_2010.csv'):
    assim_dat = mlab.csv2rec(obs)[:-1]
    gpp = assim_dat['gpp']
    gpp_qc = assim_dat['nee_qc']
    gpp_err = assim_dat['gpp_se']
    gpp_ob_position = np.where(gpp_qc == 1.0)[0]
    gpp_obs = gpp[gpp_ob_position]
    gpp_err = gpp_err[gpp_ob_position]
    return gpp_obs, gpp_ob_position, gpp_err


def extract_energy_day_data(obs='metolius_day_assim.csv'):
    assim_dat = mlab.csv2rec(obs)[:-1]
    le = assim_dat['le']
    le_qc = assim_dat['le_qc']
    le_err = assim_dat['le_unc']
    le_ob_position = np.where(le_qc == 1.0)[0]
    le_obs = le[le_ob_position]
    le_err = le_err[le_ob_position]

    h = assim_dat['h']
    h_qc = assim_dat['h_qc']
    h_err = assim_dat['h_unc']
    h_ob_position = np.where(h_qc == 1.0)[0]
    h_obs = h[h_ob_position]
    h_err = h_err[h_ob_position]

    g = assim_dat['g']
    g_qc = assim_dat['g_qc']
    g_ob_position = np.where(g_qc == 1)[0]
    g_obs = g[g_ob_position]
    return le_obs, le_ob_position, le_err, h_obs, h_ob_position, h_err, g_obs, g_ob_position


def extract_sm_day_data(obs='metolius_day_assim.csv'):
    assim_dat = mlab.csv2rec(obs)[:-1]
    sm = assim_dat['sm1']/100.
    sm_qc = assim_dat['sm1_qc']
    sm_ob_position = np.where(sm_qc == 1)[0]
    sm_obs = sm[sm_ob_position]
    return sm_obs, sm_ob_position