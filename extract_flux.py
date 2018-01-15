import numpy as np
import matplotlib.mlab as mlab
import datetime as dt


def extract_met(fname, sname):
    re_arr = mlab.csv2rec(fname)
    time_step_start = re_arr['timestamp_start']
    sw_down = re_arr['sw_in_f']
    lw_down = re_arr['lw_in_f']
    precip = (1/(30.*60.))*re_arr['p_f']
    print precip[precip>0.]
    press = 1000.*re_arr['pa_f']
    air_temp = re_arr['ta_f'] + 273.15
    wind = re_arr['ws_f']
    vpd = re_arr['vpd_f']
    q = vpd2q_air(vpd, air_temp, press)
    met_dat = np.column_stack((sw_down, lw_down, precip, air_temp, wind, press, q))
    np.savetxt(sname, met_dat, fmt=['%1.1f', '%1.1f', '%1.3e', '%1.2f', '%1.3f', '%1.1f', '%1.3e'],
               header='swdown lwdown precip airtemp wind pressure specific_humidity')
    np.savetxt('timesteps.dat', time_step_start, fmt='%1.0f', header='start of timestep')
    return 'done'


def rh2q_air(rh, t_air, pressure):
    tc = t_air - 273.15
    es = 6.112 * np.exp((17.67 *tc) / (tc + 243.5))
    e = (rh/100)*es
    p_mb = pressure/100.0
    q_air = (0.622*e) / (p_mb - (0.378*e))
    return q_air


def vpd2q_air(vpd, t_air, pressure):
    tc = t_air - 273.15
    es = 6.112 * np.exp((17.67 *tc) / (tc + 243.5))
    rh = 100*(vpd + es) / es
    e = (rh/100)*es
    p_mb = pressure/100.0
    q_air = (0.622*e) / (p_mb - (0.378*e))
    return q_air


def extract_assim_day(fname, sname, year):
    re_arr = mlab.csv2rec(fname)
    time_step = np.array([dt.datetime.strptime(str(timestamp), '%Y%m%d') for timestamp in re_arr['timestamp']])
    years = np.array([date.year for date in time_step])
    yr_idx = np.where(years==year)[0]
    times = re_arr['timestamp'][yr_idx]
    le = re_arr['le_f_mds'][yr_idx]
    le_qc = re_arr['le_f_mds_qc'][yr_idx]
    le_unc = re_arr['le_randunc'][yr_idx]
    sm1 = re_arr['swc_f_mds_1'][yr_idx]
    sm1_qc = re_arr['swc_f_mds_1_qc'][yr_idx]
    # sm2 = re_arr['swc_f_mds_2'][yr_idx]
    # sm2_qc = re_arr['swc_f_mds_2_qc'][yr_idx]
    nee = re_arr['nee_cut_ref'][yr_idx]
    nee_qc = re_arr['nee_cut_ref_qc'][yr_idx]
    gpp = re_arr['gpp_nt_cut_ref'][yr_idx]
    gpp_se = re_arr['gpp_nt_cut_se'][yr_idx]
    assim_dat = np.column_stack((times, le, le_qc, le_unc, sm1, sm1_qc, nee, nee_qc, gpp, gpp_se))
    np.savetxt(sname, assim_dat,
               header='times, le, le_qc, le_unc, sm1, sm1_qc, nee, nee_qc, gpp, gpp_se')
    return 'done'