import numpy as np
import netCDF4 as nc
import datetime as dt
import scipy.optimize as spop
import multiprocessing as mp
import itertools as itt
import py_jules as pyj
import os
import shutil as sh
import twin_obs
import glob
import subprocess
import sys
# Rewrite this as a class?


class Jules_DA:
    def __init__(self,):
        self.jules_class = pyj.jules()  # python jules wrapper class
        # Set output dirctory
        self.output_dir = "'./output/ens_run_15',"
        self.steps = []
        # mod_truth and background
        self.x_truth = np.array([5.70e-4, 6.631])
        self.xb = np.array([0.033, -10., 26., 4.])
        self.b_mat = np.eye(4)*(self.xb*0.25)**2  # 20% background error
        self.xbs = np.array([[  3.15254089e-02,  -1.51112488e+01,   3.32115984e+01,
                              2.04146839e+00],
                           [  3.75116380e-02,  -8.51150432e+00,   2.67364269e+01,
                              2.69807469e+00],
                           [  2.85951727e-02,  -7.59460748e+00,   2.98701871e+01,
                              3.43684267e+00],
                           [  2.28919224e-02,  -8.68533174e+00,   2.19438489e+01,
                              5.92473048e+00],
                           [  2.61133472e-02,  -8.44272991e+00,   2.87963898e+01,
                              4.59764351e+00],
                           [  5.48199432e-02,  -5.97917631e+00,   2.92305578e+01,
                              3.47107492e+00],
                           [  2.85453696e-02,  -1.17593473e+01,   2.80719181e+01,
                              3.61957078e+00],
                           [  1.53686223e-02,  -9.70126236e+00,   2.48889156e+01,
                              3.45158127e+00],
                           [  1.47414162e-02,  -7.14499330e+00,   2.43584233e+01,
                              3.35083403e+00],
                           [  3.23544229e-02,  -8.70239212e+00,   2.15481979e+01,
                              3.86739605e+00],
                           [  3.91390913e-02,  -1.09058027e+01,   2.93737427e+01,
                              2.68096827e+00],
                           [  3.89344194e-02,  -1.13332968e+01,   3.04539943e+01,
                              1.93499291e+00],
                           [  2.78233999e-02,  -9.86222376e+00,   2.37311911e+01,
                              3.81188239e+00],
                           [ -1.53871212e-03,  -7.57942786e+00,   2.53587519e+01,
                              5.32989346e+00],
                           [  3.43994348e-02,  -1.19694185e+01,   2.12123137e+01,
                              4.86842310e+00]])

  # array of xbi's drawn from background cov mat
        self.size_ens = 15
        # abs(np.random.multivariate_normal(self.xb, self.b_mat))f

    def make_twin_obs(self):
        # extract mod obs
        self.npp_obs, self.npp_ob_position, self.sm_obs, self.sm_ob_position = twin_obs.extract_data()
        self.npp_err = np.ones(len(self.npp_obs))*np.mean(self.npp_obs[self.npp_obs > 0])*0.1  # 10% error in mod obs
        self.sm_err = np.ones(len(self.sm_obs))*np.mean(self.sm_obs)*0.1  # 10% error in mod obs
        self.yoblist = np.append(self.npp_obs, self.sm_obs)
        self.yerr = np.append(self.npp_err, self.sm_err)
        self.rmatrix = np.eye(len(self.yerr))*self.yerr**2

    def make_obs(self):
        self.gpp_obs, self.gpp_ob_position, self.gpp_err = twin_obs.extract_day_data()
        self.yoblist = self.gpp_obs
        self.yerr = 0.1*self.yoblist
        self.rmatrix = np.eye(len(self.yerr))*self.yerr**2

    def make_twin_hxb(self):
        self.npp_xb, self.npp_ob_position, self.sm_xb, self.sm_ob_position = twin_obs.extract_data(
            'output/background/xb.daily.nc')
        self.hxb = np.append(self.npp_xb, self.sm_xb)
        self.xb_mat = (1./(np.sqrt(self.size_ens-1)))*np.array([xbi - self.xb for xbi in self.xbs])
        self.xb_mat_inv = np.linalg.pinv(self.xb_mat.T)

    def make_hxb(self):
        xb_nc = nc.Dataset('output/test/metol.daily.nc', 'r')
        self.hxb = 1000 * 60 * 60 * 24 * xb_nc.variables['gpp'][self.gpp_ob_position, 1, 0, 0]
        self.xb_mat = (1./(np.sqrt(self.size_ens-1)))*np.array([xbi - self.xb for xbi in self.xbs])
        self.xb_mat_inv = np.linalg.pinv(self.xb_mat.T)

    def run_jules(self, neff=8.00e-4, b=6.631, n_leaf=0.033, tlow=-10, tupp=26, lai=4, run_id='ens'):
        """
        Runs JULES changing soil parameters
        :param run_id: id of run as a string
        :return: location of JULES output as string
        """
        self.jules_class.output_nml.mapping["jules_output_1_run_id"] = "'" + run_id + "',"
        self.jules_class.output_nml.mapping["jules_output_1_output_dir"] = self.output_dir
        self.jules_class.jules_vegetation_nml.mapping["jules_vegetation_1_l_phenol"] = ".false.,"
        self.jules_class.timesteps_nml.mapping["jules_time_1_main_run_start"] = " '2010-01-01 00:00:00',"
        self.jules_class.timesteps_nml.mapping["jules_spinup_1_max_spinup_cycles"] = " 2"
        self.jules_class.ancillaries_nml.mapping["jules_soil_props_1_const_val"] = str(b) + ", 0.3967309, 0.0027729999, 0.45809999, " \
                                                                    "0.3283205, 0.1866, 1185786.0, 0.2269195, 0.17,"
        self.jules_class.pft_params_nml.mapping["jules_pftparm_1_neff_io"] = "8.00e-4,8.00e-4,8.00e-4,4.00e-4,8.00e-4," + str(neff) + \
                                                              ", 8.00e-4,4.00e-4,8.00e-4,"
        self.jules_class.pft_params_nml.mapping["jules_pftparm_1_tlow_io"] = " 0," + str(tlow) + ",0,13,0,0,0,13,0,"
        self.jules_class.pft_params_nml.mapping["jules_pftparm_1_tupp_io"] = " 36," + str(tupp) + ",36,45,36,46,36,45,36,"
        self.jules_class.pft_params_nml.mapping["jules_pftparm_1_nl0_io"] = " 0.046," + str(n_leaf) + ",0.073,0.06,0.06,0.073,0.073,0.06,0.073,"
        self.jules_class.pft_params_nml.mapping["jules_pftparm_1_lai_io"] = " 5," + str(lai) + ",2,4,1,2,2,4,2,"
        self.jules_class.runJules_print()
        return self.output_dir + "/" + run_id + '.outvars.nc'

    def create_ensemble(self):
        for xbi in enumerate(self.xbs):
            self.run_jules(n_leaf=xbi[1][0], tlow=xbi[1][1], tupp=xbi[1][2], lai=xbi[1][3], run_id='ens'+str(xbi[0]))
        dumps = glob.glob('output/ens_run_15/ens*.dump*')
        for f in dumps:
            os.remove(f)
        hm_xbs = []
        for xb_fname in glob.glob('output/ens_run_15/ens*.daily.nc'):
            xbi_nc = nc.Dataset(xb_fname, 'r')
            hxbi = 1000 * 60 * 60 * 24 * xbi_nc.variables['gpp'][self.gpp_ob_position, 1, 0, 0]
            hm_xbs.append(hxbi)
        self.hmx_mat = (1. / (np.sqrt(self.size_ens - 1))) * np.array([hmxb - self.hxb for hmxb in hm_xbs])

    def xvals2wvals(self, xvals):
        # return np.dot(self.xb_mat_inv, (xvals - np.mean(self.xbs, axis=0)))
        return np.dot(self.xb_mat_inv, (xvals - self.xb))
    def wvals2xvals(self, wvals):
        # return np.mean(self.xbs, axis=0) + np.dot(self.xb_mat.T, wvals)
        return self.xb + np.dot(self.xb_mat.T, wvals)

    def obcost_ens_inc(self, wvals):
        """Observational part of cost fn.
        """
        return np.dot(np.dot((np.dot(self.hmx_mat.T, wvals) + self.hxb - self.yoblist), np.linalg.inv(self.rmatrix)),
                      (np.dot(self.hmx_mat.T, wvals) + self.hxb - self.yoblist).T)

    def cost_ens_inc(self, wvals):
        modcost = np.dot(wvals, wvals.T)
        obcost = self.obcost_ens_inc(wvals)
        ret_val = 0.5 * modcost + 0.5*obcost
        return ret_val

    def gradcost_ens_inc(self, wvals):
        obcost = np.dot(self.hmx_mat, np.dot(np.linalg.inv(self.rmatrix),
                                           (np.dot(self.hmx_mat.T, wvals) + self.hxb - self.yoblist).T))
        gradcost = obcost + wvals  # + bnd_cost
        return gradcost

    def find_min_tnc_ens_inc(self, dispp=5, f_tol=1e-4):
        """Function which minimizes 4DVAR cost fn. Takes an initial state
        (pvals).
        """
        wvals = self.xvals2wvals(self.xb)
        find_min = spop.fmin_tnc(self.cost_ens_inc, wvals, fprime=self.gradcost_ens_inc, disp=dispp, ftol=f_tol)
        xa = self.wvals2xvals(find_min[0])
        return find_min, xa

    def x_err(self, truth, x):
        return np.array([100*(abs(truth[i] - x[i])/truth[i]) for i in xrange(len(truth))])

    def da_run(self):
        self.jules_class.output_nml.mapping["JULES_OUTPUT_1_output_dir"] = "'"+self.output_dir+"',"
        #self.output_dir = out_dir
        if not os.path.exists(self.output_dir):
            os.makedirs(self.output_dir)
        res = self.minimize(self.xb)
        output = open('da_out_'+str(self.lat)+'_'+str(self.lon)+'.csv', 'w')
        for item in self.steps:
            output.write(str(item).strip("[]") + "\n")
        output.close()
        sh.rmtree(self.output_dir)
        print res.x


def ens_run_setup(ens_number, xbi):
    """
    Runs JULES for specified ens_member
    :param lat_lon: tuple containing latitude and longitude coordinate
    :return: na
    """
    out_dir = 'output_ens_' + str(ens_number) + '/'
    run_file = '/home/if910917/qSub_runMe/output_ens_' + str(ens_number) + '_' +\
               os.getcwd()[-2:] + '.bash'
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    for file in glob.glob(r'*.nml'):
        sh.copy(file, out_dir)
    os.chdir(out_dir)
    try:
        os.makedirs('output')
    except OSError:
        print 'already output dir'
    try:
        Jules_DA()
        lines = []
        lines.append('cd '+os.getcwd()+'\n')
        lines.append('module load python/canopy\n')
        lines.append('python ../fourdenvar.py ens_run '+str(xbi[0])+' '+str(xbi[1])+'\n')
        f = open(run_file, 'w')
        for line in lines:
            f.write(line)
        f.close()
        os.chdir('../')
    except ValueError:
        os.chdir('../')
        sh.rmtree(out_dir)
        print 'No data at:', lat_lon


if __name__ == "__main__":
    if sys.argv[1] == 'setup_ens':
        jcda = Jules_DA()
        # lats = np.array([9.75])
        # lons = np.array([0.75])
        for xbi in enumerate(jcda.xbs):
            ens_run_setup(xbi[0], xbi[1])
    elif sys.argv[1] == 'ens_run':
        jcda = Jules_DA()
        jcda.jules_class.model_grid_nml.mapping["jules_land_frac_1_file"] = " '../WFD-EI-LandFraction2d-2.nc',"
        jcda.jules_class.model_grid_nml.mapping["jules_latlon_1_file"] = " '../WFDEI-long-lat-2d-2.nc',"
        jcda.jules_class.drive_nml.mapping["jules_drive_1_file"] = " '../drive-file-WFDEI-79-12.txt',"
        jcda.jules_class.initial_conditions_nml.mapping["JULES_INITIAL_1_file"] = " '../output/gl4.dump.spin3.19830101.0.nc',"
        jcda.jules_class.ancillaries_nml.mapping["jules_frac_1_file"] = "'../land_frac_wheat100.nc',"
        jcda.run_jules(float(sys.argv[2]), float(sys.argv[3]))
    elif sys.argv[1] == 'run_da':
        print 'reading ensemble'
        update_soil_nc_all('soil.regional.nc', 'lonlat.regional.nc')
        os.chdir('run_forecast')
        subprocess.call(['python', 'py_jules_run.py'])
        print 'forecast finished'