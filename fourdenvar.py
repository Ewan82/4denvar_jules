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
        self.xb = np.array([4.2e-4, 5.2])
        self.b_mat = np.eye(2)*(self.x_truth*0.25)**2  # 20% background error
        self.xbs_10 = np.array([[4.56360650e-04, 4.64701202e+00],
                             [3.55562411e-04, 5.96723538e+00],
                             [2.99631089e-04, 5.69409886e+00],
                             [4.44187488e-04, 5.79117784e+00],
                             [3.96493282e-04, 4.71205087e+00],
                             [3.44420306e-04, 3.15449793e+00],
                             [5.22947407e-04, 3.84212025e+00],
                             [2.11949136e-04, 2.47367972e+00],
                             [4.02132357e-04, 5.82215200e+00],
                             [3.88853493e-04, 6.61392364e+00]])
        self.xbs = np.array([[3.55623130e-04, 3.14778039e+00],
                             [4.83270080e-04, 8.49176910e+00],
                             [6.22532291e-04, 7.56028227e+00],
                             [1.99819129e-04, 3.32287783e+00],
                             [5.65031187e-04, 6.24845614e+00],
                             [4.71011968e-04, 5.78745510e+00],
                             [5.89896466e-04, 5.91959661e+00],
                             [2.63734454e-04, 4.13578686e+00],
                             [4.92816911e-04, 6.47999844e+00],
                             [5.64322728e-04, 7.52271069e+00],
                             [2.88931053e-04, 3.26590862e+00],
                             [4.22701167e-04, 4.78879962e+00],
                             [5.49437392e-04, 2.96928816e+00],
                             [4.64534190e-04, 5.53094021e+00],
                             [2.73524070e-04, 7.88589023e+00]])

  # array of xbi's drawn from background cov mat
        self.size_ens = 15
        # abs(np.random.multivariate_normal(self.xb, self.b_mat))f

    def make_obs(self):
        # extract mod obs
        self.npp_obs, self.npp_ob_position, self.sm_obs, self.sm_ob_position = twin_obs.extract_data()
        self.npp_err = np.ones(len(self.npp_obs))*np.mean(self.npp_obs[self.npp_obs > 0])*0.1  # 10% error in mod obs
        self.sm_err = np.ones(len(self.sm_obs))*np.mean(self.sm_obs)*0.1  # 10% error in mod obs
        self.yoblist = np.append(self.npp_obs, self.sm_obs)
        self.yerr = np.append(self.npp_err, self.sm_err)
        self.rmatrix = np.eye(len(self.yerr))*self.yerr**2

    def make_hxb(self):
        self.npp_xb, self.npp_ob_position, self.sm_xb, self.sm_ob_position = twin_obs.extract_data(
            'output/background/xb.daily.nc')
        self.hxb = np.append(self.npp_xb, self.sm_xb)
        self.xb_mat = (1./(np.sqrt(self.size_ens-1)))*np.array([xbi - self.xb for xbi in self.xbs])
        self.xb_mat_inv = np.linalg.pinv(self.xb_mat.T)

    def run_jules(self, neff, b, run_id='ens'):
        """
        Runs JULES changing soil parameters
        :param run_id: id of run as a string
        :return: location of JULES output as string
        """
        self.jules_class.output_nml.mapping["jules_output_1_run_id"] = "'" + run_id + "',"
        self.jules_class.output_nml.mapping["jules_output_1_output_dir"] = self.output_dir
        self.jules_class.timesteps_nml.mapping["jules_time_1_main_run_start"] = " '2012-01-01 00:00:00',"
        self.jules_class.timesteps_nml.mapping["jules_spinup_1_max_spinup_cycles"] = " 0"
        self.jules_class.ancillaries_nml.mapping["jules_soil_props_1_const_val"] = str(b) + ", 0.3967309, 0.0027729999, 0.45809999, " \
                                                                    "0.3283205, 0.1866, 1185786.0, 0.2269195, 0.17,"
        self.jules_class.pft_params_nml.mapping["jules_pftparm_1_neff_io"] = "8.00e-4,8.00e-4,8.00e-4,4.00e-4,8.00e-4," + str(neff) + \
                                                              ", 8.00e-4,4.00e-4,8.00e-4,"
        self.jules_class.runJules()
        return self.output_dir + "/" + run_id + '.outvars.nc'

    def create_ensemble(self):
        #for xbi in enumerate(self.xbs):
        #    self.run_jules(xbi[1][0], xbi[1][1], 'ens'+str(xbi[0]))
        dumps = glob.glob('output/ens_run_15/ens*.dump*')
        for f in dumps:
            os.remove(f)
        hm_xbs = []
        for xb_fname in glob.glob('output/ens_run_15/ens*.daily.nc'):
            npp, npp_obp, sm, sm_obp = twin_obs.extract_data(xb_fname)
            hxbi = np.append(npp, sm)
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