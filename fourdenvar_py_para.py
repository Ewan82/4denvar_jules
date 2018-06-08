import numpy as np
import netCDF4 as nc
import datetime as dt
import scipy.optimize as spop
import scipy.stats as sps
import scipy.linalg as splinal
import multiprocessing as mp
import matplotlib.pyplot as plt
import itertools as itt
import jules
import run_jules_da as rjda
import os
import shutil as sh
import twin_obs
import glob
import subprocess
import sys
import pickle


class JulesTwinDA:
    def __init__(self, params='sow_date,alpha_io', assim=0):
        self.params = params
        self.p_keys = params.split(',')
        # python jules wrapper class
        # Set output dirctory
        self.output_dir = "output/ens_run_150"
        self.xb_outdir = "output/background"
        self.xi_outdir = "output/xi_run"
        # mod_truth and background
        self.p_dict_truth = {'sow_date': 110, 'sen_dvi': 0.4, 'gamma_io': 17.6, 'delta_io': -0.33, 'neff': 5.7e-4,
                       'alpha_io': 0.055, 'nl0_io': 0.07, 'f0_io': 0.4, 'dq_crit_io': 0.075, 'tt_rep': 908.4,
                       'tt_veg': 846.2, 'tt_emr_io': 108.2}
        self.p_dict_prior = {'sow_date': 99.0, 'sen_dvi': 0.44, 'gamma_io': 17.6, 'delta_io': -0.33, 'neff': 5.7e-4,
                       'alpha_io': 0.044, 'nl0_io': 0.062, 'f0_io': 0.4, 'dq_crit_io': 0.075, 'tt_rep': 908.4,
                       'tt_veg': 846.2, 'tt_emr_io': 98.0}
        self.p_dict = {'sow_date': 110, 'sen_dvi': 0.4, 'gamma_io': 17.6, 'delta_io': -0.33, 'neff': 5.7e-4,
                       'alpha_io': 0.055, 'nl0_io': 0.07, 'f0_io': 0.4, 'dq_crit_io': 0.075, 'tt_rep': 908.4,
                       'tt_veg': 846.2, 'tt_emr_io': 108.2}
        self.p_bnds = {'sow_date': (90.0, 130.0), 'sen_dvi': (0.1, 1.0), 'gamma_io': (12.0, 25.0),
                       'delta_io': (-1.0, -0.1), 'neff': (0., 1.0e-2), 'alpha_io': (0, 0.15), 'nl0_io': (0, 0.15),
                       'f0_io': (0., 1.0), 'dq_crit_io': (0, 0.1), 'tt_rep': (800., 1400.), 'tt_veg': (650., 950.),
                       'tt_emr_io': (50., 130.)}
        self.val_bnds = [self.p_bnds[key] for key in self.p_keys]
        self.xbs = pickle.load(open('/home/if910917/projects/4denvar_jules/xbs_150.p', 'rb'))
        self.xb = np.mean(self.xbs, axis=0)  # np.array([float(self.p_dict_prior[key]) for key in self.p_keys])
        self.x_true = np.array([float(self.p_dict_truth[key]) for key in self.p_keys])
        self.size_ens = 150
        self.xb_sd = self.xb*0.4
        self.b_mat = np.eye(len(self.xb))*((self.xb_sd)**2)  # 20% background error
        self.make_twin_obs()
        if assim == 1:
            self.make_twin_hxb()
            self.create_twin_ensemble()
        # self.xbs=pickle.load(open('/home/if910917/flux_net_envar_jules/xbs.p', 'rb'))
        # array of xbi's drawn from background cov mat
        # abs(np.random.multivariate_normal(self.xb, self.b_mat))f
        #self.make_obs()
        #self.make_hxb()
        #self.create_ensemble()

    def set_p_dict(self, x_arr):
        for i in xrange(len(self.p_keys)):
            self.p_dict[self.p_keys[i]] = x_arr[i]

    def make_twin_obs(self):
        # extract mod obs
        self.gpp_obs, self.gpp_ob_position, self.lai_obs, self.lai_ob_position = twin_obs.extract_twin_data()
        self.gpp_err = np.ones(len(self.gpp_obs))*0.5  # *np.mean(self.gpp_obs[self.gpp_obs > 0])*0.1  # 10% error in mod obs
        self.lai_err = np.ones(len(self.lai_obs))*0.05  # *np.mean(self.lai_obs)*0.05  # 10% error in mod obs
        self.yoblist = np.append(self.gpp_obs, self.lai_obs)
        self.yerr = np.append(self.gpp_err, self.lai_err)
        self.rmatrix = np.eye(len(self.yerr))*self.yerr**2

    def make_twin_hxb(self):
        rj = rjda.RunJulesDa(self.params, values=self.xb, nml_dir='twin_xb_nml/')
        rj.run_jules_dic(output_name='xb', out_dir='../output/background')
        self.gpp_xb, self.gpp_ob_position, self.lai_xb, self.lai_ob_position = twin_obs.extract_twin_data(
            'output/background/xb.daily.nc')
        self.hxb = np.append(self.gpp_xb, self.lai_xb)
        self.xb_mat = (1./(np.sqrt(self.size_ens-1)))*np.array([xbi - self.xb for xbi in self.xbs])
        self.xb_mat_inv = np.linalg.pinv(self.xb_mat.T)

    def create_twin_ensemble(self):
        self.hm_xbs = []
        for xb_fname in glob.glob(self.output_dir+'/ens*.daily.nc'):
            xbi_nc = nc.Dataset(xb_fname, 'r')
            gpp_xi = 1000 * 60 * 60 * 24 * xbi_nc.variables['gpp'][self.gpp_ob_position, 7, 0, 0]
            lai_xi = xbi_nc.variables['lai'][self.lai_ob_position, 7, 0, 0]
            hxbi = np.append(gpp_xi, lai_xi)
            self.hm_xbs.append(hxbi)
            xbi_nc.close()
        self.hmx_mat = (1. / (np.sqrt(self.size_ens - 1))) * np.array([hmxb - self.hxb for hmxb in self.hm_xbs])

    def make_obs(self):
        self.gpp_obs, self.gpp_ob_position, self.gpp_err = twin_obs.extract_gpp_day_data()

        self.le_obs, self.le_ob_position, self.le_err, self.h_obs, self.h_ob_position, self.h_err, self.g_obs, \
        self.g_ob_position = twin_obs.extract_energy_day_data()

        self.sm_obs, self.sm_ob_position = twin_obs.extract_sm_day_data()
        self.yoblist = self.gpp_obs
        self.yerr = np.ones((len(self.yoblist)))*0.30  # 0.03*self.yoblist
        self.rmatrix = np.eye(len(self.yerr))*self.yerr**2

    def make_hxb(self):
        self.run_jules_dic(self.xb_outdir, run_id='xb_exp')
        xb_nc = nc.Dataset(self.xb_outdir+'/xb_exp.daily.nc', 'r')
        self.gpp_xb = 1000 * 60 * 60 * 24 * xb_nc.variables['gpp_gb'][self.gpp_ob_position, 0, 0]
        self.le_xb = xb_nc.variables['latent_heat'][self.le_ob_position, 0, 0]
        self.he_xb = xb_nc.variables['ftl_gb'][self.h_ob_position, 0, 0]
        self.sm_xb = xb_nc.variables['smcl'][self.sm_ob_position, 0, 0, 0]/100.
        self.hxb = 1000 * 60 * 60 * 24 * xb_nc.variables['gpp'][self.gpp_ob_position, 1, 0, 0]
        self.xb_mat = (1./(np.sqrt(self.size_ens-1)))*np.array([xbi - self.xb for xbi in self.xbs])
        self.xb_mat_inv = np.linalg.pinv(self.xb_mat.T)

    def make_hxi(self, x_arr):
        print 'cost', x_arr
        self.set_p_dict(x_arr)
        self.run_jules_dic(self.xi_outdir, run_id='xi_test')
        xi_nc = nc.Dataset(self.xi_outdir+'/xi_test.daily.nc', 'r')
        self.hxi = 1000 * 60 * 60 * 24 * xi_nc.variables['gpp'][self.gpp_ob_position, 1, 0, 0]
        self.hmx_mat = (1. / (np.sqrt(self.size_ens - 1))) * np.array([hmxb - self.hxi for hmxb in self.hm_xbs])

    def create_ensemble(self):
        #for xbi in enumerate(self.xbs):
        #    self.run_jules(lai=xbi[1][0], tlow=xbi[1][1], tupp=xbi[1][2], run_id='ens'+str(xbi[0]))
        #dumps = glob.glob(self.output_dir+'/ens*.dump*')
        #for f in dumps:
        #    os.remove(f)
        self.hm_xbs = []
        for xb_fname in glob.glob(self.output_dir+'/ens*.daily.nc'):
            xbi_nc = nc.Dataset(xb_fname, 'r')
            hxbi = 1000 * 60 * 60 * 24 * xbi_nc.variables['gpp'][self.gpp_ob_position, 1, 0, 0]
            self.hm_xbs.append(hxbi)
            xbi_nc.close()
        self.hmx_mat = (1. / (np.sqrt(self.size_ens - 1))) * np.array([hmxb - self.hxb for hmxb in self.hm_xbs])

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

    def obcost_ens(self, wvals):
        """Observational part of cost fn.
        """
        xvals = self.wvals2xvals(wvals)
        self.make_hxi(xvals)
        return np.dot(np.dot((self.hxi - self.yoblist), np.linalg.inv(self.rmatrix)),
                      (self.hxi - self.yoblist).T)

    def cost_ens_inc(self, wvals):
        modcost = np.dot(wvals, wvals.T)
        obcost = self.obcost_ens_inc(wvals)
        ret_val = 0.5 * obcost + 0.5 * modcost
        return ret_val

    def cost_ens(self, wvals):
        modcost = np.dot(wvals, wvals.T)
        obcost = self.obcost_ens(wvals)
        ret_val = 0.5 * obcost + 0.5 * modcost
        return ret_val

    def gradcost_ens_inc(self, wvals):
        obcost = np.dot(self.hmx_mat, np.dot(np.linalg.inv(self.rmatrix),
                                           (np.dot(self.hmx_mat.T, wvals) + self.hxb - self.yoblist).T))
        gradcost = obcost + wvals  # + bnd_cost
        return gradcost

    def gradcost_ens(self, wvals):
        xvals = self.wvals2xvals(wvals)
        print 'grad', xvals
        obcost = np.dot(self.hmx_mat, np.dot(np.linalg.inv(self.rmatrix),
                                           (self.hxi - self.yoblist).T))
        gradcost = obcost + wvals  # + bnd_cost
        return gradcost

    def a_cov(self):
        a_cov = np.linalg.inv(splinal.sqrtm(np.eye(self.size_ens) + np.dot(self.hmx_mat,
                        np.dot(np.linalg.inv(self.rmatrix), self.hmx_mat.T))))
        return a_cov

    def a_ens(self, xa):
        a_cov = self.a_cov()
        xa_mat = np.dot(self.xb_mat.T, a_cov).T
        a_ens = np.array([xa + np.sqrt(self.size_ens-1)*xbi for xbi in xa_mat])
        return a_ens

    def find_min_tnc_ens_inc(self, dispp=5):
        """Function which minimizes 4DVAR cost fn. Takes an initial state
        (pvals).
        """
        wvals = self.xvals2wvals(self.xb)
        find_min = spop.fmin_tnc(self.cost_ens_inc, wvals, fprime=self.gradcost_ens_inc, disp=dispp)
        xa = self.wvals2xvals(find_min[0])
        return find_min, xa

    def find_min_tnc_ens(self, dispp=5):
        """Function which minimizes 4DVAR cost fn. Takes an initial state
        (pvals).
        """
        wvals = self.xvals2wvals(self.xb)
        find_min = spop.fmin_tnc(self.cost_ens, wvals, fprime=self.gradcost_ens, disp=dispp, maxfun=20000)
        xa = self.wvals2xvals(find_min[0])
        return find_min, xa

    def find_min_ens_inc(self, dispp=5):
        """Function which minimizes 4DVAR cost fn. Takes an initial state
        (pvals).
        """
        wvals = self.xvals2wvals(self.xb)
        find_min = spop.fmin_ncg(self.cost_ens_inc, wvals, fprime=self.gradcost_ens_inc, disp=dispp, full_output=1,
                                 maxiter=20000)
        xa = self.wvals2xvals(find_min[0])
        return find_min, xa

    def find_min_glob_ens_inc(self, dispp=1):
        """Function which minimizes 4DVAR cost fn. Takes an initial state
        (pvals).
        """
        wvals = self.xvals2wvals(self.xb)
        find_min = spop.basinhopping(self.cost_ens_inc, wvals,
                                     minimizer_kwargs={'method': 'BFGS', 'jac': self.gradcost_ens_inc}, disp=dispp)
        xa = self.wvals2xvals(find_min.x)
        return find_min, xa

    def x_err(self, truth, x):
        return np.array([100*(abs(truth[i] - x[i])/truth[i]) for i in xrange(len(truth))])

    def test_bnds(self, xbi):
        for xi in enumerate(xbi):
            if self.val_bnds[xi[0]][0] < xi[1] < self.val_bnds[xi[0]][1]:
                continue
            else:
                return False
        return True

    def generate_param_ens(self, size_ens):
        ens = []
        i = 0
        while len(ens) < size_ens:
            xbi = np.random.multivariate_normal(self.xb, self.b_mat)
            i += 1
            if self.test_bnds(xbi) is True:
                ens.append(xbi)
            else:
                continue
        print i
        return np.array(ens)

    def perturb_true_state(self, x_truth, err=0.05):
        i = 0
        while i < 1:
            xbi = np.random.multivariate_normal(x_truth, (err*self.x_true)**2 * np.eye(len(self.x_true)))
            if self.test_bnds(xbi) is True:
                i += 1
            else:
                continue
        return xbi


def write_qsub():
    run_file = '/home/if910917/qSub_runMe/jules_da_py_para.bash'
    lines = []
    lines.append('cd '+os.getcwd()+'\n')
    lines.append('module load python/canopy-2.1.3\n')
    lines.append('python fourdenvar_py_para.py\n')
    f = open(run_file, 'w')
    for line in lines:
        f.write(line)
    f.close()
    return 'done'


def ens_member_run(ens_number_xbi):
    """
    Runs JULES for specified ens_member
    :param lat_lon: tuple containing latitude and longitude coordinate
    :return: na
    """
    ens_number = ens_number_xbi[0]
    xbi = ens_number_xbi[1]
    out_dir = 'output_ens_' + str(ens_number) + '/'
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    for file in glob.glob(r'twin_xb_nml/*.nml'):
        sh.copy(file, out_dir)
    try:
        rj = rjda.RunJulesDa(values=xbi, nml_dir=out_dir)
        rj.run_jules_dic(output_name='ens'+str(ens_number), out_dir="../output/ens_run_150")
        dumps = glob.glob('/home/if910917/projects/4denvar_jules/output/ens_run_150/ens'+str(ens_number)+'.dump*')
        for f in dumps:
            os.remove(f)
        sh.rmtree(out_dir)
    except ValueError:
        sh.rmtree(out_dir)
        print 'Something went wront at: ' + str(ens_number)


def ens_run(x_ens):
    print 'Running ensemble'
    mp.freeze_support()
    pool = mp.Pool()
    res = pool.map(ens_member_run, enumerate(x_ens))
    print res
    pool.close()
    pool.join()


def array2str(arr):
    arr_str = ''
    for val in arr:
        arr_str += str(val)+','
    arr_str = arr_str.strip(',')
    return arr_str


def get_truncated_normal(mean=0, sd=1, low=0, upp=10):
    return sps.truncnorm((low - mean) / sd, (upp - mean) / sd, loc=mean, scale=sd)


def plot_gpp(xb_nc, xa_nc, obs=None, ob_pos=None, key='gpp'):
    datxb = nc.Dataset(xb_nc, 'r')
    datxa = nc.Dataset(xa_nc, 'r')
    times = nc.num2date(datxa.variables['time'][:], datxa.variables['time'].units)
    if key == 'gpp':
        plt.plot(1000 * 60 * 60 * 24 * datxb.variables[key][:, 7, 0, 0])
        plt.plot(1000 * 60 * 60 * 24 * datxa.variables[key][:, 7, 0, 0])
    else:
        plt.plot(datxb.variables[key][:, 7, 0, 0])
        plt.plot(datxa.variables[key][:, 7, 0, 0])
    if obs is not None:
        plt.plot(ob_pos, obs, 'o')
    #plt.plot(times[jcda.gpp_ob_position], jcda.gpp_obs, 'o')
    plt.show()


if __name__ == "__main__":
    jcda = JulesTwinDA()
    print jcda.xbs
    ens_run(jcda.xbs)
    print 'ensemble has been run'
