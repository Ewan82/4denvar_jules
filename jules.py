# !/home/db903833//dataLand01/enthoughtDistros/epd-7.2-2-rh5-x86_64/bin/python
# !/usr/bin/env python

# core python modules:
import subprocess
import os
# 3rd party modules:
import numpy as np
import sys
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import netCDF4 as nc
import f90nml
import glob


class Jules():
    """Class to read in nml files and run JULES.

    :param nml_dir: location of JULES nml file directory.
    :type nml_dir: str
    :param jules_exe: location of JULES executable.
    :type jules_exe: str

    .. note:: You must have JULES installed on local system with a version of 4.8 or higher.

    """
    def __init__(self, nml_dir, jules_exe='/home/if910917/jules/models/jules4.8/build/bin/jules.exe'):
        self.nml_dir = nml_dir
        self.jules = jules_exe
        self.nml_dic = self.read_nml(nml_dir)

    def read_nml(self, nml_dir):
        """
        Reading all nml files in specified directory and stores them in a dictionary.
        :param nml_dir: directory to read nml files from
        :return: dictionary of nml files
        """
        nml_dic = {}
        for f_nml in glob.glob(nml_dir + '/*.nml'):
            nml_dic[f_nml.split('/')[-1][:-4]] = f90nml.read(f_nml)
        return nml_dic

    def write_nml(self):
        """
        Function to write dictionary of stored nml data to nml files in current working dir.
        :return: n/a
        """
        for key in self.nml_dic.keys():
            self.nml_dic[key].write(key + '.nml', force=True)

    def run_jules(self):
        """Write all NML files to disk. Run JULES in a subprocess. Check output for fatal errors.

        :return: stdout and stderr output from JULES model run.
        :rtype: tuple
        """

        # write all the nml files here so the
        # user doesn't have to remember to...

        # run JULES
        dir_path = os.path.dirname(os.path.realpath(__file__))
        cwd = os.getcwd()
        os.chdir(dir_path+'/'+self.nml_dir)
        self.write_nml()
        cmd = []
        cmd.append(self.jules)

        p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        out = p.stdout.readlines()
        err = p.stderr.readlines()
        p.wait()
        os.chdir(cwd)

        # catch "fatal" errors
        for line in out:
            if len(line.split()) == 0: continue
            if line.split()[0] == "[FATAL":
                print >> sys.stderr, "*** runJules: caught fatal error in JULES run:"
                print >> sys.stderr, line,
                sys.exit()

        # catch anything in stderr
        if len(err) > 0:
            for line in err:
                print >> sys.stderr, "*** runJules: caught output on stderr in JULES run:"
                print >> sys.stderr, line,
                # sys.exit()
        return (out, err)

    def run_jules_print(self):
        """Write all NML files to disk. Run JULES in a subprocess. Check output for fatal errors.

        :return: stdout and stderr output from JULES model run.
        :rtype: str
        """

        # write all the nml files here so the
        # user doesn't have to remember to...
        dir_path = os.path.dirname(os.path.realpath(__file__))
        cwd = os.getcwd()
        os.chdir(dir_path+'/'+self.nml_dir)
        self.write_nml()

        # run JULES
        cmd = []
        cmd.append(self.jules)
        proc = subprocess.Popen(cmd, shell=False)
        proc.communicate()
        os.chdir(cwd)
        return 'Done', 'Done'


def plot_class_var(output_nc, var, level=0, line_type='-', ax='None'):
    """Plot specified variable.

    :param output_nc: Location of JULES output netCDF file.
    :type output_nc: str
    :param var: Variables from JULES to plot.
    :type var: str
    :return: Figure.
    :rtype: object
    """
    dat = nc.Dataset(output_nc, 'r')
    date_lst = nc.num2date(dat.variables['time'][:], dat.variables['time'].units)
    if len(dat.variables[var]) == 4:
        var_dat = dat.variables[var][:, level, 0, 0]
    else:
        var_dat = dat.variables[var][:, 0, 0]
    plt.plot(date_lst, var_dat, line_type)
    if var == 'croplai':
        plt.ylabel(r'Crop LAI (m$^2$ m$^{-2}$)')
    elif var == 'smcl':
        plt.ylabel(r'Soil Moisture (kg m$^{-2}$ s$^{-1}$)')
    elif var == 'cropcanht':
        plt.ylabel(r'Crop canopy height (m)')
    else:
        plt.ylabel(dat.variables[var].long_name + ' (' + dat.variables[var].units + ')')
    plt.xlabel('Date')
    plt.title('JULES output for Wallerfing')
    myFmt = mdates.DateFormatter('%B')
    plt.gca().xaxis.set_major_formatter(myFmt)
    plt.gcf().autofmt_xdate()
    # plt.xaxis.set_major_formatter(myFmt)
    # plt.legend(loc=2)
    # plt.show()
    return 'plot finished'


if __name__ == "__main__":

    # j = jules()
    # j.timesteps_nml.mapping["jules_time_1_main_run_start"] = " '2011-01-01 00:00:00',"
    # j.runJules()

    crop_run()
    for x in xrange(100):
        print x
        sow_date = int(np.random.normal(110, 0.1*110))
        b = np.random.normal(6.631, 0.1*6.631)
        smwilt = np.random.normal(0.1866, 0.1*0.1866)
        neff = np.random.normal(5.7e-4, 0.1*5.7e-4)
        crop_run(sow_date, b, smwilt, neff)
