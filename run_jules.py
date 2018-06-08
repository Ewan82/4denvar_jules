import jules


def crop_run(sow_date=110, b=6.631, smwilt=0.1866, neff=4.0e-4, pft_idx=7, nml_dir='example_nml', output_name='test',
             out_dir='../output/test'):
    """
    Function that runs JULES with crop model turned on and given user defined parameters at Wallerfing site. Output is
    saved in folder and file specified within function.

    :param sow_date: Sow date, between 90 and 150.
    :type sow_date: int.
    :param b: Brooks-Corey exponent factor.
    :type b: float.
    :param smwilt: Soil moisture wilting point.
    :type smwilt: float.
    :param neff: Nitrogen use efficiency of crop (Vcmax).
    :type neff: float.
    :param pft_idx: index for pft parameters to update.
    :type pft_idx: int
    :param nml_dir: Location of JULES nml file directory.
    :type nml_dir: str.
    :param output_name: Name to use for outputted JULES netCDF file.
    :type output_name: str.
    :param output_dir: Directory for writing JULES output.
    :type output_dir: str.
    :return: 'Done' to notify used JULES run has finished.
    :rtype: str
    """
    j = jules.Jules(nml_dir)
    # j.drive_nml.mapping['file']='path/to/your/drivers/metData'+n+'.dat'  # unnecessary here as using WFD for jules
    j.nml_dic['output']['jules_output']['run_id'] = output_name
    j.nml_dic['output']['jules_output']['output_dir'] = out_dir
    year=2010
    j.nml_dic['timesteps']['jules_time']['main_run_start'] = str(year) + '-01-01 00:00:00'
    j.nml_dic['timesteps']['jules_time']['main_run_end'] = str(year) + '-12-31 23:30:00'
    j.nml_dic['timesteps']['jules_spinup']['max_spinup_cycles'] = 2
    j.nml_dic['ancillaries']['jules_crop_props']['const_val'][2] = sow_date
    j.nml_dic['ancillaries']['jules_soil_props']['const_val'][0] = b
    j.nml_dic['ancillaries']['jules_soil_props']['const_val'][5] = smwilt
    j.nml_dic['pft_params']['jules_pftparm']['neff_io'][pft_idx] = neff
    j.run_jules()
    return 'Done'


def plot_class_var(output_nc, var, level=0, line_type='-', ax='None'):
    """Plot specified variable.

    :param output_nc: Location of JULES output netCDF file.
    :type output_nc: str
    :param var: Variables from JULES to plot.
    :type var: str
    :return: Figure.
    :rtype: object
    """
    sns.set_context('poster', font_scale=1.2, rc={'lines.linewidth': 2, 'lines.markersize': 6})
    #ax.xaxis_date()
    sns.set_style('whitegrid')
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
