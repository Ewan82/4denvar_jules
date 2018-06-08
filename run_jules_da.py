import jules


class RunJulesDa:
    def __init__(self, params='sow_date,neff', values='default', nml_dir='example_nml',
                 year=2010):
        self.p_keys = params.split(',')
        self.nml_dir = nml_dir
        self.year = year
        # mod_truth and background
        self.pft_idx = 7
        self.p_dict = {'sow_date': 110, 'sen_dvi': 0.4, 'gamma_io': 17.6, 'delta_io': -0.33, 'neff': 5.7e-4,
                       'alpha_io': 0.055, 'nl0_io': 0.07, 'f0_io': 0.4, 'dq_crit_io': 0.075, 'tt_rep': 908.4,
                       'tt_veg': 846.2, 'tt_emr_io': 108.2}
        if values is not 'default':
            for i in xrange(len(self.p_keys)):
                print values[i]
                self.p_dict[self.p_keys[i]] = values[i]

    def run_jules_dic(self, output_name='test', out_dir='../output/test'):
        """
        Function that runs JULES with crop model turned on and given user defined parameters. Output is saved in folder
        and file specified within function.

        :param output_name: Name to use for outputted JULES netCDF file.
        :type output_name: str.
        :param output_dir: Directory for writing JULES output.
        :type output_dir: str.
        :return: 'Done' to notify used JULES run has finished.
        :rtype: str
        """
        j = jules.Jules(self.nml_dir)
        # j.drive_nml.mapping['file']='path/to/your/drivers/metData'+n+'.dat'  # unnecessary here as using WFD for jules
        j.nml_dic['output']['jules_output']['run_id'] = output_name
        j.nml_dic['output']['jules_output']['output_dir'] = out_dir
        j.nml_dic['timesteps']['jules_time']['main_run_start'] = str(self.year) + '-01-01 00:00:00'
        j.nml_dic['timesteps']['jules_time']['main_run_end'] = str(self.year) + '-12-31 23:30:00'
        j.nml_dic['timesteps']['jules_spinup']['max_spinup_cycles'] = 4
        # j.nml_dic['ancillaries']['jules_soil_props']['const_val'][0] = self.p_dict['b']
        # j.nml_dic['ancillaries']['jules_soil_props']['const_val'][5] = self.p_dict['sm_wilt']
        j.nml_dic['pft_params']['jules_pftparm']['neff_io'][self.pft_idx] = self.p_dict['neff']
        j.nml_dic['pft_params']['jules_pftparm']['alpha_io'][self.pft_idx] = self.p_dict['alpha_io']
        j.nml_dic['pft_params']['jules_pftparm']['nl0_io'][self.pft_idx] = self.p_dict['nl0_io']
        j.nml_dic['pft_params']['jules_pftparm']['f0_io'][self.pft_idx] = self.p_dict['f0_io']
        j.nml_dic['pft_params']['jules_pftparm']['dqcrit_io'][self.pft_idx] = self.p_dict['dq_crit_io']
        if self.pft_idx > 4:
            j.nml_dic['crop_params']['jules_cropparm']['sen_dvi_io'][self.pft_idx - 5] = self.p_dict['sen_dvi']
            j.nml_dic['crop_params']['jules_cropparm']['gamma_io'][self.pft_idx - 5] = self.p_dict['gamma_io']
            j.nml_dic['crop_params']['jules_cropparm']['delta_io'][self.pft_idx - 5] = self.p_dict['delta_io']
            j.nml_dic['crop_params']['jules_cropparm']['tt_emr_io'][self.pft_idx - 5] = self.p_dict['tt_emr_io']
            j.nml_dic['ancillaries']['jules_crop_props']['const_val'][0] = self.p_dict['tt_veg']
            j.nml_dic['ancillaries']['jules_crop_props']['const_val'][1] = self.p_dict['tt_rep']
            j.nml_dic['ancillaries']['jules_crop_props']['const_val'][2] = self.p_dict['sow_date']
        j.run_jules_print()
        return output_name, 'done'