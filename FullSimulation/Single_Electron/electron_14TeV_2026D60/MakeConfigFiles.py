#!/usr/bin/env python3

import os
import subprocess

class RunInstance:
    def __init__(self):
        self.E_list = [50, 100, 250, 500]
        self.eta_list = [0., 0.6, 1.2, 1.8, 2.4, 3.0, 3.6]
        self.nevents = 500
        self.top_save_dir = os.getcwd()
        self.dir_run = os.path.abspath(__file__ + '/../run/')
        
        self.preformat_cfgfile = '{dir_step}/step{step}_2026D60_14TeV_electron_E{E}_eta{eta}{suffix}_cfg.py'
        self.preformat_annotation = '\'step{save_dir}_2026D60_14TeV_electron_E{step}_eta:{E} nevts:{eta}\''
        self.preformat_rootfile = '\'file:{save_dir}/step{step}_electron_E{E}_eta{eta}{suffix}.root\''
        self.preformat_aging = ''
        self.aging_suffix = ''
    
    def set_E (self, *args):
        """ Sets list of E values that will be passed on to the generator.
        
        Parameters
        ----------
        *args
            Variable length list of E values.
            If any of the values entered is not integer, it will raise an exception.
        
        Returns
        -------
        (none)
            If no value is provided in args, it will be set to default.
            Otherwise, it will save the values in the order they were provided.
        
        """
        
        if any(type(arg)!=int for arg in args):
            raise Exception("Please provide int only in args.")
        elif len(args) > 0:
            self.E_list = args


    def set_eta (self, *args):
        """ Sets list of eta values that will be passed on to the generator.

        Parameters
        ----------
        *args
	        Variable length list of eta values.
	        If any of the values entered is not int or float, it will raise an exception.

        Returns
        -------
        (none)
	        If no value is provided in args, it will be set to default.
            Otherwise, it will save the values in the order they were provided.

        """

        if any(type(arg)!=float and type(arg)!=int for arg in args):
            raise Exception("Please provide int or float only in args.")
        elif len(args) != 0:
            self.eta_list = args
                
                
    def set_nevents (self, nevents):
        """ Sets the number of events to generate.
        
        Parameters
        ----------
        nevents
            Number of events in integer.
        
        Returns
        -------
        (none)
            If arg is set to 0 or below, it will be set to default.
        
        """
        
        if type(nevents) != int:
            raise Exception("Please provide int only in arg.")
        elif nevents > 0:
            self.nevents = nevents


    def set_top_save_dir (self, dir):
        """ Set and create the top directory to save the files.
        
        Parameters
        ----------
        dir
            Name of the top saving directory in string.
        
        Returns
        -------
        (none)
            If arg is an empty string, it will be set to default.
        
        """
            
        if type(dir) != str:
            raise Exception("Please provide string only in arg.")
        elif len(dir) > 0:
            self.top_save_dir = dir
            if not os.path.exists(self.top_save_dir):
                os.makedirs(self.top_save_dir)

        self.top_save_dir = self.top_save_dir + self.aging_suffix
                
                
    def set_aging (self, lumi, tag=''):
        """ Set custom aging parameters.
        
        Parameters
        ----------
        lumi
            Luminosity in inverse femtobarns.
        
        tag
            Tag of custom aging process.
        
        Returns
        -------
        (none)
        """
        aging_template = ""
        
        if lumi == 0:
            return
        elif lumi == 3000:
            if tag == '':
                aging_template = "custom_aging_3000.txt"
            elif tag == 'ultimate':
                aging_template = "custom_aging_3000_ultimate.txt"
        elif lumi == 4500:
            if tag == '':
                aging_template = "custom_aging_4500.txt"
                
        if aging_template == "":
            raise Exception('Please provide valid arguments.')
        else:
            if tag != '':
                tag = '_' + tag
            self.aging_suffix = '_aging_{lumi}{tag}'.format(lumi=lumi, tag=tag)
            with open('cfg_templates/{}'.format(aging_template), 'r') as template_file:
                self.preformat_aging = template_file.read()
            
            
    def makeAllConfigFiles (self):
        """ Create config files
        
        Parameters
        ----------
        (none)
        
        Returns
        -------
        None
              
        """
        self.makeStep1ConfigFiles()
        self.makeStep2ConfigFiles()
        self.makeStep3ConfigFiles()
        
        return None

        
    def runSteps (self, *steps):
        """ Run the step cfg.py files for the step numbers specified.
        
        All config files in the same steps will run in parallel.
        
        Parameters
        ----------
        
        *steps: int or str
            Steps to run. It will automatically be reordered in increasing value.
            
        Returns
        -------
        None
        
        """
#        
#        if not os.path.exists(self.dir_run):
#            os.makedirs(self.dir_run)
            
        steps = sorted(s for s in steps)
        
        for step in steps:
            cmd_list = []
            dir_step = self.dir_run + '/step{}_config'.format(step)
            for E in self.E_list:
                for eta in self.eta_list:
                    cmd_list.append("cmsRun " + self.preformat_cfgfile.format(dir_step=dir_step, step=step, suffix=self.aging_suffix, E=E, eta=eta))
            dir_cfg = os.path.abspath(self.dir_run + "/step{}_config".format(step))
            bash_command = " & ".join(cmd_list)
            p = subprocess.Popen(bash_command, shell=True)
            p.wait()
            
        return None
        
        
    def purge (self):
        """ Clear the run directory of all config files.
        
        Parameters
        ----------
        (none)
        
        Returns
        -------
        None
        
        """
        
        subprocess.run("rm {}/*".format(self.dir_run))
        
        return None


    def makeStep1ConfigFiles (self):
        """ Generates CMSSW cfg python scripts to run step 1 generation.
        
        Parameters
        ----------
        None
            
        Returns
        -------
        None
        """
        
        # What to write in file
        with open('cfg_templates/step1_template.txt', 'r') as template_file:
            template_text = template_file.read()
    
        # Set directory to put cfg.py files
        dir_step1 = self.dir_run + '/step1_config'
        if not os.path.exists(self.dir_run + '/step1_config'):
            os.makedirs(dir_step1)
        
        # Save the cfg file for all E and eta values
        for E in self.E_list:
            for eta in self.eta_list:
            
                preformat_dict = dict()
                preformat_dict['save_dir'] = self.top_save_dir + '/electron_E{0}_eta{1}'.format(E, eta)
                preformat_dict['E'] = E
                preformat_dict['eta'] = eta

                # Content of cfg file
                format_dict = dict()
                format_dict['maxevents'] = self.nevents
                format_dict['annotation'] = self.preformat_annotation.format(step=1, **preformat_dict)
                format_dict['outfile'] = self.preformat_rootfile.format(step=1, suffix=self.aging_suffix + '', **preformat_dict)
                format_dict['maxEta'] = eta + 0.001
                format_dict['minEta'] = eta - 0.001
                format_dict['maxE'] = E + 0.001
                format_dict['minE'] = E - 0.001
                format_dict['psethack'] = '\'single electron E {0} eta {1}\''.format(E, eta)
                template_text = template_text.format(**format_dict)
                
                # Name of cfg file
                cfg_filename = self.preformat_cfgfile.format(dir_step=dir_step1, step=1, suffix=self.aging_suffix, **preformat_dict)

                # Set directory to save root files
                if not os.path.exists(preformat_dict['save_dir']):
                    os.makedirs(preformat_dict['save_dir'])
                
                with open(cfg_filename, 'w') as cfg_file:
                    cfg_file.write(template_text)
        
        return

        

    def makeStep2ConfigFiles (self):
        """ Generates CMSSW cfg python files to run step 2.
        
        Parameters
        ----------
        None
            
        Returns
        -------
        None
        """
        
        # What to write in file
        with open('cfg_templates/step2_template.txt', 'r') as template_file:
            template_text = template_file.read()
    
        # Set directory to put cfg.py files
        dir_step2 = self.dir_run + '/step2_config'
        if not os.path.exists(self.dir_run + '/step2_config'):
            os.makedirs(dir_step2)
        
        # Save the cfg file for all E and eta values
        for E in self.E_list:
            for eta in self.eta_list:
                    
                preformat_dict = dict()
                preformat_dict['save_dir'] = self.top_save_dir + '/electron_E{0}_eta{1}'.format(E, eta)
                preformat_dict['E'] = E
                preformat_dict['eta'] = eta
            
                # Content of cfg file
                format_dict = dict()
                format_dict['maxevents'] = self.nevents
                format_dict['inputfile'] = self.preformat_rootfile.format(step=1, suffix=self.aging_suffix + '', **preformat_dict)
                format_dict['annotation'] = self.preformat_annotation.format(step=2, **preformat_dict)
                format_dict['outfile'] = self.preformat_rootfile.format(step=2, suffix=self.aging_suffix + '', **preformat_dict)
                format_dict['custom_aging'] = self.preformat_aging
                
                template_text = template_text.format(**format_dict)
                
                # Name of cfg file
                cfg_filename = self.preformat_cfgfile.format(dir_step=dir_step2, step=2, suffix=self.aging_suffix, **preformat_dict)

                # Set directory to save root files
                if not os.path.exists(preformat_dict['save_dir']):
                    os.makedirs(preformat_dict['save_dir'])
                    
                with open(cfg_filename, 'w') as cfg_file:
                    cfg_file.write(template_text)
        
        return


    def makeStep3ConfigFiles (self):
        """ Generates CMSSW cfg python scriEs to run step 3 reconstruction.
        
        If there are N E values given, N files will be generated in the sub-directory run/step3_config.
        
        Parameters
        ----------
        None
            
        Returns
        -------
        None
        """
        
        # What to write in file
        with open('cfg_templates/step3_template.txt', 'r') as template_file:
            template_text = template_file.read()
    
        # Set directory to put cfg.py files
        dir_step3 = self.dir_run + '/step3_config'
        if not os.path.exists(self.dir_run + '/step3_config'):
            os.makedirs(dir_step3)
        
        # Save the cfg file for all E and eta values
        for E in self.E_list:
            for eta in self.eta_list:
            
                preformat_dict = dict()
                preformat_dict['save_dir'] = self.top_save_dir + '/electron_E{0}_eta{1}'.format(E, eta)
                preformat_dict['E'] = E
                preformat_dict['eta'] = eta
            
                # Content of cfg file
                format_dict = dict()
                
                format_dict['maxevents'] = self.nevents
                format_dict['inputfile'] = self.preformat_rootfile.format(step=2, suffix=self.aging_suffix + '', **preformat_dict)
                format_dict['annotation'] = self.preformat_annotation.format(step=3, **preformat_dict)
                format_dict['outfile'] = self.preformat_rootfile.format(step=3, suffix=self.aging_suffix + '', **preformat_dict)
                format_dict['outfile_miniaodsim'] = self.preformat_rootfile.format(step=3, suffix=self.aging_suffix + '_inMINIAODSIM', **preformat_dict)
                format_dict['outfile_dqm'] = self.preformat_rootfile.format(step=3, suffix=self.aging_suffix + '_inDQM', **preformat_dict)
                
                template_text = template_text.format(**format_dict)
                
                # Name of cfg file
                cfg_filename = self.preformat_cfgfile.format(dir_step=dir_step3, step=3, suffix=self.aging_suffix, **preformat_dict)

                # Set directory to save root files
                if not os.path.exists(preformat_dict['save_dir']):
                    os.makedirs(preformat_dict['save_dir'])
                
                with open(cfg_filename, 'w') as cfg_file:
                    cfg_file.write(template_text)
                    
        return

