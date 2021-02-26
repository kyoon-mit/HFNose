#!/usr/bin/env python2.7

import os
import subprocess

class RunInstance:
    def __init__(self):
        self.E_list = [50, 100, 250, 500]
        self.eta_list = [0., 0.6, 1.2, 1.8, 2.4]
        self.nevents = 500
        self.top_save_dir = os.getcwd()
        self.dir_run = os.path.abspath(__file__ + '/../run/')
        
        self.preformat_cfgfile = '{0}/step{1}_2026D60_HGCcenter_14TeV_photon_E{2}_eta{3}_cfg.py'
        self.preformat_annotation = '\'step{0}_2026D60_14TeV_photon_E{1}_eta:{2} nevts:{3}\''
        self.preformat_rootfile = '\'file:{0}/step{1}_photon_E{2}_eta{3}{4}.root\''  
    ''
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
	    If any of the values entered exceeds 3.0, it will raise an exception

        Returns
        -------
        (none)
	        If no value is provided in args, it will be set to default.
            Otherwise, it will save the values in the order they were provided.

        """

        if any(type(arg)!=float and type(arg)!=int for arg in args):
            raise Exception("Please provide int or float only in args.")
        elif any(arg>=3.0 for arg in args):
            raise Exception("Please provide eta values below 3.0.")
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
                    cmd_list.append("cmsRun " + self.preformat_cfgfile.format(dir_step, step, E, eta))
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
        if not os.path.exists(self.dir_run + '/step2_config'):
            os.makedirs(dir_step1)
        
        # Save the cfg file for all E and eta values
        for E in self.E_list:
            for eta in self.eta_list:
            
                # Content of cfg file
                save_dir = self.top_save_dir + '/photon_E{0}_eta{1}'.format(E, eta)
                maxevents = self.nevents
                annotation = self.preformat_annotation.format(1, E, eta, self.nevents)
                outfile = self.preformat_rootfile.format(save_dir, 1, E, eta, "")
                maxEta = eta + 0.001
                minEta = eta - 0.001
                maxE = E + 0.001
                minE = E - 0.001
                psethack = '\'single photon E {0} eta {1}\''.format(E, eta)
                template_text = template_text.format(maxevents, annotation, outfile, maxEta, minEta, maxE, minE, psethack)
                
                # Name of cfg file
                cfg_filename = self.preformat_cfgfile.format(dir_step1, 1, E, eta)

                # Set directory to save root files
                if not os.path.exists(save_dir):
                    os.makedirs(save_dir)
                
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
            
                # Content of cfg file
                save_dir = self.top_save_dir + '/photon_E{0}_eta{1}'.format(E, eta)
                maxevents = self.nevents
                inputfile = self.preformat_rootfile.format(save_dir, 1, E, eta, "")
                annotation = self.preformat_annotation.format(2, E, eta, self.nevents)
                outfile = self.preformat_rootfile.format(save_dir, 2, E, eta, "")
                template_text = template_text.format(maxevents, inputfile, annotation, outfile)
                
                # Name of cfg file
                cfg_filename = self.preformat_cfgfile.format(dir_step2, 2, E, eta)

                # Set directory to save root files
                if not os.path.exists(save_dir):
                    os.makedirs(save_dir)
                    
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
        if not os.path.exists(self.dir_run + '/step2_config'):
            os.makedirs(dir_step3)
        
        # Save the cfg file for all E and eta values
        for E in self.E_list:
            for eta in self.eta_list:
            
                # Content of cfg file
                save_dir = self.top_save_dir + '/photon_E{0}_eta{1}'.format(E, eta)
                maxevents = self.nevents
                inputfile = self.preformat_rootfile.format(save_dir, 2, E, eta, "")
                annotation = self.preformat_annotation.format(3, E, eta, self.nevents)
                outfile = self.preformat_rootfile.format(save_dir, 3, E, eta, "")
                outfile_miniaodsim = self.preformat_rootfile.format(save_dir, E, eta, 3, '_inMINIAODSIM')
                outfile_dqm = self.preformat_rootfile.format(save_dir, E, eta, 3, '_inDQM')
                template_text = template_text.format(maxevents, inputfile, annotation, outfile, outfile_miniaodsim, outfile_dqm)
                
                # Name of cfg file
                cfg_filename = self.preformat_cfgfile.format(dir_step3, 3, E, eta)

                # Set directory to save root files
                if not os.path.exists(save_dir):
                    os.makedirs(save_dir)
                
                with open(cfg_filename, 'w') as cfg_file:
                    cfg_file.write(template_text)
                    
        return

