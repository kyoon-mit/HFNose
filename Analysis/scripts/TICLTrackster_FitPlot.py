#!/usr/bin/env python2

# For reference, see: https://root.cern.ch/root/htmldoc/guides/users-guide/FittingHistograms.html

from ROOT import *
from math import sqrt, sinh
import numpy as np
import os

# PyConfig.IgnoreCommandLineOptions = True # Don't pass argv into TApplication
gROOT.SetBatch(True) # Enable batch mode (i.e. suppress graphic windows)

# ---- Function definitions ----

def addQuad (*args):
    # type: float, [float, ] -> float
    " Addition in quadrature. "
    ans_squared = 0
    for arg in args:
        ans_squared += arg**2
    return sqrt(ans_squared)
    

def addQuadExpression (*args):
    # type: str, [str, ] -> str
    " Math expression (string) of addition in quadrature. "
    ans_squared = ""
    nterms = len(args)
    for i in range(nterms):
        if i < nterms - 1:
            ans_squared += "(%s)**2 + " % (args[i])
        else:
            ans_squared += "(%s)**2" % (args[i])
    return "sqrt(%s)" % (ans_squared)
    
    
def getMathEResolution (graph):
    # type: TGraphErrors -> TFitResult
    " Fit dE/E to TGraphErrors. "
    
    Xaxis = graph.GetXaxis()
    Yaxis = graph.GetYaxis()
    
    function = addQuadExpression("[0]/sqrt(x)", "[1]")
    fit = TF1('EResFit', function, 0, Xaxis.GetXmax())
    fit.SetParameters(7,7)
    fitresult = graph.Fit('EResFit', 'QS')
    
    n = dict()
    n['stochastic_term'] = fitresult.Parameter(0)
    n['constant_term'] = fitresult.Parameter(1)
    
    return n


def fitGaussian (hist, mode='num'):
    # type: TH1 -> TCanvas (mode='draw')
    # type: TH1 -> dict('key': 'float') (mode='data') <- useful for Pandas DataFrame conversion
    # type: TH1 -> (TH1, TF1) (mode='raw')
    # type: TH1 -> (TH1, TF1, dict('key': 'float')) (mode='all')
    " Fit un-normalized Gaussian to histogram. "
    
    Xaxis = hist.GetXaxis()    
    Yaxis = hist.GetYaxis()
    fit = TF1('fit', 'gaus', Xaxis.GetXmin(), Xaxis.GetXmax())
    fitresult = hist.Fit('fit', 'QS')
    
    if mode=='draw':
        c = TCanvas(hist.GetTitle())
        hist.Draw()
        fit.Draw('SAME')
        return c
    
    elif mode=='num':
        n = dict()
        n['scale'], n['scale_error'] = fitresult.Parameter(0), fitresult.ParError(0)
        n['mean'], n['mean_error'] = fitresult.Parameter(1), fitresult.ParError(1)
        n['sigma'], n['sigma_error'] = fitresult.Parameter(2), fitresult.ParError(2)
        return n
    
    elif mode=='raw':
        fit.GetHistogram().GetXaxis().SetTitle(Xaxis.GetTitle())
        fit.GetHistogram().GetYaxis().SetTitle(Yaxis.GetTitle())
        return hist, fit
    
    elif mode=='all':        
        n = dict()
        n['scale'], n['scale_error'] = fitresult.Parameter(0), fitresult.ParError(0)
        n['mean'], n['mean_error'] = fitresult.Parameter(1), fitresult.ParError(1)
        n['sigma'], n['sigma_error'] = fitresult.Parameter(2), fitresult.ParError(2)
        
        return fit, n
        
    else:
        raise Exception(" Correct mode not provided in function \'fitGaussian\'. Valid modes are \'draw\', \'num\', \'raw\', and \'all\'.")


        
def setTGraphStyle (graph, histtitle, preset):
    # type: (TGraphErrors, str, int) -> None
    " Set style, title, labels, etc. Should add more styles here for readability. "
    
    if preset=='ERes_1':
        graph.SetTitle(histtitle) 
        graph.GetXaxis().SetTitle('E_{caloParticle} [GeV]')
        graph.GetXaxis().SetLimits(0, 550)
        graph.GetYaxis().SetTitle('{#sigma(E_{trackster})}/{E_{trackster} (%)')
        graph.SetMinimum(0)
        graph.SetLineColor(14)
        graph.SetMarkerColor(9)
        graph.SetLineWidth(1)
        graph.SetLineStyle(1)
        graph.SetMarkerSize(1)
        graph.SetMarkerStyle(21)
        
    elif preset=='ERes_2':
        graph.SetTitle(histtitle)
        graph.GetXaxis().SetTitle('E [GeV]')
        graph.GetXaxis().SetLimits(0, 150)
        graph.GetYaxis().SetTitle('#sigma_{E}/E (%)')
        graph.SetMinimum(0)
        graph.SetLineColor(14)
        graph.SetMarkerColor(46)
        graph.SetLineWidth(1)
        graph.SetLineStyle(1)
        graph.SetMarkerSize(1)
        graph.SetMarkerStyle(20)
        
    elif preset=='EMean_1':
        graph.SetTitle(histtitle)
        graph.GetXaxis().SetTitle('E [GeV]')
        graph.GetXaxis().SetLimits(0, 550)
        graph.GetYaxis().SetTitle('E [GeV]')
        graph.SetMinimum(0)
        graph.SetMaximum(150)
        graph.SetLineColor(14)
        graph.SetMarkerColor(9)
        graph.SetLineWidth(1)
        graph.SetLineStyle(1)
        graph.SetMarkerSize(1)
        graph.SetMarkerStyle(21)
        
    elif preset=='EMean_2':
        graph.SetTitle(histtitle)
        graph.GetXaxis().SetTitle('E [GeV]')
        graph.GetXaxis().SetLimits(0, 600)
        graph.GetYaxis().SetTitle('E [GeV]')
        graph.SetMinimum(0)
        graph.SetMaximum(150)
        graph.SetLineColor(14)
        graph.SetMarkerColor(46)
        graph.SetLineWidth(1)
        graph.SetLineStyle(1)
        graph.SetMarkerSize(1)
        graph.SetMarkerStyle(20)
        
    elif preset=='EMean_3':
        graph.SetTitle(histtitle)
        graph.GetXaxis().SetTitle('E [GeV]')
        graph.GetXaxis().SetLimits(0, 600)
        graph.GetYaxis().SetTitle('E [GeV]')
        graph.SetMinimum(0)
        graph.SetMaximum(150)
        graph.SetLineColor(14)
        graph.SetMarkerColor(14)
        graph.SetLineWidth(1)
        graph.SetLineStyle(1)
        graph.SetMarkerSize(2)
        graph.SetMarkerStyle(33)
        
    elif preset=='EScale_1':
        graph.SetTitle(histtitle)
        graph.GetXaxis().SetTitle('E_{caloParticle} [GeV]')
        graph.GetXaxis().SetLimits(0, 550)
        graph.GetYaxis().SetTitle('E_{caloParticle} - E_{trackster} [GeV]')
        graph.SetMinimum(0)
        graph.SetMaximum(200)
        graph.SetLineColor(14)
        graph.SetMarkerColor(9)
        graph.SetLineWidth(1)
        graph.SetLineStyle(1)
        graph.SetMarkerSize(1)
        graph.SetMarkerStyle(21)
        
    elif preset=='EtaScale_1':
        graph.SetTitle(histtitle)
        graph.GetXaxis().SetTitle('E_{caloParticle} [GeV]')
        graph.GetXaxis().SetLimits(0, 550)
        graph.GetYaxis().SetTitle('|E_{caloParticle} - E_{trackster}| [GeV]')
        graph.SetMinimum(0)
        graph.SetMaximum(0.05)
        graph.SetLineColor(14)
        graph.SetMarkerColor(9)
        graph.SetLineWidth(1)
        graph.SetLineStyle(1)
        graph.SetMarkerSize(1)
        graph.SetMarkerStyle(21)
    
    return

    
def getPlotsEResolutionComprehensive (tuple_filenames, tuple_E):
    # type: ((str), (str)) -> dict('key': 'TGraph' / 'list(hist)')
    " Plot resolution for different energy values. "
    
    num_points = len(tuple_filenames)

    # Containers
    truth_E_array = np.zeros(num_points)
    trackster_ERes_array, trackster_ERes_error_array = np.zeros(num_points), np.zeros(num_points)
    trackster_EScale_array, trackster_EScale_error_array = np.zeros(num_points), np.zeros(num_points)
    trackster_EtaScale_array = np.zeros(num_points)
    
    for i in range(num_points):
    
        filename = tuple_filenames[i]
        E = tuple_E[i]
        infile = TFile.Open(filename, 'READ')
        tree = infile.Get('Analysis_SingleElectron')
        
        ### Get summary histograms
        truthE = tree.Get('truthE')
        truthEta = tree.Get('truthEta')
        tracksterRawEDist = tree.Get('tracksterRawEDist')
        tracksterAbsEtaDist = tree.Get('tracksterAbsEtaDist')
        tracksterRawEScale = tree.Get('tracksterRawEScale')
        tracksterAbsEtaScale = tree.Get('tracksterAbsEtaScale')
        
        # Extend scope of histograms
        truthE.SetDirectory(0)
        truthEta.SetDirectory(0)
        tracksterRawEDist.SetDirectory(0)
        tracksterAbsEtaDist.SetDirectory(0)
        tracksterRawEScale.SetDirectory(0)
        tracksterAbsEtaScale.SetDirectory(0)
        
        ### Get stats
        
        # Truth Energy
        truth_E_array[i] = truthE.GetMean(1)
        
        # Trackster "Raw Energy" Resolution
        trackster_raw_energy_stats = fitGaussian(tracksterRawEDist, mode='num')
        trackster_mean_raw_energy, trackster_mean_error_raw_energy = trackster_raw_energy_stats['mean'], trackster_raw_energy_stats['mean_error']
        trackster_sigma_raw_energy, trackster_sigma_error_raw_energy = trackster_raw_energy_stats['sigma'], trackster_raw_energy_stats['sigma_error']
        
        trackster_raw_energy_scale_stats = fitGaussian(tracksterRawEScale, mode='num')
        trackster_mean_raw_energy_scale, trackster_sigma_raw_energy_scale = trackster_raw_energy_scale_stats['mean'], trackster_raw_energy_scale_stats['sigma']
        
        trackster_ERes_array[i] = trackster_sigma_raw_energy/trackster_mean_raw_energy
        trackster_ERes_error_array[i] = trackster_ERes_array[i] * addQuad(trackster_sigma_error_raw_energy/trackster_sigma_raw_energy, trackster_mean_error_raw_energy/trackster_mean_raw_energy)
        
        # Trackster "Raw Energy" Scale
        trackster_EScale_array[i] = trackster_mean_raw_energy_scale
        trackster_EScale_error_array[i] = trackster_sigma_raw_energy_scale
        
        # Trackster Eta Scale
        trackster_EtaScale_array[i] = tracksterAbsEtaScale.GetMean(1)
        
        infile.Close()
        
    # Set dictionary for containing the return plots
    dict_plots = dict()
    
    dict_plots['trackster_RawERes_graph'] = TGraphErrors(num_points, truth_E_array, trackster_ERes_array*100, np.zeros(num_points), trackster_ERes_error_array*100)
    
    dict_plots['trackster_RawEScale_graph'] = TGraphErrors(num_points, truth_E_array, trackster_EScale_array, np.zeros(num_points), trackster_EScale_error_array)
    
    dict_plots['trackster_AbsEtaScale_graph'] = TGraph(num_points, truth_E_array, trackster_EtaScale_array)
    
    return dict_plots
    

def plotEResolutionFit (hits_ERes, clusters_ERes, top_save_dir):
    # (TGraph, TGraph, string) -> None
    " Plot energy resolution + fit"
        
    if not os.path.exists(top_save_dir + '/ERes_plots'):
        os.makedirs(top_save_dir + '/ERes_plots')
        
    # Set style
    group = 1 # group will later disappear
    
    setTGraphStyle(hits_ERes, "Energy Resolution (hits): HGCNose, single #gamma, |#eta| = 3.5", preset='ERes_%d' % (group))
    setTGraphStyle(clusters_ERes, "Energy Resolution (clusters): HGCNose, single #gamma, |#eta| = 3.5", preset='ERes_%d' % (group))
    
    # Get FitResult
    clusters_fit = getMathEResolution (clusters_ERes)
        
    s, c = clusters_fit['stochastic_term'], clusters_fit['constant_term']   
    formula_display = TPaveText(.6, .76, .88, .88, 'NDC')
    formula_display.SetFillColor(kWhite)
    formula_display.SetTextFont(43)
    formula_display.SetTextSize(14) # in pixels
    formula_display.SetTextColor(kRed+1)
    formula_display.AddText("#frac{#sigma_{E}}{E} = #frac{%.1f%%}{#sqrt{E}} #oplus %.1f%%" % (s, c))
        
    c1 = TCanvas("ERes_hits_TGraphErrors")
    hits_ERes.Draw('AP')
    
    c2 = TCanvas("ERes_clusters_TGraphErrors")
    c2.cd()
    clusters_ERes.Draw('AP')
    #clusters_fit.Draw('SAME')
    formula_display.Draw('SAME')
    c2.Update()
    
    c1.SaveAs(top_save_dir + '/ERes_plots/EResolution_hits.png')
    c2.SaveAs(top_save_dir + '/ERes_plots/EResolution_clusters.png')
    
    return
    
    
def plotGraphs (dict_plots, top_save_dir):
    # (dict(TGraphs), string) -> None
    " Plot all graphs"
    
    if not os.path.exists(top_save_dir + '/SingleElectron_plots'):
        os.makedirs(top_save_dir + '/SingleElectron_plots')

    setTGraphStyle(dict_plots["trackster_RawERes_graph"], "Trackster \"raw energy\" resolution", preset='ERes_1')
    setTGraphStyle(dict_plots["trackster_RawEScale_graph"], "Trackster \"raw energy\" scale", preset='EScale_1')
    setTGraphStyle(dict_plots["trackster_AbsEtaScale_graph"], "Trackster eta deviation", preset='EtaScale_1')
    
    c1 = TCanvas("RawERes")
    dict_plots["trackster_RawERes_graph"].Draw('ALP')
    
    c2 = TCanvas("RawEScale")
    c2.cd()
    dict_plots["trackster_RawEScale_graph"].Draw('ALP')
    
    c3 = TCanvas("AbsEtaScale")
    c3.cd()
    dict_plots["trackster_AbsEtaScale_graph"].Draw('ALP')
    
    c1.SaveAs(top_save_dir + '/SingleElectron_plots/RawERes.png')
    c2.SaveAs(top_save_dir + '/SingleElectron_plots/RawEScale.png')
    c3.SaveAs(top_save_dir + '/SingleElectron_plots/AbsEtaScale.png')
    
    return

 
# ---- Main -----

def main():
    # Set open and save diretories
    if not 'DIRANALYSIS_HGCNOSE' in os.environ:
        raise Exception ('DIRANALYSIS_HGCNOSE not set. Please first run config.sh which is in your top-level directory.')
    input_dir = os.environ['DIRANALYSIS_HGCNOSE'] + '/output'
    top_save_dir = os.environ['DIRANALYSIS_HGCNOSE'] + '/plots/SingleElectron'

    # File names
    tuple_E = ('50', '100', '150', '200', '250', '450', '500')
    tuple_filenames = tuple(input_dir + '/Single_Electron_E%s.root' % (i) for i in tuple_E)
        
    # Get the comprehensive dictionary of plots
    dict_plots = getPlotsEResolutionComprehensive(tuple_filenames, tuple_E)
    
    # Save Energy Resolution plots
    plotGraphs (dict_plots, top_save_dir)


if __name__=='__main__':
    main()
