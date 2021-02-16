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
        
    elif preset=='EFraction_1':
        graph.SetTitle(histtitle)
        graph.GetXaxis().SetTitle('E_{caloParticle} [GeV]')
        graph.GetXaxis().SetLimits(0, 550)
        graph.GetYaxis().SetTitle('E_{trackster} / E_{caloParticle} [%]')
        graph.SetMinimum(0)
        graph.SetMaximum(100)
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
        graph.GetYaxis().SetTitle('|#eta_{caloParticle} - #eta_{trackster}| [GeV]')
        graph.SetMinimum(0)
        graph.SetMaximum(0.05)
        graph.SetLineColor(14)
        graph.SetMarkerColor(9)
        graph.SetLineWidth(1)
        graph.SetLineStyle(1)
        graph.SetMarkerSize(1)
        graph.SetMarkerStyle(21)
        
    elif preset=='RScale_1':
        graph.SetTitle(histtitle)
        graph.GetXaxis().SetTitle('E_{caloParticle} [GeV]')
        graph.GetXaxis().SetLimits(0, 550)
        graph.GetYaxis().SetTitle('|#Delta R_{caloParticle - trackster}| [GeV]')
        graph.SetMinimum(0)
        graph.SetMaximum(0.05)
        graph.SetLineColor(14)
        graph.SetMarkerColor(9)
        graph.SetLineWidth(1)
        graph.SetLineStyle(1)
        graph.SetMarkerSize(1)
        graph.SetMarkerStyle(21)
    
    return

    
def getPlotsComprehensive (filenames, energies, eta, trackster_name, top_save_dir):
    # type: ((str), (str)), str) -> dict('key': 'TGraph' / 'TH1')
    " Plot resolution for different energy values. "
    
    num_points = len(filenames)

    # Containers for graph points
    truth_E_array = np.zeros(num_points)
    trackster_ERes_array, trackster_ERes_error_array = np.zeros(num_points), np.zeros(num_points)
    trackster_EScale_array, trackster_EScale_error_array = np.zeros(num_points), np.zeros(num_points)
    trackster_EFraction_array, trackster_EFraction_error_array = np.zeros(num_points), np.zeros(num_points)
    trackster_EtaScale_array = np.zeros(num_points)
    trackster_RScale_array = np.zeros(num_points)
    
    # Set dictionary for the plots
    dict_plots = dict()
    
    for i in range(num_points):
    
        filename = filenames[i]
        E = energies[i]
        infile = TFile.Open(filename, 'READ')
        tree = infile.Get('TICLAnalyzer')
        
        ### Get summary histograms
        truthE = tree.Get('truthE')
        truthEta = tree.Get('truthEta')
        tracksterRawEDist = tree.Get('RawEDist_{}'.format(trackster_name))
        tracksterRawEScale = tree.Get('RawEScale_{}'.format(trackster_name))
        tracksterDeltaRDist = tree.Get('DeltaR_{}'.format(trackster_name))
        clusterSumE = tree.Get('EDist_layerClusters_scalar_sum')
        recHitsSumE = tree.Get('EDist_recHits')
        
        truthE.SetDirectory(0)
        truthEta.SetDirectory(0)
        tracksterRawEDist.SetDirectory(0)
        tracksterAbsEtaDist.SetDirectory(0)
        tracksterRawEScale.SetDirectory(0)
        tracksterDeltaRDist.SetDirectory(0)
        clusterSumE.SetDirectory(0)
        recHitsSumE.SetDirectory(0)
        
        ### Get stats
        
        # Truth Energy
        truth_E_array[i] = truthE.GetMean(1)
        
        # Trackster "Raw Energy" Resolution
        trackster_raw_energy_stats = fitGaussian(tracksterRawEDist, mode='num')
        trackster_mean_raw_energy, trackster_mean_error_raw_energy = trackster_raw_energy_stats['mean'], trackster_raw_energy_stats['mean_error']
        trackster_sigma_raw_energy, trackster_sigma_error_raw_energy = trackster_raw_energy_stats['sigma'], trackster_raw_energy_stats['sigma_error']
        
        trackster_ERes_array[i] = trackster_sigma_raw_energy/trackster_mean_raw_energy
        trackster_ERes_error_array[i] = trackster_ERes_array[i] * addQuad(trackster_sigma_error_raw_energy/trackster_sigma_raw_energy, trackster_mean_error_raw_energy/trackster_mean_raw_energy)
        
        # Trackster "Raw Energy" Scale
        trackster_raw_energy_scale_stats = fitGaussian(tracksterRawEScale, mode='num')
        trackster_mean_raw_energy_scale, trackster_sigma_raw_energy_scale = trackster_raw_energy_scale_stats['mean'], trackster_raw_energy_scale_stats['sigma']
        
        trackster_EScale_array[i] = trackster_mean_raw_energy_scale
        trackster_EScale_error_array[i] = trackster_sigma_raw_energy_scale
        
        # Trackster "Raw Energy" / Truth Energy
        trackster_EFraction_array[i] = trackster_mean_raw_energy / truthE.GetMean(1)
        trackster_EFraction_error_array[i] = trackster_mean_error_raw_energy / truthE.GetMean(1)
        
        # Trackster Delta R Distribution
        trackster_DeltaR_array[i] = tracksterDeltaRDist.GetMean(1)
        
        # RecHits + Clusters + Trackster comparison plots
        plotEDistHistograms(recHitsSumE, clusterSumE, tracksterRawEDist, E, eta, top_save_dir)
        
        infile.Close()
        
    
    dict_plots['trackster_RawERes_graph'] = TGraphErrors(num_points, truth_E_array, trackster_ERes_array*100, np.zeros(num_points), trackster_ERes_error_array*100)
    
    dict_plots['trackster_RawEScale_graph'] = TGraphErrors(num_points, truth_E_array, trackster_EScale_array, np.zeros(num_points), trackster_EScale_error_array)
    
    dict_plots['trackster_RawEFraction_graph'] = TGraphErrors(num_points, truth_E_array, trackster_EFraction_array*100, np.zeros(num_points), trackster_EFraction_error_array)
    
    dict_plots['trackster_DeltaR_graph'] = TGraph(num_points, truth_E_array, trackster_DeltaR_array)
    
    return dict_plots
    
    
def plotEDistHistograms (hits_EDist, clusters_EDist, trackster_EDist, E, eta, top_save_dir):
    # (TH1, TH1, TH1, string) -> None
    " Plot histograms"
    
    if not os.path.exists(top_save_dir + '/EDistribution_plots'):
        os.makedirs(top_save_dir + '/EDistribution_plots')
        
    c = TCanvas("EDistribution_E{}".format(E))
    hits_EDist.Draw()
    clusters_EDist.Draw('SAME')
    trackster_EDist.Draw('SAME')
    
    hits_EDist.SetLineColor(kBlue)
    clusters_EDist.SetLineColor(kGreen+1)
    trackster_EDist.SetLineColor(kRed+1)
    
    c.Update()
    
    c.SaveAs(top_save_dir + '/EDistribution_plots/comparison_E{}_eta{}.jpg'.format(E, eta))
    
    

def plotEResolutionFit (hits_EGraph, clusters_EGraph, trackster_EGraph, top_save_dir):
    # (TGraph, TGraph, TGraph, string) -> None
    " Plot energy resolution + fit"
        
    if not os.path.exists(top_save_dir + '/EResolution_plots'):
        os.makedirs(top_save_dir + '/EResolution_plots')
        
    # Set style
    group = 1
    
    setTGraphStyle(hits_EGraph, "Energy Resolution (hits): HGCNose, single #gamma, |#eta| = 3.5", preset='ERes_%d' % (group))
    setTGraphStyle(clusters_EGraph, "Energy Resolution (clusters): HGCNose, single #gamma, |#eta| = 3.5", preset='ERes_%d' % (group))
    setTGraphStyle(trackster_EGraph, "Energy Resolution (trackster): HGCNose, single #gamma, |#eta| = 3.5", preset='ERes_%d' % (group))
    
    # Get FitResult
    formula_display_list = []
    for EGraph in (hits_EGraph, clusters_EGraph, trackster_EGraph):
        fit = getMathEResolution (EResPlot)
        
        s, c = clusters_fit['stochastic_term'], clusters_fit['constant_term']   
        formula_display = TPaveText(.6, .76, .88, .88, 'NDC')
        formula_display.SetFillColor(kWhite)
        formula_display.SetTextFont(43)
        formula_display.SetTextSize(14) # in pixels
        formula_display.SetTextColor(kRed+1)
        formula_display.AddText("#frac{#sigma_{E}}{E} = #frac{%.1f%%}{#sqrt{E}} #oplus %.1f%%" % (s, c))
        formula_display_list.append(formula_display)
        
    c1 = TCanvas("ERes_hits_TGraphErrors")
    hits_ERes.Draw('AP')
    formula_display_list[0].Draw('SAME')
    
    c2 = TCanvas("ERes_clusters_TGraphErrors")
    c2.cd()
    clusters_ERes.Draw('AP')
    formula_display_list[1].Draw('SAME')
    c2.Update()
    
    c3 = TCanvas("ERes_trackster_TGraphErrors")
    c3.cd()
    clusters_ERes.Draw('AP')
    formula_display_list[2].Draw('SAME')
    c3.Update()
    
    c1.SaveAs(top_save_dir + '/EResolution_plots/EResolution_hits.jpg')
    c2.SaveAs(top_save_dir + '/EResolution_plots/EResolution_clusters.jpg')
    c3.SaveAs(top_save_dir + '/EResolution_plots/EResolution_trackster.jpg')
    
    return
    
    
def plotGraphs (dict_plots, top_save_dir):
    # (dict(TGraphs), string) -> None
    " Plot all graphs"
    
    if not os.path.exists(top_save_dir + '/TICL_plots'):
        os.makedirs(top_save_dir + '/TICL_plots')

    setTGraphStyle(dict_plots["trackster_RawERes_graph"], "Trackster raw energy resolution", preset='ERes_1')
    setTGraphStyle(dict_plots["trackster_RawEScale_graph"], "Trackster raw energy scale", preset='EScale_1')
    setTGraphStyle(dict_plots["trackster_RawEFraction_graph"], "Trackster (raw energy) / (truth energy)", preset='EFraction_1')
    setTGraphStyle(dict_plots["trackster_DeltaR_graph"], "#Delta |R_{caloParticle} - R_{trackster}|", preset='RScale_1')
    
    c1 = TCanvas("RawERes")
    dict_plots["trackster_RawERes_graph"].Draw('ALP')
    
    c2 = TCanvas("RawEScale")
    c2.cd()
    dict_plots["trackster_RawEScale_graph"].Draw('ALP')
    
    c3 = TCanvas("RawEFraction")
    c3.cd()
    dict_plots["trackster_RawEFraction_graph"].Draw('ALP')
    
    c4 = TCanvas("DeltaR")
    c4.cd()
    dict_plots["trackster_DeltaR_graph"].Draw('ALP')
    
    c1.SaveAs(top_save_dir + '/graph_plots/RawERes.png')
    c2.SaveAs(top_save_dir + '/graph_plots/RawEScale.png')
    c3.SaveAs(top_save_dir + '/graph_plots/RawEFraction.png')
    c4.SaveAs(top_save_dir + '/graph_plots/DeltaR.png')
    
    return

 
# ---- Main -----

def main():
    # Set open and save diretories
    if not 'DIRANALYSIS_HGCNOSE' in os.environ:
        raise Exception ('DIRANALYSIS_HGCNOSE not set. Please first run config.sh which is in your top-level directory.')
    input_dir = os.environ['DIRANALYSIS_HGCNOSE'] + '/output/pid_0.5'
    top_save_dir = os.environ['DIRANALYSIS_HGCNOSE'] + '/plots/HGCNose'

    # File names and energies
    energies = ('10', '50', '100', '200', '300', '400', '500')
    eta = 3.5
    filenames = tuple(input_dir + '/Single_Photon_E{}_eta{}.root'.format(E, eta) for E in energies)
        
    # Get the comprehensive dictionary of plots
    trackster_name = "tracksterHFNoseEM"
    dict_plots = getPlotsComprehensive(filenames, energies, eta, trackster_name, top_save_dir)
    
    # Save Energy Resolution plots
    plotGraphs (dict_plots, top_save_dir)


if __name__=='__main__':
    main()
