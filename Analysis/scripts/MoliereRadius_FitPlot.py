#!/usr/bin/env python2

# For reference, see: https://root.cern.ch/root/htmldoc/guides/users-guide/FittingHistograms.html

from ROOT import *
import numpy as np
import os

# PyConfig.IgnoreCommandLineOptions = True # Don't pass argv into TApplication
gROOT.SetBatch(True) # Enable batch mode (i.e. suppress graphic windows)

# ---- Function definitions ----

def fitGaussian (hist, mode='num'):
    # type: TH1 -> TCanvas (mode='draw')
    # type: TH1 -> dict('key': float) (mode='data') <- useful for Pandas DataFrame conversion
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
    
        
def drawOverlayPlots (list_th1, legends, stack_title, save_title, TOP_DIR, pavetext=""):
    # type: [TH1], [str], str -> TCanvas
    " Draw multiple th1 histograms on the same canvas. "
    
    if not os.path.exists(TOP_DIR + '/pt_plots'):
        os.makedirs(TOP_DIR + '/pt_plots')

    hs = THStack("hs", stack_title)
    
    if not len(list_th1) == len(legends):
        raise Exception("Length of lists are different in input parameters.")
    
    gStyle.SetPalette(kRainBow)
    
    legend = TLegend(0.75, 0.50, 0.87, 0.77)
    legend.SetTextFont(43)
    legend.SetTextSize(13)
    legend.SetTextColor(kGray+2)
    
    for i in range(len(list_th1)):
        hist = list_th1[i]
        fit = hist.GetFunction('fit')
        hist.GetListOfFunctions().Remove(fit)
        hs.Add(hist)
        # hs.SetDirectory(0)
        
        legend.AddEntry(hist, legends[i], 'fl')
        
    c = TCanvas()
    hs.Draw("PLC NOSTACK")
    hs.GetXaxis().SetTitle("R (mm)")
    hs.GetYaxis().SetTitle("count")
    
    # PaveText
    display = TPaveText(.58, .80, .88, .88, 'NDC')
    display.SetFillColor(kWhite)
    display.SetTextFont(43)
    display.SetTextSize(14) # in pixels
    display.SetTextColor(kBlue)
    display.AddText(pavetext)
    display.Draw('SAME')
    
    # Draw legend
    legend.Draw('SAME')
    
    c.Update()
    c.SaveAs(TOP_DIR + '/pt_plots/' + save_title)
        
        
def getPlotsMoliereRadiusComprehensive (tuple_filenames, tuple_pt, Nlayers, TOP_DIR):
    # type: ((str), (str), int) -> dict('key': <TGraph> or <[THist]>)
    " Fit linear function to get relationship of jet energy containment and radius. frac_E / E = [0] + [1] * dR"
    
    eta = 3.5
    
    num_files = len(tuple_filenames)
    
    for i in range(num_files): # per file
    
        filename = tuple_filenames[i]
        pt = tuple_pt[i]
        infile = TFile.Open(filename, 'READ')
        tree = infile.Get('Analysis_EMShowers')
        
        ### Get detailed histograms
        # Container for 90% radii PER LAYER
        ninety_percent_draw_gauss = [] # list of histograms
        ninety_percent_mean_array, ninety_percent_mean_error_array = np.zeros(Nlayers), np.zeros(Nlayers)
        ninety_percent_sigma_array, ninety_percent_sigma_error_array = np.zeros(Nlayers), np.zeros(Nlayers)
        
        for layer in range(1, Nlayers+1):
            # Fit Gaussian to layer histogram
            ninety_percent_R_layer = tree.Get("Ninety_Percent_R_layer%d" % (layer))
            gauss_stats = fitGaussian(ninety_percent_R_layer)
            
            # Extend scope of histograms
            ninety_percent_R_layer.SetDirectory(0)
            
            # Fill container
            ninety_percent_draw_gauss.append(ninety_percent_R_layer)
            ninety_percent_mean_array[layer-1] = gauss_stats['mean']
            ninety_percent_sigma_array[layer-1] = gauss_stats['sigma']
        
        # Draw layer plots overlayed
        legends = ["layer%d" % (i) for i in range(1, Nlayers+1)]
        drawOverlayPlots(ninety_percent_draw_gauss, legends, "Ninety Percent R (layers stacked)", "pt%s_ninety_percent_R_overlay.png" % (pt), TOP_DIR, "Single #gamma, pt=%s GeV, E=%.2f GeV" % (pt, float(pt)*np.sinh(eta)))
        
        ### Get summary histograms
        ninety_percent_layer = tree.Get("Ninety_Percent_Layer")
        nevents_withjets_perlayer = tree.Get("Nevents_withJets_PerLayer")
        
        # Extend scope of histograms
        ninety_percent_layer.SetDirectory(0)
        nevents_withjets_perlayer.SetDirectory(0)
        
        # Draw other histograms
        c1 = TCanvas()
        c2 = TCanvas()
        
        c1.cd()
        ninety_percent_layer.Draw()
        c1.SaveAs(TOP_DIR + '/pt_plots/' + 'pt%s_ninety_percent_layer.png' % (pt))
        
        c2.cd()
        nevents_withjets_perlayer.Draw()
        c2.SaveAs(TOP_DIR + '/pt_plots/' + 'pt%s_nevents_withjets_perlayer.png' % (pt))

        
# ---- Main -----

def main():
    # File names
    tuple_pt = ('.33', '.66', '1', '2', '3', '4', '5', '6', '7', '8', '9')
    tuple_filenames = tuple('/home/kyoon/CMSSW_11_1_0_pre7_RECHIT/src/HGCNose/Analysis/output/EMShowers_pt%s.root' % (i) for i in tuple_pt)
    
    # Top save directory
    TOP_DIR = '/home/kyoon/CMSSW_11_1_0_pre7_RECHIT/src/HGCNose/Analysis/plots/MoliereRadius_Single_Photon'
    
    getPlotsMoliereRadiusComprehensive (tuple_filenames, tuple_pt, 8, TOP_DIR)
    
    """
    # Save TCanvases that show linear-fitted jet containment radii
    for i in range(len(tuple_pt)):
        
        pt = tuple_pt[i]
        save_dir = TOP_DIR + '/pt_plots/pt_%s' % (pt)
        
        if not os.path.exists(save_dir):
            os.makedirs(save_dir)
        
        filename = tuple_filenames[i]
        infile = TFile.Open(filename, 'READ')
        tree = infile.Get('Analysis_EMShowers')
        
        for layer in range(1, 9):
            TH2hist_R_Efrac = tree.Get('R_frac_containment_layer%s' % (layer))
            TH2hist_R_E = tree.Get('R_containment_layer%s' % (layer))
        
            c1 = TCanvas()
            c2 = TCanvas()
            
            c1.cd()
            fitresult1, fit1, graph1 = fitLinearJetContainment(TH2hist_R_Efrac)
            graph1.Draw()
            fit1.Draw('SAME')
            c1.SaveAs(save_dir + '/pt%s_R_frac_containment_layer%s.png' % (pt, layer))
            
            c2.cd()
            fitresult2, fit2, graph2 = fitLinearJetContainment(TH2hist_R_E)
            graph2.Draw()
            fit2.Draw('SAME')
            c2.SaveAs(save_dir + '/pt%s_R_containment_layer%s.png' % (pt, layer))
        """

if __name__=='__main__':
    main()
