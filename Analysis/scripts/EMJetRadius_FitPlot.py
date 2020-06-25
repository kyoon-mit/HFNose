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


        
def fitGaussianProjectionY (TH2hist, mode='num'):
    # type: TH2 -> [(TCanvas, float)] (mode='draw')
    # type: TH2 -> dict('key': float) (mode='data')
    # type: TH2 -> [(TH1, TF1, float)] (mode='raw')
    # type: TH2 -> [(TH1, TF1, float)], dict('key': float) (mode='all')
    " Fits Gaussian for each Y projection along the X bins of a TH2. Return type always paired with float that denotes the X bin content. "
    
    nbins = TH2hist.GetXaxis().GetNbins()
    
    if mode=='draw':
        canvases = []
        for i in range(nbins):
            TH1hist = TH2hist.ProjectionY("proj{}".format(i+1), i+1, i+2, 'e')
            x = TH2hist.GetXaxis().FindBin(i+1)
            canvases.append((fitGaussian(TH1hist, 'draw'), x))
        return canvases
    
    elif mode=='num' or mode=='raw' or mode=='all':
        numbers = dict()
        numbers['radii'] = []
        numbers['means'], numbers['mean_errors'] = [], []
        numbers['sigmas'], numbers['sigma_errors'] = [], []
        
        raws = []
        
        for i in range(nbins):
            TH1hist = TH2hist.ProjectionY("proj{}".format(i+1), i+1, i+2, 'e')
            x = TH2hist.GetXaxis().GetBinLowEdge(i+1)
            numbers['radii'].append(x)
            print (x)    
            
            n = fitGaussian(TH1hist, 'num')
            numbers['means'].append((n['mean']))
            print (n['mean'])
            numbers['mean_errors'].append((n['mean_error']))
            numbers['sigmas'].append((n['sigma']))
            numbers['sigma_errors'].append((n['sigma_error']))
            
            th1, tf1 = fitGaussian(TH1hist, 'raw')
            raws.append((th1, tf1, x))
            
        if mode=='num':
            return numbers
        elif mode=='raw':
            return raws
        elif mode=='all':
            return raws, numbers

    else:
        raise Exception(" Correct mode not provided in function \'fitGaussian\'. Valid modes are \'draw\', \'num\', \'raw\', and \'all\'.")
        
        
def fitLinearJetContainment (TH2hist):
    # type: TH2 -> TFitResult, TF1,  TGraphErrors
    " Fit linear function to get relationship of jet energy containment and radius. frac_E / E = [0] + [1] * dR"
    
    numbers = fitGaussianProjectionY(TH2hist, 'num')
    
    radii = numbers['radii']
    means = numbers['means']
    uncertainties = numbers['sigma_errors']
    npoints = len(radii)
    
    graph = TGraphErrors(npoints, np.array(radii), np.array(means), np.zeros(npoints), np.array(uncertainties))
    fit = TF1('linearfit', "[0] + [1] * x", 0, max(radii) + 0.05)
    fitresult = graph.Fit('linearfit', 'QS')
    
    return fitresult, fit, graph
    
    
def getNinetyPercentRadiusperLayer (*TH2hist):
    # type: TH2, [TH2, ] -> (float)
    " Get 90% containment radius per layer. "
    

        
# ---- Main -----

def main():
    # File names
    tuple_pt = ('5')
    tuple_filenames = tuple('/home/kyoon/CMSSW_11_1_0_pre7_RECHIT/src/HGCNose/Analysis/output/EMShowers_pt%s.root' % (i) for i in tuple_pt)
    
    # Top save directory
    TOP_DIR = '/home/kyoon/CMSSW_11_1_0_pre7_RECHIT/src/HGCNose/Analysis/plots/MoliereRadius_Single_Photon'
        
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

if __name__=='__main__':
    main()
