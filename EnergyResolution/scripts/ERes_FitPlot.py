#!/usr/bin/env python2

# For reference, see: https://root.cern.ch/root/htmldoc/guides/users-guide/FittingHistograms.html

from ROOT import *
from math import sqrt
from array import array

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
    
    
def defFitEResolution (tgrapherrors, mode='num'):
    # type: TGraphErrors -> dict('key': 'float') (mode='num')
    # type: TGraphErrors -> TGraphErrors (mode='obj')
    # type: TGraphErrors -> (TGraphErrors, dict) (mode='all')
    " Fit dE/E to TGraphErrors. "

def getTGraphErrors (list_filenames, treename, branchname, histname):
    # type: (str, str, str, str) -> TGraphErrors
    " Plot resolution for different energy values. "
    
    mean_array, mean_error_array = array('f'), array('f')
    ERes_array, ERes_error_array = array('f'), array('f')
    
    for filename in list_filenames:
    
        infile = TFile.Open(filename, 'READ')
        tree = infile.Get(treename)
        
        if branchname == '': # possible? Get(branchname + histname)
            hist = tree.Get(histname)
        else:
            hist = tree.Get(branchname).Get(histname)
            
        stats = fitGaussian(hist, mode='num')
        sigma, sigma_error = stats['sigma'], stats['sigma_error']
        mean, mean_error = stats['mean'], stats['mean_error']
        
        mean_array.append(mean)
        mean_error_array.append(mean_error)
        ERes_array.append(sigma/mean)
        ERes_error_array.append( (sigma/mean)*addQuad(sigma_error/sigma, mean_error/mean) )
    
    ERes_graph = TGraphErrors(len(mean_array), mean_array, ERes_array, mean_error_array, ERes_error_array)
    ERes_graph.SetLineColor(14)
    ERes_graph.SetLineWidth(1)
    ERes_graph.SetLineStyle(1)
    ERes_graph.SetMarkerColor(9)
    ERes_graph.SetMarkerSize(1)
    ERes_graph.SetMarkerStyle(21)
    ERes_graph.GetXaxis().SetTitle('E [GeV]')
    ERes_graph.GetXaxis().SetLimits(0, 150)
    ERes_graph.GetYaxis().SetTitle('#Delta{E}/E')
    ERes_graph.GetYaxis().SetLimits(0, 1)
    ERes_graph.SetTitle('Energy Resolution: HGCNose, single #gamma, |#eta| = 3.5')
    #ERes_graph.Draw('ACP')
    
    return ERes_graph
        
    
def fitGaussian (hist, mode='num'):
    # type: TH1 -> TCanvas (mode='draw')
    # type: TH1 -> dict('key': 'float') (mode='data') <- useful for Pandas DataFrame conversion
    # type: TH1 -> (TH1, TF1) (mode='obj')
    # type: TH1 -> (TCanvas, TH1, TF1, TFitResultPtr) (mode='all')
    " Fit un-normalized Gaussian to histogram. "
    
    Xaxis = hist.GetXaxis()    
    Yaxis = hist.GetYaxis()
    fit = TF1('fit', 'gaus', Xaxis.GetXmin(), Xaxis.GetXmax())
    fitresult = hist.Fit('fit', 'QS')
    
    if mode=='draw':
        c = TCanvas('Gaussian Fit')
        hist.Draw()
        fit.Draw('SAME')
        return c
    
    elif mode=='num':
        n = dict()
        n['scale'], n['scale_error'] = fitresult.Parameter(0), fitresult.ParError(0)
        n['mean'], n['mean_error'] = fitresult.Parameter(1), fitresult.ParError(1)
        n['sigma'], n['sigma_error'] = fitresult.Parameter(2), fitresult.ParError(2)
        return n
    
    elif mode=='obj':
        fit.GetHistogram().GetXaxis().SetTitle(Xaxis.GetTitle())
        fit.GetHistogram().GetYaxis().SetTitle(Yaxis.GetTitle())
        return hist, fit
    
    elif mode=='all':
        fit.GetHistogram().GetXaxis().SetTitle(Xaxis.GetTitle())
        fit.GetHistogram().GetYaxis().SetTitle(Yaxis.GetTitle())
        
        n = dict()
        n['scale'], n['scale_error'] = fitresult.Parameter(0), fitresult.ParError(0)
        n['mean'], n['mean_error'] = fitresult.Parameter(1), fitresult.ParError(1)
        n['sigma'], n['sigma_error'] = fitresult.Parameter(2), fitresult.ParError(2)
        
        c = TCanvas('Gaussian Fit')
        hist.Draw()
        fit.Draw('SAME')
        
        return c, hist, fit, n
        
    else:
        raise Exception(" Correct mode not provided in function \'fitGaussian\'. Valid modes are \'draw\', \'num\', \'obj\', and \'all\'.")
        
# ---- Main -----

def main():
    list_filenames = ['output/ERes_pt{}.root'.format(i) for i in range(1,10)]
    treename = 'analysis'
    histname = 'EDist'
    graph = getTGraphErrors(list_filenames, treename, '', histname)
    canvas = TCanvas('ERes_TGraphErrors', 'Energy Resolution Trend')
    graph.Draw('ACP')
    canvas.SaveAs('ERes_Trend.png')


if __name__=='__main__':
    main()
