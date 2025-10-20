#!/usr/bin/env python3

# print('loading first')
# import heppyy
# fj = heppyy.load_cppyy('fastjet', True)
# print(fj)
# print(fj.contrib)
# print('loading first again')
# fj = heppyy.load_cppyy('fastjet', True)
# print(fj)
# print(fj.contrib)
# from alian.io import data_io
# from alian.utils import data_fj
# import numpy as np
# from alian.analysis.base.selections import EvSel, TrgSel
# from alian.analysis.base import BaseAnalysis
from alian.analysis.base.histogram import Histogram
import ROOT as r
from array import array as arr
import numpy as np
# r.TH1.AddDirectory(False)
# r.TH1.SetDefaultSumw2()
r.TH2.SetDefaultSumw2()
# h = r.TH1F('a', 'a', 100, 0, 100)
# h = r.TH1F(name = 'a', title = 'a', nbinsx= 2, xbins= arr('f', [0.2, 1.1, 2.3]))
# class th1sub(r.TH1F):
#     def __init__(self, *args, **kwargs):
#         super(type(self), self).__init__(*args, **kwargs)

# kw = {'name': 'a', 'title': 'a', 'nbinsx': 2, 'xbins': arr('f', [0.20000000298023224, 1.100000023841858, 2.299999952316284])}
# args = ()
# print(r.TH1F)
# args = ('a', 'a', 2, arr('f', [0.2, 1.1, 2.3]))
# h = r.TH1F(*args)
# h2 = r.TH1F(*args)
# h = th1sub(*args)
# h = th1sub('a', 'a', 2, arr('f', [0.2, 1.1, 2.3]))
# h = th1sub(name = 'a', title = 'a',nbinsx= 2,xbins= arr('f', [0.2, 1.1, 2.3]))
# h = r.TH1F(name = 'a', title = 'a',nbinsx= 2,xbins= arr('f', [0.2, 1.1, 2.3]))
# h = r.TH1F(name = 'a', title = 'a',nbinsx= 2,xbins= arr('f', [0.2, 1.1, 2.3]))
# h = r.TH1F(**kw)
# h = th1sub(**kw)
def linbins(xmin, xmax, nbins, dtype = None):
    return np.linspace(xmin, xmax, nbins+1, dtype = dtype)
def logbins(xmin, xmax, nbins, dtype = None):
    return arr('f', np.logspace(np.log10(xmin), np.log10(xmax), nbins+1))

h = Histogram(htype = '1f', name = 'a', title = 'b;x;y',
              nbinsx= 8,
              xbins= logbins(1e-2, 10, 8),
              )
print(type(h))
# h2 = Histogram(htype = 1, name = 'a222', title = 'b;x;y',
#               nbinsx= 2,
#               xbins= arr('f', [0.2, 1.1, 2.3]),
#               )
# h = Histogram(1, 'a', title = 'a',nbinsx= 2,xbins= arr('f', [0.2, 1.1, 2.3]), sumw2 = False)
# h2 = Histogram(1, 'a', title = 'a',nbinsx= 2,xbins= arr('f', [0.2, 1.1, 2.3]), sumw2 = False)
# h = Histogram(htype = 1, **kw, sumw2 = False)
# h = Histogram(1, *args)
# h = Histogram(1, 'a', 'a', 2,
#               arr('f', [0.2, 1.1, 2.3]))
# print(h.GetSumw2())
h.Fill(0.3, 4)
h.Fill(0.3, 3)
h.Fill(1.4, 9)
h.Fill(1.4, 12)
# print(h.is_sumw2)
# print(bool(h.GetSumw2()))
# with r.TFile('test.root', 'RECREATE') as f:
    # f.mkdir("mydir")
    # h.save(f)
    # h.save(f, "RECREATE")
    # h.save(f, 'NEW22', mode = "RECREATE", path = "sub1/sub2")
    # h2.save(f, 'dddd', mode = "RECREATE", path = "sub1/sub4")
# print(len(h.GetSumw2()))
# r.TH1.SetDefaultSumw2()
# r.TH2.SetDefaultSumw2()
# h.Sumw2(False)
# print(h.GetSumw2())
# print(bool(h.GetSumw2()))
# print(len(h.GetSumw2()))