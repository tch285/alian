from functools import singledispatchmethod
from io import TextIOBase

import numpy as np
import ROOT as r
import yaml

"""
Defines a Histogram class, adding additional functionality on top of
ROOT histograms. Define a Histogram like
    h = Histogram("1f", "name", "title", 10, 0, 10)
which will create a TH1F with the given attributes. Histograms are
automatically detached from gDirectory, and a custom save function is
provided, i.e. h.save(filename). Pass a keyword argument sumw2 to the
constructor to specify whether to store proper errors; if unspecified,
the ROOT default is used. One can set the default with e.g.
    ROOT.TH1.SetDefaultSumw2(True)
within your script.

Also defines a HistogramCollection class
"""

r.TH1.AddDirectory(False)
# r.TH1.SetDefaultSumw2()
# r.TH2.SetDefaultSumw2()

def linbins(xmin, xmax, nbins):
    return np.linspace(xmin, xmax, nbins+1)

def logbins(xmin, xmax, nbins):
    return np.logspace(np.log10(xmin), np.log10(xmax), nbins+1)

class HistogramCollection:
    def __init__(self, cfg = None):
        self._bins = {}
        self._hists = {}
        if cfg is not None:
            self.load(cfg)

    @singledispatchmethod
    def load(cls, file):
        raise NotImplementedError(f'Cannot extract histograms from type: {type(file)}.')
    @load.register
    def _load(self, file: str):
        if file.endswith(".yaml"):
            with open(file) as f:
                self._from_yaml(f)
        else:
            raise NotImplementedError(f'Cannot extract histograms from file: {file}.')
    @load.register
    def _load(self, file: TextIOBase):
        if file.name.endswith(".yaml"):
            return self._from_yaml(file)
        else:
            raise NotImplementedError(f'Cannot extract histograms from file: {file.name}.')

    def _from_yaml(self, cfg):
        with open(cfg) as stream:
            config = yaml.safe_load(stream)
        self._setup_bins(config['bins'])
        self._add_hists(config['hists'])

    def __getitem__(self, name): return self._hists[name]
    def __setitem__(self, name, h): self._hists[name] = h
    @property
    def nbins(self):
        return {name: len(bins) - 1 for name, bins in self._bins.items()}

    def _setup_bins(self, bins):
        for name, arglist in bins.items():
            self._bins[name] = self._calc_bins(*arglist)
    def _calc_bins(self, binning, *edge_info):
        match binning:
            case 'logarithmic' | 'logarithm' | 'log' | 'lg':
                binmin, binmax, nbins = edge_info
                return logbins(binmin, binmax, nbins)
            case 'linear' | 'line' | 'lin' | 'ln':
                binmin, binmax, nbins = edge_info
                return linbins(binmin, binmax, nbins)
            case 'custom' | 'cst' | 'c':
                return np.array(edge_info[0], dtype = float)
            case _:
                raise ValueError(f"Binning type {binning} invalid.")

    def _validate_bin_names(self, names):
        """Check bin names from Histogram configs are defined."""
        invalid_names = [name for name in names if name not in self._bins]
        if invalid_names:
            raise KeyError(f"Requested undefined binning: {', '.join(invalid_names)}")
    def _add_hists(self, hists):
        for name, hist_cfg in hists.items():
            htype, title, *bin_names = hist_cfg
            self._validate_bin_names(bin_names)
            binnings = [val for name in bin_names for val in (self.nbins[name], self._bins[name])]
            self[name] = Histogram(htype, name, title, *binnings)

    def save(self, file, mode = "RECREATE"):
        match file:
            case r.TFile():
                [file.WriteObject(h) for h in self._hists.values()]
            case str() if file.endswith(".root"):
                with r.TFile(file, mode) as f:
                    [f.WriteObject(h) for h in self._hists.values()]
            case _:
                raise TypeError(f'Unsupported file type {type(file)} for saving')


def save(self, file, *args, **kwargs):
    """
    There are two uses cases:
        If trying to save one histogram to a file, pass
        file as a string, and make sure to specify the file
        opening mode with the `mode` kwarg (UPDATE by default).

        If trying to save multiple histograms, open a `with`
        context block and do `Histogram.save(file)`. If the
        mode is specified, a warning will be shown.
    """
    path = kwargs.pop("path", None)
    if isinstance(file, str):
        # if mode unspecified, set to UPDATE
        mode = kwargs.pop("mode", "UPDATE")
        print(f"Saving histogram to {file} ({mode})")
        if ":" in file:
            filename, subdirs = file.split(":")
        else:
            filename, subdirs = file, None
        with r.TFile(filename, mode) as f:
            if subdirs is not None:
                file.mkdir(path, "")
                # HACK: can't get subsubdirectories with just above prior to 6.36
                tdir = file.Get(path)
                tdir.WriteObject(self, *args, **kwargs)
            else:
                f.WriteObject(self, *args, **kwargs)
    elif isinstance(file, r.TFile):
        mode = kwargs.pop("mode", None)
        if mode is not None:
            print(f"Mode {mode} will be ignored.")
        if file.IsOpen():
            if isinstance(path, str):
                file.mkdir(path, "")
                # HACK: can't get subsubdirectories with just above prior to 6.36
                tdir = file.Get(path)
                tdir.WriteObject(self, *args, **kwargs)
            else:
                file.WriteObject(self, *args, **kwargs)
        else:
            raise RuntimeError(f"TFile {file} is not open?")
    else:
        raise TypeError(f"Cannot save to {file} with type {type(file)}")
def is_sumw2(self):
    return bool(self.GetSumw2())

class HistogramMeta(type):
    def __init__(cls, name, bases, attrs):
        super().__init__(name, bases, attrs)
        cls._hsub_cache = {}

    def __call__(cls, htype, *args, **kwargs):
        sumw2 = kwargs.pop("sumw2", None)
        if htype not in cls._hsub_cache:
            hclass = cls.get_histogram_class(htype)

            def __init__(self, sumw2, *args, **kwargs):
                if sumw2 is not None:
                    # if sumw2 specified, then use it, otherwise keep ROOT default
                    self.Sumw2(sumw2)
            def __new__(subcls, sumw2, *args, **kwargs):
                # workaround to use constructor kwargs (some cppyy BS)
                instance = hclass(*args, **kwargs)
                instance.__class__ = subcls
                return instance

            hsubclass = type(
                f"{hclass.__name__}_sub",
                (hclass,),
                {   '__init__': __init__,
                    '__new__': __new__,
                    'save': save,
                    'is_sumw2': property(is_sumw2),
                }
            )

            cls._hsub_cache[htype] = hsubclass

        return cls._hsub_cache[htype](sumw2, *args, **kwargs)

class Histogram(metaclass = HistogramMeta):
    @staticmethod
    def get_histogram_class(htype):
        if htype not in HISTOGRAM_REGISTRY:
            raise KeyError(f"Could not parse {htype} into a valid ROOT histogram type.")
        return HISTOGRAM_REGISTRY.get(htype)

registry = {
    r.TH1F: [1, '1F', '1f'],
    r.TH2F: [2, '2F', '2f'],
    r.TH3F: [3, '3F', '3f'],
    r.TH1D: ['1D', '1d'],
    r.TH2D: ['2D', '2d'],
    r.TH3D: ['3D', '3d'],
    r.TH1I: ['1I', '1i'],
    r.TH2I: ['2I', '2i'],
    r.TH3I: ['3I', '3i'],
}

HISTOGRAM_REGISTRY = {
    key: value
    for value, keys in registry.items()
    for key in keys
}
