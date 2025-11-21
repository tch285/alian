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

Also defines a HistogramCollection class that handles a group of Histograms.
"""

# r.TH1.AddDirectory(False)
r.TH1.SetDefaultSumw2()
r.TH2.SetDefaultSumw2()

def linbins(xmin, xmax, nbins):
    return np.linspace(xmin, xmax, nbins+1)

def logbins(xmin, xmax, nbins):
    return np.logspace(np.log10(xmin), np.log10(xmax), nbins+1)

class HistogramCollection:
    def __init__(self, cfg = None, sumw2  = None):
        self._bins = {}
        self._hists = {}
        self._sumw2 = sumw2
        if cfg is not None:
            self.load(cfg)

    def load(self, file):
        match file:
            case TextIOBase() if file.name.endswith(".yaml"):
                config = yaml.safe_load(file)
                self._from_dict(config)
            case str() if file.endswith(".yaml"):
                with open(file) as f:
                    config = yaml.safe_load(f)
                self._from_dict(config)
            case _:
                raise NotImplementedError(f'Cannot extract histograms from type: {type(file)}.')

    def _from_dict(self, config):
        self._setup_bins(config['bins'])
        self._add_hists(config['hists'])

    def __getitem__(self, name): return self._hists[name]
    def __setitem__(self, name, h): self._hists[name] = h
    @property
    def nbins(self):
        return {name: len(bins) - 1 for name, bins in self._bins.items()}
    @property
    def hists(self):
        def _yield_hists(hists):
            for key, hist in hists.items():
                if isinstance(hist, dict):
                    yield from _yield_hists(hist)
                else:
                    yield hist

        return _yield_hists(self._hists)


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
        paths = self._find_hist_paths(hists)
        for path in paths:
            *dirpath, name = path
            hist_cfg = hists
            for subdir in path:
                hist_cfg = hist_cfg[subdir]

            htype, title, *bin_names = hist_cfg
            self._validate_bin_names(bin_names)
            # interlace number of bins and bin arrays for each axis
            binnings = [val for name in bin_names for val in (self.nbins[name], self._bins[name])]
            hist_loc = self._hists
            # dirpath is empty if histogram saved to lowest directory
            for subdir in dirpath:
                # mimic nested dict behavior, since we don't want it after this
                try:
                    hist_loc = hist_loc[subdir]
                except KeyError:
                    hist_loc[subdir] = {}
                    hist_loc = hist_loc[subdir]
            hist_loc[name] = Histogram(htype, name, title, *binnings, path = "/".join(dirpath), sumw2 = self._sumw2)
        # print(self._hists)

    def _find_hist_paths(self, hists, current_key=[]):
        list_keys = []
        for key, value in hists.items():
            new_key = current_key + [key]
            if isinstance(value, dict):
                list_keys.extend(self._find_hist_paths(value, new_key))
            elif isinstance(value, list):
                list_keys.append(new_key)
        return list_keys

    def save(self, file, mode = "RECREATE"):
        match file:
            case r.TFile():
                print(f"Mode {mode} will be ignored")
                self._save_to_tfile(file)
            case str() if file.endswith(".root"):
                with r.TFile(file, mode) as f:
                    self._save_to_tfile(f)
            case _:
                raise TypeError(f'Unsupported file type {type(file)} for saving')
    def _save_to_tfile(self, f):
        for h in self.hists:
            h.save(f)
    def set_sumw2(self, sumw2):
        for h in self.hists:
            h.Sumw2(sumw2)

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
    # override self.path with new path if given
    if path is not None:
        self.path = path
    if isinstance(file, str):
        # if mode unspecified, set to UPDATE
        mode = kwargs.pop("mode", "UPDATE")
        print(f"Saving histogram to {file} ({mode})")
        with r.TFile(file, mode) as f:
            self._save_to_tfile(f, *args, **kwargs)
    elif isinstance(file, r.TFile):
        mode = kwargs.pop("mode", None)
        if mode is not None:
            print(f"Mode {mode} will be ignored.")
        self._save_to_tfile(file, *args, **kwargs)
    else:
        raise TypeError(f"Cannot save to {file} with type {type(file)}")
def _save_to_tfile(self, f, *args, **kwargs):
    # save to open ROOT file
    if self.path:
        tdir = f.Get(self.path)
        if not tdir:
            f.mkdir(self.path)
            tdir = f.Get(self.path)
        tdir.WriteObject(self, *args, **kwargs)
    else:
        f.WriteObject(self, *args, **kwargs)
def is_sumw2(self):
    return bool(self.GetSumw2())
def fill(self, *args, **kwargs):
    # either the weight is specified by keyword
    if 'w' in kwargs:
        fillw = kwargs.pop('w')
    # or all arguments were passed without keywords
    else:
        try:
            # weight is last arg
            fillw = args[self.GetDimension()]
        except IndexError:
            # or weight not given
            fillw = 1

    kwargs['w'] = fillw * self.event_weight
    super().Fill(*args, **kwargs)  # Calls ROOT's Fill


class HistogramMeta(type):
    _hsub_cache = {}

    def __call__(cls, htype, *args, **kwargs):
        sumw2 = kwargs.pop("sumw2", None)
        path = kwargs.pop("path", None)
        if htype not in cls._hsub_cache:
            hclass = cls.get_histogram_class(htype)

            def __init__(self, sumw2, path, *args, **kwargs):
                if sumw2 is not None:
                    # if sumw2 specified, then use it, otherwise keep ROOT default
                    self.Sumw2(sumw2)
                self.path = path
                self.event_weight = 1
            def __new__(subcls, sumw2, path, *args, **kwargs):
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
                    '_save_to_tfile': _save_to_tfile,
                    'is_sumw2': property(is_sumw2),
                    'fill': fill,
                }
            )

            cls._hsub_cache[htype] = hsubclass

        return cls._hsub_cache[htype](sumw2, path, *args, **kwargs)

class Histogram(metaclass = HistogramMeta):
    @staticmethod
    def get_histogram_class(htype):
        if htype not in HISTOGRAM_REGISTRY:
            valid = "\n".join([f"\t{cls.__name__}: [{', '.join(map(repr, valid_codes))}]" for cls, valid_codes in registry.items()])
            raise KeyError(f"Could not parse {htype} into a valid ROOT histogram type. Valid codes:\n{repr(valid)}"
            )
        return HISTOGRAM_REGISTRY.get(htype)

registry = {
    r.TH1F: [1, '1', '1F', '1f'],
    r.TH2F: [2, '2', '2F', '2f'],
    r.TH3F: [3, '3', '3F', '3f'],
    r.TH1D: ['1D', '1d'],
    r.TH2D: ['2D', '2d'],
    r.TH3D: ['3D', '3d'],
    r.TH1I: ['1I', '1i'],
    r.TH2I: ['2I', '2i'],
    r.TH3I: ['3I', '3i'],
}

# reverse registry mapping to map key to specific ROOT class
HISTOGRAM_REGISTRY = {
    key: value
    for value, keys in registry.items()
    for key in keys
}
