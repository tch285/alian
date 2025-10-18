from enum import Flag, auto
from numpy import number

# selections from O2Physics::daily-20240910-0200

class Selection(Flag):
    @classmethod
    def _missing_(cls, value):
        if isinstance(value, str):
            res = cls(0)
            for sel in value.split(','):
                res |= cls[sel.strip()]
        elif isinstance(value, int) and value >= 0:
            return super()._missing_(value)
        elif isinstance(value, number):
            res = cls(value.item())
        elif isinstance(value, cls):
            return value
        else:
            raise TypeError(f"Cannot parse {cls.__name__} '{value}' of type {type(value)}")
        return res

class TrgSel(Selection):
    noTrigSel = auto()
    JetChLowPt = auto()
    JetChHighPt = auto()
    TrackLowPt = auto()
    TrackHighPt = auto()
    JetD0ChLowPt = auto()
    JetD0ChHighPt = auto()
    JetLcChLowPt = auto()
    JetLcChHighPt = auto()
    EMCALReadout = auto()
    JetFullHighPt = auto()
    JetFullLowPt = auto()
    JetNeutralHighPt = auto()
    JetNeutralLowPt = auto()
    GammaVeryHighPtEMCAL = auto()
    GammaVeryHighPtDCAL = auto()
    GammaHighPtEMCAL = auto()
    GammaHighPtDCAL = auto()
    GammaLowPtEMCAL = auto()
    GammaLowPtDCAL = auto()
    GammaVeryLowPtEMCAL = auto()
    GammaVeryLowPtDCAL = auto()

class EvSel(Selection):
    sel8 = auto()
    sel8Full = auto()
    sel8FullPbPb = auto()
    selMC = auto()
    selMCFull = auto()
    selMCFullPbPb = auto()
    selUnanchoredMC = auto()
    selTVX = auto()
    sel7 = auto()
    sel7KINT7 = auto()

class TrackSel(Selection):
    trackSign = auto() # 1 = positive, 0 = negative
    globalTrack = auto()
    qualityTrack = auto()
    qualityTrackWDCA = auto()
    hybridTrack = auto()