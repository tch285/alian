from enum import IntFlag, auto
from numpy import number

class InvalidSelectionError(Exception):
    """Exception raised for unrecognized selection."""
    pass

class Selection(IntFlag):
    @classmethod
    def max_value(cls):
        return 2 ** len(cls) - 1

    def __invert__(self):
        return self.__class__(self.__class__.max_value() - int(self))

    @classmethod
    def _missing_(cls, value):
        if isinstance(value, str):
            res = cls(0)
            for sel in value.split(','):
                res |= cls[sel.strip()]
        elif isinstance(value, number):
            res = cls(value.item())
        elif isinstance(value, int) and value >= 0:
            if value <= cls.max_value():
                return super()._missing_(value)
            else:
                raise InvalidSelectionError(f"{cls.__name__}: {value} is out of range ({cls.max_value()}).")
        elif isinstance(value, cls):
            return value
        else:
            raise InvalidSelectionError(f"{cls.__name__}: {value} is invalid. Valid selections: {cls._member_names_}")
        return res

# selections from O2Physics::daily-20240910-0200
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

# selections from O2Physics::daily-20240910-0200
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

# selections from O2Physics::daily-20240910-0200
class TrackSel(Selection):
    trackSign = auto() # 1 = positive, 0 = negative
    globalTrack = auto()
    qualityTrack = auto()
    qualityTrackWDCA = auto()
    hybridTrack = auto()