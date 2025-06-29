from enum import Flag, auto

class TrgSel(Flag):
    noTrigSel = auto(),
    JetChLowPt = auto(),
    JetChHighPt = auto(),
    TrackLowPt = auto(),
    TrackHighPt = auto(),
    JetD0ChLowPt = auto(),
    JetD0ChHighPt = auto(),
    JetLcChLowPt = auto(),
    JetLcChHighPt = auto(),
    EMCALReadout = auto(),
    JetFullHighPt = auto(),
    JetFullLowPt = auto(),
    JetNeutralHighPt = auto(),
    JetNeutralLowPt = auto(),
    GammaVeryHighPtEMCAL = auto(),
    GammaVeryHighPtDCAL = auto(),
    GammaHighPtEMCAL = auto(),
    GammaHighPtDCAL = auto(),
    GammaLowPtEMCAL = auto(),
    GammaLowPtDCAL = auto(),
    GammaVeryLowPtEMCAL = auto(),
    GammaVeryLowPtDCAL = auto()

class EvSel(Flag):
    sel8 = auto(),
    sel8Full = auto(),
    sel7 = auto(),
    selMC = auto(),
    selUnanchoredMC = auto(),
    sel7KINT7 = auto()

class TrackSel(Flag):
    trackSign = auto(),
    globalTrack = auto(),
    qualityTrack = auto(),
    hybridTrack = auto(),
    uniformTrack = auto(),
    uniformTrackWoDCA = auto()