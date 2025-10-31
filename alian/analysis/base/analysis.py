import heppyy

class BaseAnalysis(heppyy.GenericObject):
    _defaults = {}

    def __init__(self, **kwargs):
        super(BaseAnalysis, self).__init__(**kwargs)
        for k, val in self.__class__._defaults.items():
            if not hasattr(self, k) or getattr(self, k) is None:
                setattr(self, k, val)
