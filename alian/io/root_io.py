
import ROOT

from heppyy.util.logger import Logger
log = Logger()

# singleton class for for analysis output root file
class SingleRootFile(object):

    __instance = None

    def __new__(cls, fname = 'output.root'):
        if SingleRootFile.__instance is None:
            SingleRootFile.__instance = object.__new__(cls)
            SingleRootFile.__instance.filename = fname
            SingleRootFile.__instance.root_file = ROOT.TFile(fname, 'RECREATE')            
            SingleRootFile.__instance.objects = []
            log.info(f'new root file: {SingleRootFile.__instance.root_file.GetName()}')
        return SingleRootFile.__instance

    def add(self, obj):
        if obj:
            self.objects.append(obj)

    def __del__(self):
        # Safely get the instance using getattr to avoid attribute errors during shutdown
        inst = getattr(SingleRootFile, '__instance', None)
        if inst is None:
            return
        
        # Safely get root_file using getattr
        root_file = getattr(inst, 'root_file', None)
        if root_file is None:
            return
        
        # Use try/except block to safely check if file is open and get name
        try:
            is_open_func = getattr(root_file, 'IsOpen', None)
            if callable(is_open_func) and is_open_func():
                get_name_func = getattr(root_file, 'GetName', None)
                file_name = get_name_func() if callable(get_name_func) else ''
                # Wrap logging and close in try/except for shutdown safety
                try:
                    log.info(f'closing {file_name}')
                except:
                    pass
                try:
                    inst.close()
                except:
                    pass
        except (AttributeError, TypeError):
            # Silently handle errors during interpreter shutdown
            pass

    def write(self):
        _rfile = self.root_file
        log.info(f'SingleRootFile:write root_file is {_rfile}')
        if _rfile:
            _rfile.cd()
            _ = [o.Write() for o in self.objects]
            _rfile.Write()
            log.info(f'SingleRootFile: wrote {_rfile.GetName()}')
                
    def close(self):
        _rfile = self.root_file
        if _rfile:
            log.info(f'SingleRootFile: purging {_rfile.GetName()}')
            self.write()
            _rfile = ROOT.TFile(self.filename, 'UPDATE')
            _rfile.Purge()
            _rfile.Write()
            _rfile.Close()
        else:
            log.debug('SingleRootFile: no root file to write')
        _rfile = None
        log.info(f'done {self.filename}')
