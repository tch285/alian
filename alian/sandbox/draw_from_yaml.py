#!/usr/bin/env python3

import yaml
import ROOT
import numpy as np
import array
import sys
import os
import argparse
from jinja2 import Environment, FileSystemLoader

def logbins(xmin, xmax, nbins):
    xmin = max(xmin, 1e-10)
    xmax = max(xmax, 1e-10)
    lspace = np.logspace(np.log10(xmin), np.log10(xmax), nbins+1)
    arr = array.array('d', lspace)
    return arr

def logbins_from_config(hist_config):
    return logbins(float(hist_config['xrange'][0]), float(hist_config['xrange'][1]), int(hist_config['xnbins']))

def get_hist_clone_from_file(hist_config, hist_new_name, output_file):
    hist_config_split = hist_config.split(':')
    if len(hist_config_split) != 2:
        print('Invalid hist_config:', hist_config)
        sys.exit(1)
    file_name = hist_config_split[0]
    hist_name = hist_config_split[1]
    root_file = ROOT.TFile.Open(file_name)
    hist = root_file.Get(hist_name)
    output_file.cd()
    hist_clone = hist.Clone()
    hist_clone.SetName(hist_new_name)
    hist_clone.Reset()
    hist_clone.Write()
    return hist_clone


import jinja2
def process_yaml_file(yaml_file, defines=None):

    env = Environment(loader=FileSystemLoader(os.path.dirname(yaml_file)))
    template = env.get_template(os.path.basename(yaml_file))
    output = template.render(**context)
    # process the yaml file using jinja2
    print('defines:', defines)
    new_file = yaml_file.replace('.yaml', '_processed.yaml')
    with open(new_file, 'w') as f:
        f.write(output)
    return new_file

def evaluate_min_max(sexpr1dim, tree, condition=''):
    _htmp_name = f'htmp_{sexpr1dim}'
    tree.Draw(f'{sexpr1dim}>>{_htmp_name}', condition, 'goff')
    h = ROOT.gDirectory.Get(_htmp_name)
    rval = (h.GetXaxis().GetXmin(), h.GetXaxis().GetXmax())
    h.Delete()
    return rval

def hsetup_from_config(hist_config, axis='x', hsetup=None, tree=None):
    if hsetup is None:
        hsetup = {}
    hsetup['clone'] = False  # default
    if len(hist_config['var'].split(':')) > 2:
        print('Invalid var:', hist_config['var'])
        print(' - note: only 1D and 2D histograms are supported')
        sys.exit(1)
    if len(hist_config['var'].split(':')) > 1:
        hsetup['y'] = True
    else:
        hsetup['y'] = False
        if axis == 'y':        
            return hsetup
    k_bins = axis + 'bins'
    k_nbins = axis + 'nbins'
    k_range = axis + 'range'
    k_log = axis + 'log'
    k_svar = hist_config['var'].split(':')[0]
    k_cond = axis + 'cond'
    if 'cond' in hist_config:
        hsetup[k_cond] = hist_config['cond']
    if axis == 'x' and len(hist_config['var'].split(':')) == 2:
        k_svar = hist_config['var'].split(':')[1]
    if k_nbins in hist_config:
        if isinstance(hist_config[k_nbins], int):
            hsetup[k_nbins] = hist_config[k_nbins]
    else:
        hsetup[k_nbins] = 20
    if k_range in hist_config:
        hsetup[k_range] = (hist_config[k_range][0], hist_config[k_range][1])
    else:
        hsetup[k_range] = evaluate_min_max(k_svar, tree)
    hsetup[k_bins] = None
    if k_bins in hist_config:
        hsetup[axis] = hist_config[k_bins]
        if isinstance(hist_config[k_bins], list):
            hsetup[k_nbins] = len(hist_config[k_bins])-1
            hsetup[k_bins] = hist_config[k_bins]
            return hsetup
        if ':' in hist_config[k_bins]:
            hsetup['clone'] = hist_config[k_bins]
            return hsetup
        if 'auto' in hist_config[k_bins]:
            # hsetup[k_range] = (tree.GetMinimum(k_svar), tree.GetMaximum(k_svar))
            hsetup[k_range] = evaluate_min_max(k_svar, tree)
        if 'log' in hist_config[k_bins]:
            hsetup[k_log] = True
            _nbins = hsetup[k_nbins]
            hsetup[k_bins] = logbins(float(hsetup[k_range][0]), float(hsetup[k_range][1]), int(_nbins))
            hsetup[k_nbins] = _nbins
    return hsetup


def build_hist(hist_config, hist_name, tree, output_file):
    hsetup = {}
    hsetup = hsetup_from_config(hist_config, 'x', hsetup, tree)
    hsetup = hsetup_from_config(hist_config, 'y', hsetup, tree)
    print('hsetup:', hsetup)
    if hsetup['clone']:
        hist = get_hist_clone_from_file(hsetup['clone'], hist_name, output_file)
        return hist
    if hsetup['y']:
        if hsetup['xbins'] and hsetup['ybins']:
            hist = ROOT.TH2F(hist_name, hist_name, hsetup['xnbins'], hsetup['xbins'], hsetup['ynbins'], hsetup['ybins'])
        if hsetup['xbins'] and not hsetup['ybins']:
            print(hist_name, hist_name, hsetup['xnbins'], hsetup['xbins'], hsetup['ynbins'], hsetup['yrange'][0], hsetup['yrange'][1])
            hist = ROOT.TH2F(hist_name, hist_name, hsetup['xnbins'], hsetup['xbins'], hsetup['ynbins'], hsetup['yrange'][0], hsetup['yrange'][1])
        if not hsetup['xbins'] and hsetup['ybins']:
            hist = ROOT.TH2F(hist_name, hist_name, hsetup['xnbins'], hsetup['xrange'][0], hsetup['xrange'][1], hsetup['ynbins'], hsetup['ybins'])
        if not hsetup['xbins'] and not hsetup['ybins']:
            hist = ROOT.TH2F(hist_name, hist_name, hsetup['xnbins'], hsetup['xrange'][0], hsetup['xrange'][1], hsetup['ynbins'], hsetup['yrange'][0], hsetup['yrange'][1])
    else:
        if hsetup['xbins']:
            hist = ROOT.TH1F(hist_name, hist_name, hsetup['xnbins'], hsetup['xbins'])
        else:
            hist = ROOT.TH1F(hist_name, hist_name, hsetup['xnbins'], hsetup['xrange'][0], hsetup['xrange'][1])
    return hist


def main(args):
    config = None
    # Load the YAML file
    yaml_file = 'tdraw_config.yaml'
    if args.config:
        yaml_file = args.config
    if os.path.exists(yaml_file) == False:
        print('[e] Could not find file:', yaml_file)
        sys.exit(1)
    yaml_file_processed = process_yaml_file(yaml_file, defines=args.define)
    with open(yaml_file_processed, 'r') as stream:
        try:
            config = yaml.safe_load(stream)
        except yaml.YAMLError as exc:
            print(exc)
    if config is None:
        print('Could not load config from', yaml_file)
        sys.exit(1)
    print('config:', config)
    # Loop over all sections in the config
    for section in config:
        # Get the output file and mode from the first item in the current section
        output_file_name = config[section][0]['output']
        output_mode = config[section][0]['mode']

        if output_file_name == '' or output_file_name is None:
            print('No output file specified for section', section)
            continue
        # Open the output file
        output_file = ROOT.TFile(output_file_name, output_mode)
        known_hists = []

        # Loop over the rest of the items in the current section
        for item in config[section][1:]:
            # Each item is a dictionary where the key is the histogram name and the value is the histogram configuration
            for hist_name, hist_config in item.items():
                hist_name_config = hist_name
                hist_name = hist_name + '_' + section
                # Open the ROOT file and get the tree
                if hist_config['file'] == '' or hist_config['file'] is None:
                    print('No file specified for', hist_name)
                    continue
                root_file = ROOT.TFile.Open(hist_config['file'])
                tree = root_file.Get(hist_config['tree'])
                print('tree:', type(tree), 'hist_config:', hist_config)
                if isinstance(tree, ROOT.TNtuple) or isinstance(tree, ROOT.TTree) or isinstance(tree, ROOT.TChain):
                    print(f'type check for {tree} passed')
                    pass
                else:
                    print('Could not find tree:', hist_config['tree'])
                    continue

                # Create the histogram
                output_file.cd()
                hist = None
                if 'ybins' in hist_config:
                    print ('ybins:', hist_config['ybins'])
                    
                hist = build_hist(hist_config, hist_name, tree, output_file)

                if hist is None:
                    print('[w] Could not create histogram:', hist_name)
                    continue
                known_hists.append((hist_name, hist_name_config, hist))

                # Draw the variable with the specified condition
                draw_string = "{}>>{}".format(hist_config['var'], hist_name)
                tree.Draw(draw_string, hist_config['cond'])
                if hist:
                    print('draw_string:', draw_string, 'w/ cond:', hist_config['cond'], 'entries:', hist.GetEntries())
                else:
                    print('[w] no histogram:', hist_name)
                    continue

                # If 'scale' is specified, scale the histogram
                if 'scale' in hist_config:
                    print('scaling by', hist_config['scale'])
                    _scale_eval_str = hist_config['scale']
                    locals = {}
                    locals_str = {}
                    if isinstance(hist_config['scale'], str): 
                        for hx in known_hists:
                            # print(hx[1], '->?', hx[0])
                            if hx[1] in hist_config['scale']:
                                locals[hx[1]] = hx[2]
                                locals_str[hx[1]] = hx[0]
                        scale = 1.
                        try:
                            scale = eval(_scale_eval_str, {'__builtins__': None}, locals)
                            # scale = eval(hist_config['scale'], {'__builtins__': None}, {'h_jet_pt': ROOT.gDirectory.Get(hist_name)})
                        except Exception as e:
                            print('Could not evaluate scale:', e, 'for', _scale_eval_str, 'using', locals_str)
                            sys.exit(1)
                    else:
                        scale = hist_config['scale']
                    print(f'scaling {hist.GetName()} by', scale, 'using', locals_str)
                    hist.Scale(scale)

                if 'bw' in hist_config:
                    if hist_config['bw'] == True:
                        hist.Scale(1., 'width')
                print('---------------------------------')
                # Save the histogram to the output file
                hist.Write()

        # Close the output file
        output_file.Close()
    
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Draw histograms from a YAML file')
    parser.add_argument('-c', '--config', help='YAML config file')
    parser.add_argument('-d', '--define', help='Define a variable x=y to replace {{x}} by y in the config', nargs='+')
    parser.add_argument('-x', '--context', help='the context file', default=None)
    args = parser.parse_args()
    print(args)

    context = {}
    if args.context:
        with open(args.context, 'r') as f:
            context = eval(f.read())

    if args.define:
        for define in args.define:
            key, value = define.split('=')
            context[key] = value

    args.define = context
    print(args)
    main(args)