#!/usr/bin/env python

import os
import argparse
from jinja2 import Environment, FileSystemLoader

def main():
    parser = argparse.ArgumentParser(description='Process a Jinja2 template.')
    parser.add_argument('template', help='the template file')
    parser.add_argument('-o', '--output', help='the output file', default=None)
    parser.add_argument('-c', '--context', help='the context file', default=None)
    parser.add_argument('-d', '--define', help='Define a variable x=y to replace {{x}} by y in the config', nargs='+')
    args = parser.parse_args()

    context = {}
    if args.context:
        with open(args.context, 'r') as f:
            context = eval(f.read())

    if args.define:
        for define in args.define:
            key, value = define.split('=')
            context[key] = value

    env = Environment(loader=FileSystemLoader(os.path.dirname(args.template)))
    template = env.get_template(os.path.basename(args.template))
    output = template.render(**context)

    if args.output:
        with open(args.output, 'w') as f:
            f.write(output)
    else:
        print(output)

if __name__ == '__main__':
    main()
