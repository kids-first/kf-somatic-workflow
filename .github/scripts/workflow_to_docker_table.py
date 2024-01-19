#!/usr/bin/env python3

import argparse
import pathlib

def parse_command():
    """Function to parse the command line input
    Args:
        None
    Return:
        Namespace: returns the args as a standard Namespace object
    """
    parser = argparse.ArgumentParser(description='Generate Table of Tools and Their Dockers from a CWL Workflow')
    parser.add_argument('input_file')
    parser.add_argument('--output_file',
        default='TABLE.md',
        help='Path to output file.')
    args = parser.parse_args()
    if not args.output_file.endswith('.md'): args.output_file = args.output_file + ".md"
    return args

def parse_workflow(filepath: str, pathdict: dict | None = None) -> dict:
    """Function to recursively build dict of unique subworkflows and tools from a file.
    Args:
        filepath: Path to the input cwl
        pathdict: Dict to store the unique subworkflows and tools
    Return:
        pathdict: Dict containing two lists. One of unique subworkflows; one of unique tools
    """ 
    if not pathdict: pathdict = {}
    if "workflows" not in pathdict: pathdict["workflows"] = set()
    if "tools" not in pathdict: pathdict["tools"] = set()
    inpath = pathlib.Path(filepath).resolve()
    with open(inpath) as fh:
        for line in fh:
            if not line.strip().startswith('run:'): continue
            if not line.strip().endswith('.cwl'): continue
            relpath = inpath.parent / line.strip().split(' ')[-1]
            respath = relpath.resolve()
            if respath in pathdict["workflows"] or respath in pathdict["tools"]: continue
            if not respath.match('*/tools/*'):
                pathdict["workflows"].add(respath)
                pathdict = parse_workflow(respath, pathdict)
            else:
                pathdict["tools"].add(respath)
    return pathdict

def get_docker(filepath: str) -> str:
    """
    Args:
        filepath: Path to the input cwl
    Return:
        string: Name of item in dockerPull field or None
    """
    inpath = pathlib.Path(filepath)
    with open(inpath) as toolfile:
        for line in toolfile:
            if "dockerPull" in line:
                return line.strip().split(' ')[-1].strip("\"'")
    return "None"

def main():
    args = parse_command()
    appdict: dict = parse_workflow(args.input_file)
    payload: list[tuple[str]] = sorted([(cwltool.name, get_docker(cwltool)) for cwltool in appdict['tools']])
    with open(args.output_file, 'w') as outfile:
        print(f"# Dockers of {pathlib.Path(args.input_file).name}\n", file=outfile)
        print("TOOL|DOCKER\n-|-", file=outfile)
        for item in payload:
            print('|'.join(item), file=outfile)

if __name__ == '__main__':
    main()
