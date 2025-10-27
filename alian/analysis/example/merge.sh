#!/usr/bin/bash

show_help() {
    cat << EOF
Usage:
  $0 -i|--input FILE -o|--output OUTPUT [-n|--name FILENAME] [-f|--force] [-h|--help]

Searches for files in a given directory and a given filename,
and merges them with ROOT hadd into a specified output file.
The list of files that were merged will be placed in a text
file in the same directory as the output file.

Options:
  -i, --input INPUT    Input directory to search for files
  -o, --output OUTPUT  Filename for merged output
  -n, --name NAME      Search for NAME in INPUT (default: results.root)
  -j, --jobs [N]       Parallelize merging with N jobs (default: system max)
  -f, --force          Overwrite existing output file (default: false)
  -h, --help           Show this help message

EOF
}

if [ $# -eq 0 ]; then
    show_help && exit 0
fi

PARSEDARGS=$(getopt -o i:o:n:j:fh \
                    --long input:,output:,name:,jobs:,force,help \
                    -n 'merge' -- "$@")
PARSE_EXIT=$?
if [ $PARSE_EXIT -ne 0 ] ; then exit $PARSE_EXIT ; fi

# NOTE: quotes around PARSEDARGS essential to parse spaces between parsed options
eval set -- "$PARSEDARGS"

input=
output=
filename=results.root
force_opt=
njobs=$(nproc)

while true; do
    case "$1" in
        -i | --input )    input="$2"; shift 2;;
        -o | --output )   output="$2"; shift 2;;
        -n | --name )     filename="$2"; shift 2;;
        -j | --jobs )     njobs="$2"; shift 2;;
        -f | --force )    force_opt="-f"; shift ;;
        -h | --help )     show_help; exit 0 ;;
        -- ) shift; break ;;
        * ) break ;;
    esac
done

if [[ -z "$input" ]]; then
    echo "Input directory must be passed:"
    show_help
    exit 1
fi
if [[ -z "$output" ]]; then
    echo "Output filename must be passed:"
    show_help
    exit 1
fi

filelist="$(dirname -- "$(realpath "$0")")/../list.txt"

if [[ -z "$force_opt" ]]; then
    if [[ -f $filelist ]]; then
        echo "Filelist $filelist already exists; pass -f or --force to overwrite."
    fi
    if [[ -f $output ]]; then
        echo "Output file $output already exists; pass -f or --force to overwrite."
    fi
fi

echo "Searching in $input for files named $filename..."
find "$output" -name "$filename" -type f > "$filelist"
echo "Found $(wc -l < "$filelist") files to merge. Merging..."

hadd -v 99 $force_opt -j "$njobs" "$output" @"$filelist"
exitcode=$?
if [[ $exitcode -ne 0 ]]; then
    echo "Merging failed."
    exit $exitcode
fi
echo "Merging completed."