#!/bin/bash

# Function to display help information
display_help() {
    echo "Usage: $0 [OPTION]... [DIRECTORY]"
    echo "Convert all PDF files in a specified DIRECTORY to PNG with optional transparent background."
    echo ""
    echo "  -h, --help         display this help and exit"
    echo "  -v, --verbose      print names of processed files"
    echo "  -t, --transparent  make white background transparent"
    echo ""
    echo "Example:"
    echo "  $0 -t /path/to/directory   Converts all PDF files in '/path/to/directory' with transparent background."
    echo ""
    echo "Report bugs to: your-email@example.com"
}

# Variables
verbose=0
transparent=0
DIRECTORY=""

# Parse command line options
while [ "$#" -gt 0 ]; do
    case "$1" in
        -h|--help)
            display_help
            exit 0
            ;;
        -v|--verbose)
            verbose=1
            shift
            ;;
        -t|--transparent)
            transparent=1
            shift
            ;;
        *)
            DIRECTORY=$1  # Assume the first argument is the directory path
            break
            ;;
    esac
done

# Set the default directory to current directory if none is provided
if [ -z "$DIRECTORY" ]; then
    DIRECTORY=$(pwd)
fi

# Check if the provided argument is a valid directory
if [ ! -d "$DIRECTORY" ]; then
    echo "Error: $DIRECTORY is not a valid directory."
    exit 1
fi

# Change to the specified directory
cd "$DIRECTORY" || exit

# Process each PDF file in the directory
for pdf in *.pdf; do
    # Check if there are no PDF files
    if [ ! -e "$pdf" ]; then
        [ "$verbose" -eq 1 ] && echo "No PDF files found in the directory."
        exit 0
    fi

    # Convert PDF to PNG using Ghostscript
    gs -sDEVICE=pngalpha -o "${pdf%.pdf}.png" -r400 "$pdf"

    # Make white background transparent using ImageMagick, if requested
    if [ "$transparent" -eq 1 ]; then
        convert "${pdf%.pdf}.png" -transparent white "${pdf%.pdf}.png"
    fi

    [ "$verbose" -eq 1 ] && echo "Processed: $pdf"
done

[ "$verbose" -eq 1 ] && echo "All PDF files have been processed."
