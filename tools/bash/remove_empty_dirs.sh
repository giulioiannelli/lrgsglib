#!/bin/bash

# Function to display help information
display_help() {
    echo "Usage: $0 [OPTION]... [DIRECTORY]"
    echo "Recursively remove all empty directories from a specified DIRECTORY."
    echo ""
    echo "  -h, --help      display this help and exit"
    echo "  -v, --verbose   print names of deleted directories"
    echo ""
    echo "Example:"
    echo "  $0 /path/to/directory   Removes all empty directories starting from '/path/to/directory'."
    echo ""
    echo "Report bugs to: your-email@example.com"
}

# Variables
verbose=0

# Check for no arguments
if [ "$#" -eq 0 ]; then
    echo "Error: No directory provided."
    echo "Try '$0 --help' for more information."
    exit 1
fi

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
        *)
            DIRECTORY=$1  # Assume the first argument is the directory path
            break
            ;;
    esac
done

# Check if the provided argument is a directory
if [ ! -d "$DIRECTORY" ]; then
    echo "Error: Directory does not exist."
    exit 1
fi

# Recursively find and delete empty directories
if [ "$verbose" -eq 1 ]; then
    find "$DIRECTORY" -type d -empty -print -delete
    echo "Verbose mode: Printed names of all deleted directories."
else
    find "$DIRECTORY" -type d -empty -delete
fi

# Inform the user of completion
echo "Empty directories have been removed from $DIRECTORY."
