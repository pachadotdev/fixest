#!/usr/bin/bash

# Function to debug R scripts
# Adapted from https://mpadge.github.io/blog/blog012.html
# Copy the function to ~/.bashrc source ~/.bashrc again
# Then run debugr from the command line

function debugr () {
    # if no argument is provided, ask for a file
    if [ -z "$1" ]; then
        read -p "Enter the script to debug: " script
    else
        script=$1
    fi

    # if the file does not exist, exit
    if [ ! -f "$script" ]; then
        echo "File $script does not exist"
        return 1
    fi

    # if the file does not end in .R, exit
    if [[ "$script" != *.R ]]; then
        echo "File $script does not end in .R"
        return 1
    fi

    # run R in debug mode, but before that we compiled with debug symbols
    # see https://reside-ic.github.io/blog/debugging-memory-errors-with-valgrind-and-gdb/
    Rscript -e "pkgbuild::compile_dll()"
    R -d 'valgrind -s --leak-check=full --show-leak-kinds=all --track-origins=yes' -f $script
}
