#!/usr/bin/bash

# Function to debug R scripts
# Adapted from https://mpadge.github.io/blog/blog012.html

function r_debug_symbols () {
    # if src/Makevars does not exist, exit
    if [ ! -f src/Makevars ]; then
        echo "File src/Makevars does not exist"
        return 1
    fi

    # if src/Makevars contains a line that says "PKG_CPPFLAGS"
    # but there is no "-UDEBUG -g" on it
    # then add "PKG_CPPFLAGS += -UDEBUG -g" at the end
    if grep -q "PKG_CPPFLAGS" src/Makevars; then
        if ! grep -q "PKG_CPPFLAGS.*-UDEBUG.*-g" src/Makevars; then
            echo "PKG_CPPFLAGS += -UDEBUG -g" >> src/Makevars
        fi
    fi

    # if src/Makevars does not contain a line that reads
    # PKG_CPPFLAGS ...something... -UDEBUG -g ...something...
    # then add PKG_CPPFLAGS = -UDEBUG -g to it
    if ! grep -q "PKG_CPPFLAGS.*-UDEBUG.*-g" src/Makevars; then
        echo "PKG_CPPFLAGS = -UDEBUG -g" >> src/Makevars
    fi
}

function r_valgrind () {
    # if no argument is provided, ask for a file
    if [ -z "$1" ]; then
        read -p "Enter the script to debug: " script
    else
        script=$1
    fi

    # if no output file is provided, use dev/valgrind.txt
    if [ -z "$2" ]; then
        output="dev/valgrind.txt"
    else
        output=$2
    fi

    # if the file does not exist, exit
    if [ ! -f "$script" ]; then
        echo "File $script does not exist"
        return 1
    fi

    # if the file does not end in .R/.r, exit
    shopt -s nocasematch
    if [[ "$script" != *.R ]]; then
        echo "File $script does not end in .R or .r"
        return 1
    fi
    shopt -u nocasematch

    # run R in debug mode, but after that we compiled with debug symbols
    # see https://reside-ic.github.io/blog/debugging-memory-errors-with-valgrind-and-gdb/
    # R -d 'valgrind -s --leak-check=full --show-leak-kinds=all --track-origins=yes' -f $script 2>&1 | tee valgrind.txt
    R -d 'valgrind -s --track-origins=yes' -f $script 2>&1 | tee $output
} 

r_debug_symbols

r_valgrind dev/test_memory_leaks.r dev/valgrind.txt
