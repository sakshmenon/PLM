#!/bin/bash

# Check if the word is provided as an argument
if [ -z "$1" ]; then
    echo "Please provide a word as an argument."
    exit 1
fi

# Pass the word to the Python script
python print_word.py "$1"
