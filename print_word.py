import sys

# Get the word from command-line arguments
if len(sys.argv) > 1:
    word = sys.argv[1]
    print(f"The word you provided is: {word}")
else:
    print("No word was provided.")
