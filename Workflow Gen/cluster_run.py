import pandas as pd
import csv
import sys
import argparse

from bio_embeddings.embed import SeqVecEmbedder as SeqVecEmbedder
from bio_embeddings.embed import ESM1bEmbedder as Esm1bEmbedder
from bio_embeddings.embed import ESMEmbedder as EsmEmbedder
from bio_embeddings.embed import ProtTransT5XLU50Embedder as Prott5Embedder
from bio_embeddings.embed import OneHotEncodingEmbedder as OneHotEmbedder

import warnings
warnings.filterwarnings("ignore")

parser = argparse.ArgumentParser(description="Run embedding models (esm1b, esm, seq, prott5, onehot) at either `protein` level or `residue` level. Example command line: python cluster_run.py -ff PATH/fasta_file.txt -e seq -of PATH/embedding_filename -l residue")
parser.add_argument('--file','-ff', type=str, help='input fasta filename to run embedding models on')
parser.add_argument('--embed', '-e', type=str, help='Embedding models args: esm1b, esm, seq, prott5, onehot.', choices=['esm1b', 'esm', 'seq', 'prott5', 'onehot'])
parser.add_argument('--out_file', '-of', type=str, help='Output embedding filename (do not input whole .csv, for example, path/folder/filename is enough. This is for creating multiple file for both protein and residue level.')
parser.add_argument('--level', '-l', type=str, default='protein', help='arg: protein, residue, or both, default protein. Return embeddings at either protein level or residue level.', choices=['protein','residue', 'both'])

# takes file and returns dict of {name: seq}
def read_file(filepath):
    '''
    Take a fasta file and return a dictionary of {pID: sequence}
    Fasta file should each have 1 line for pid and 1 line for sequence
    Input: str filepath - file to read
    Output: dictionary {str: str}
    '''
    dictionary = {}
    with open(filepath) as fasta_file:
        seq = ''
        for line in fasta_file:
            line=line.rstrip()
            if line.startswith(">"):
                if seq.__len__():
                    dictionary[name] = seq
                name = line
                seq = ''
            else:
                seq = seq+line
        dictionary[name] = seq
        
    dic2=dict(sorted(dictionary.items(),key= lambda x:len(x[1]), reverse=True))
    return dic2

# initialize save file
def file_initialize(target_file, vect_len, is_residue=False):
    '''
    Initialize file with the right column names
    If residue mode is on, initialized both protein and residue
    level file
    If not, initialized only protein file
    Input: str target_file - file name to initialized
    int vect_len - vector length (specific to different models)
    Output: True if successfully run
    '''
    random_list = ['sv'+str(x) for x in range(1, vect_len+1)]
    random_list.append('ProteinID')
    if is_residue:
        random_list.append('ResidueID')
    # writing to csv file 
    with open(target_file, 'w') as csvfile: 
        # creating a csv writer object 
        csvwriter = csv.writer(csvfile) 

        # writing the fields 
        csvwriter.writerow(random_list)
    return True

# writing to csv file
def write_csv(target_file, row):
    '''
    A function that take embeddings dict of name: embeddings
    and update into a new file
    Input: str target_file - file to write on
    Output: True if successfully run
    '''
    with open(target_file, 'a') as csvfile: 
        # creating a csv writer object 
        csvwriter = csv.writer(csvfile) 
        # writing the fields 
        csvwriter.writerow(row)
    return True

# takes dict (name, seq) and returns dict (name, raw embed)
def embed_one_file(sequences, embed, out_file, vect_len, is_residue=False):
    '''
    A function that return the embeddings of a dictionary of sequences
    Input: dictionary of {name: sequence}, str embed - model type,
    str out_file - file to write embeddings to,
    int vect_len - vect_len specific to model type
    bool is_residue - if it is residue, generate 2 files for both protein and residue embeddings,
    generate 1 file for protein embeddings only if not
    Output: 0 if run successfully
    '''
    if embed == 'esm1b':
        embedder = Esm1bEmbedder()
    elif embed == 'esm':
        embedder = EsmEmbedder()
    elif embed == 'seq':
        embedder = SeqVecEmbedder()
    elif embed == 'prott5':
        embedder = Prott5Embedder()
    elif embed == 'onehot':
        embedder = OneHotEmbedder()
        
    for protein in sequences:
        reduced_embeddings = []
        
        embeddings = embedder.embed(sequences[protein]) # using single embed since seach element is a STRING
        if not(is_residue==1):
            # protein level
            reduced_embeddings = embedder.reduce_per_protein(embeddings)
            reduced_embeddings = reduced_embeddings.tolist()
            reduced_embeddings.append(protein)
            write_csv(out_file+'_protein.csv', reduced_embeddings)
            reduced_embeddings.clear()
        
        if not(is_residue==0):
        # residue level
            r_embeddings = embeddings
            if embed == 'seq':
                r_embeddings = r_embeddings.sum(0)
            residue_embeddings = r_embeddings.tolist()
            for idx, char in enumerate(sequences[protein]):
                res_embed = residue_embeddings[idx]
                res_embed.append(protein)
                res_embed.append(char+str(idx+1))
                write_csv(out_file+'_residue.csv', res_embed)
                # res_embed.clear()
    return 0

if __name__ == "__main__":
    # passing arguments
    args = parser.parse_args()
    if not len(sys.argv) > 1:
        parser.print_help()
        sys.exit()
    # if not enough arguments, then return help with error
    if not len(sys.argv) > 4:
        print('Missing argument(s)')
        parser.print_help()
        sys.exit()
    # reading sequence
    seq_dict = read_file(filepath=args.file)
    if args.embed == 'onehot':
        vect_len = 21
    elif ((args.embed == 'esm') | (args.embed == 'esm1b')):
        vect_len = 1280
    else:
        vect_len = 1024
    # checking for level
    if args.level=='protein':
        l = 0
    elif args.level=='residue':
        l = 1
    elif args.level=='both':
        l = 2
    # file initialization
    if not(l==0): #create two files for both option
        file_initialize(args.out_file+'_residue.csv',vect_len,l)
    if not(l==1):
        file_initialize(args.out_file+'_protein.csv',vect_len)
    # Running embeddings/writing to csv row by row
    embed_one_file(sequences=seq_dict, embed=args.embed, out_file=args.out_file, vect_len=vect_len, is_residue=l)

    