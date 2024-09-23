# The script makes search for similar annotated proteins using sequence embeddings.
# THe goal of this script is to use different search algorithms to benhcmark their performance.
# Make sure the following libraries installed:
# pip install numpy pandas scikit-learn scipy umap-learn hnswlib faiss-cpu
# load anaconda module before running the script:
# module load anaconda3

import os
import sys
import argparse
import time

import numpy as np
import pandas as pd
from sklearn.neighbors import BallTree
from sklearn.metrics import pairwise_distances
from sklearn.metrics.pairwise import cosine_similarity
from scipy.spatial.distance import cdist
# import hnswlib
# import faiss
from scipy import stats
import json
import re

import warnings
warnings.filterwarnings("ignore")

parser = argparse.ArgumentParser(description="Run the protein similarity search (brute force, ball tree, graph, faiss). Example command line: python search_similar_protein.py -rf PATH/reference_embeddings.csv -qf PATH/query_embedding.csv -m bfs -k 3")
parser.add_argument('--ref_file','-rf', type=str, help='input csv filename with reference embeddings')
parser.add_argument('--q_file','-qf', type=str, help='input csv filename with query embeddings')
parser.add_argument('--method', '-m', type=str, help='search method: Brute Force, Ball Tree, HNSW Graph, FAISS.', choices=['bfs', 'bts', 'graph', 'faiss'], default='bfs')
parser.add_argument('--out_file', '-of', type=str, help='Output embedding filename (do not input whole .csv, for example, path/folder/filename is enough. This is for creating multiple file for both protein and residue level.', default='bfs')
parser.add_argument('--top_k', '-k', type=int, default=3, help='retrieve top k similar proteins.')
parser.add_argument('--format', '-f', type=str, help='output format.', choices=['json', 'txt'])

# Load protein embedding data
def load_protein_embeddings(file_path):
    data = pd.read_csv(file_path, skiprows=1, header=None)
    embeddings = data.iloc[:, :-1].values  # All columns except the last (protein ID)
    protein_ids = data.iloc[:, -1].values  # The last column (protein ID)
    return embeddings, protein_ids


# Brute-Force Search
def brute_force_search(query_embedding, reference_embeddings, top_k):
    distances = pairwise_distances(query_embedding.reshape(1, -1), reference_embeddings, metric='euclidean')
    top_indices = np.argsort(distances)[0][:top_k]
    return top_indices, distances[0][top_indices]


# Ball Tree Search
def ball_tree_init(reference_embeddings):
    return BallTree(reference_embeddings, metric='euclidean')

def ball_tree_search(ball_tree, query_embedding, top_k):
    distances, indices = ball_tree.query(query_embedding.reshape(1, -1), k=top_k)
    return indices[0], distances[0]


# HNSW Graph Search
def hnsw_init(reference_embeddings):
    dim = reference_embeddings.shape[1]
    hnsw_index = hnswlib.Index(space='l2', dim=dim)
    hnsw_index.init_index(max_elements=reference_embeddings.shape[0], ef_construction=200, M=16)
    hnsw_index.add_items(reference_embeddings)
    return hnsw_index

def hnsw_search(hnsw_index, query_embedding, top_k):
    hnsw_index.set_ef(top_k * 10)  # Sets the recall effort parameter
    indices, distances = hnsw_index.knn_query(query_embedding, k=top_k)
    return indices[0], distances[0]


# FAISS Search
def faiss_init(reference_embeddings):
    dim = reference_embeddings.shape[1]
    faiss_index = faiss.IndexFlatL2(dim)
    faiss_index.add(reference_embeddings)
    return faiss_index

def faiss_search(faiss_index, query_embedding, top_k):
    distances, indices = faiss_index.search(query_embedding.reshape(1, -1), top_k)
    return indices[0], distances[0]


# Compute confidence intervals for distances
def compute_confidence_interval(distances):
    mean_dist = np.mean(distances)
    std_dist = np.std(distances, ddof=1)
    conf_int = stats.t.interval(0.95, len(distances) - 1, loc=mean_dist, scale=std_dist / np.sqrt(len(distances)))
    return mean_dist, conf_int

def query_ref_df_gen(query_file, ref_file):
    query_df = pd.read_csv(query_file)
    qindex_df = query_df['ProteinID']

    pattern = r">([\w.]+)"
    for index in qindex_df.index:
        query_id = re.findall(pattern, query_df['ProteinID'][index])
        qindex_df[index] = query_id[0]
    
    ref_df = pd.read_csv(ref_file)
    rindex_df = ref_df['ProteinID']

    pattern = r">([\w.]+)"
    for index in rindex_df.index:
        ref_id = re.findall(pattern, ref_df['ProteinID'][index])
        rindex_df[index] = ref_id[0]

    return query_df, qindex_df, ref_df, rindex_df

# Main similarity search function
def similarity_search(query_embeddings, query_protein_ids, reference_embeddings, reference_protein_ids, method, top_k, format, reference_file, query_file, query_df, qindex_df, ref_df, rindex_df):
    # Initialize the search model
    search_model = None
    if method == 'bts':
        search_model = ball_tree_init(reference_embeddings)
    elif method == 'graph':
        search_model = hnsw_init(reference_embeddings)
    elif method == 'faiss':
        search_model = faiss_init(reference_embeddings)
    else:
        search_model = None
    
    methods = {
        'bfs': {"Brute-Force Search": brute_force_search},
        'bts': {"Ball Tree": ball_tree_search},
        'graph': {"HNSW Graphs": hnsw_search},
        'faiss': {"FAISS": faiss_search}
    }

    pat = r">([\w.]+)\s.*?"
    if format == 'json':
        flag = 1
        output_file = reference_file + '_' + query_file + '_' + method + '.json'

    else:
        output_file = reference_file + '_' + query_file + '_' + method + '.txt'
        os.chdir('/Users/sakshmenon/Desktop/PLM/Workflow Gen/Report Gen')
        output_file_obj = open(output_file, 'w')
        flag = 0

    method_report = {method : []}
    
    for method_name, search_func in methods[method].items():
        if format == 'txt':
            output_file_obj.write(f"Method:\t{method_name}\n")
        total_time = 0

        score_report = []

        for query_id, query_embedding in zip(query_protein_ids, query_embeddings):
            sub_report = {'Query' : {"ID": '', "Time": '', "Matches" : []}}
            start_time = time.time()
            if search_model is None:
                top_indices, distances = search_func(query_embedding, reference_embeddings, top_k)
            else:
                top_indices, distances = search_func(search_model, query_embedding, top_k)
            end_time = time.time()
            
            query_time = end_time - start_time
            total_time += query_time
#            mean_dist, conf_int = compute_confidence_interval(distances)
            query = re.findall(pat, query_id)[0]
            if flag:
                sub_report['Query']['ID'] = query
            else:
                output_file_obj.write(f"\nQuery:\t{query_id}\n")

#            print(f"Top {top_k} similar proteins (distance,match):")
            for i, (idx, dist) in enumerate(zip(top_indices, distances)):
                match = re.findall(pat, reference_protein_ids[idx])[0]


                qindex = qindex_df[qindex_df == query].index[0]
                query_vector = query_df.loc[qindex].values[:-1]

                rindex = rindex_df[rindex_df == match].index[0]
                ref_vector = ref_df.loc[rindex].values[:-1]

                score = cosine_similarity([query_vector], [ref_vector])

                if flag:
                    sub_report['Query']['Matches'].append(({"ID": match, "Score": round(float(score[0]), ndigits=5)}))

                else:
                    output_file_obj.write(f"Top{i+1}:\t{dist:.3f}\t{reference_protein_ids[idx]}\n")

#            print(f"  Mean distance: {mean_dist:.3f}, Confidence interval: {conf_int}")

            if flag:
                sub_report['Query']['Time'] = round(query_time, 6)
                score_report.append(sub_report)
            else:
                output_file_obj.write(f"Time:\t{query_time:.5f}\n")
        
        if flag:
            json_object = json.dumps(score_report, indent=4)
            os.chdir('/Users/sakshmenon/Desktop/PLM/Workflow Gen/Report Gen')
            with open(output_file, "w") as outfile:
                outfile.write(json_object)
        
        else:
            output_file_obj.write(f"\nTotal query proteins={len(query_embeddings)}, reference proteins={len(reference_embeddings)}")
            output_file_obj.write(f"Total search time for {method_name}: {total_time:.5f} seconds")
            avg_time = total_time / len(query_embeddings)
            output_file_obj.write(f"Average search time for {method_name}: {avg_time:.5f} seconds")

            output_file_obj.close()

# Example Usage
if __name__ == "__main__":
    # passing arguments 
    args = parser.parse_args()
    if len(sys.argv) < 2:
        parser.print_help()
        sys.exit()
        
    # Check if files were provided as arguments; if not, prompt for them
    if args.ref_file:
        reference_file = args.ref_file
    else:
        print("Error: the reference embeddings file not specified")
        parser.print_help()
        sys.exit(1)
    if args.q_file:
        query_file = args.q_file
    else:
        print("Error: the query embeddings file not specified")
        parser.print_help()
        sys.exit(1)

    # Ensure the files exist before proceeding
    if not os.path.exists(reference_file):
        print(f"Error: the reference embeddings file '{reference_file}' not found")
        sys.exit(1)
    if not os.path.exists(query_file):
        print(f"Error: the query embeddings file '{query_file}' not found")
        sys.exit(1)

    # Process other optional parameters
    search_method = args.method
    top_k = args.top_k
    if args.format:
        format = args.format
    else:
        print("Error: output format not specified")
        parser.print_help()
        sys.exit(1)

    # Load reference and query data
    reference_embeddings, reference_protein_ids = load_protein_embeddings(reference_file)
    query_embeddings, query_protein_ids = load_protein_embeddings(query_file)
    
    # Perform similarity search
    query_df, qindex_df, ref_df, rindex_df = query_ref_df_gen(query_file, reference_file)
    similarity_search(query_embeddings, query_protein_ids, reference_embeddings, reference_protein_ids, search_method, top_k, format, reference_file, query_file, query_df, qindex_df, ref_df, rindex_df)

    # reference_embeddings, reference_protein_ids = load_protein_embeddings('/Users/sakshmenon/Desktop/results/workflow files/s_cerevisiae_protein.csv')
    # query_embeddings, query_protein_ids = load_protein_embeddings('/Users/sakshmenon/Desktop/results/workflow files/p_carinii_protein.csv')
    
    # query_df, qindex_df, ref_df, rindex_df = query_ref_df_gen('/Users/sakshmenon/Desktop/results/workflow files/p_carinii_protein.csv', '/Users/sakshmenon/Desktop/results/workflow files/s_cerevisiae_protein.csv')
    # similarity_search(query_embeddings, query_protein_ids, reference_embeddings, reference_protein_ids, 'bfs', 3, 'txt', '/Users/sakshmenon/Desktop/results/workflow files/s_cerevisiae_protein.csv', '/Users/sakshmenon/Desktop/results/workflow files/p_carinii_protein.csv', query_df, qindex_df, ref_df, rindex_df)
