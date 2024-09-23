import pandas as pd
import os
from sklearn.metrics.pairwise import cosine_similarity
import re
import json
from pathlib import Path

file_path = "/Users/sakshmenon/Desktop/results/workflow files/"
os.chdir(file_path)

ref_file = file_path + "s_cerevisiae_protein.csv"
query_file = file_path + "p_carinii_protein.csv"
top_match_file = file_path + "p_carinii2s_cerevisiae.bfs.k3.txt"

def extract_queries(file):
    pattern = r"Query:\s+>([\w.]+)\s.*?\nTop1:\s+[\d.]+\s+>([\w.]+)\s.*?\nTop2:\s+[\d.]+\s+>([\w.]+)\s.*?\nTop3:\s+[\d.]+\s+>([\w.]+)\s"
    
    with open(file) as obj:
        queries = {}
        file_content = obj.read()
        query_top1_pairs = re.findall(pattern, file_content)
        query_top_dict = v = {i[0]: i[1:] for i in query_top1_pairs}
        for id in query_top_dict:
            if id not in queries:
                queries.update({id: query_top_dict[id]})
            else:
                queries[id].append(query_top_dict[id])

        query_IDs = list(queries.keys())
        t_pattern = r"Time:\s+([\d.]+)"
        times = re.findall(t_pattern, file_content)

        pairs = zip(query_IDs, times)
        q_times = {}
        for pair in pairs:
            q_times[pair[0]] = pair[1]

    return queries, q_times

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

def score_gen(qindex_df, rindex_df, query_df, ref_df, queries, pairs):
    score_report = []

    for query in qindex_df:
        qindex = qindex_df[qindex_df == query].index[0]
        query_vector = query_df.loc[qindex].values[:-1]
        matches = []
        for inst in queries[query]:
            rindex = rindex_df[rindex_df == inst].index[0]
            ref_vector = ref_df.loc[rindex].values[:-1]
            score = cosine_similarity([query_vector], [ref_vector])
            matches.append(({"ID": inst, "Score": round(float(score[0]), ndigits=5)}))

        sub_report = {'Query': ({"ID": query, "Time" : float(pairs[query]), 'Matches': matches})}
        score_report.append(sub_report)
    
    return score_report

def run_score_report():
    match_file = Path(top_match_file)
    query_f = Path(query_file)
    ref_f = Path(ref_file)
    if match_file.is_file():
        if query_f.is_file():
            if ref_f.is_file():
                queries, pairs = extract_queries(top_match_file)
                query_df, qindex_df, ref_df, rindex_df = query_ref_df_gen(query_file, ref_file)
                score_report = score_gen(qindex_df, rindex_df, query_df, ref_df, queries, pairs)
                score_report

                os.chdir('/Users/sakshmenon/Desktop/PLM/Workflow Gen/Report Gen')

                json_object = json.dumps(score_report, indent=4)
                with open("sample.json", "w") as outfile:
                    outfile.write(json_object)
                print('Json File Dumped At: /Users/sakshmenon/Desktop/PLM/Workflow Gen/Report Gen')
            else:
                print('Reference File Not Found')
        else:
            print('Query File Not Found')
    else:
        print('Top Match File Not Found')