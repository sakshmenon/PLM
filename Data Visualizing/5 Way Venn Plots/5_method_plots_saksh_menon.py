import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde
from venn import venn
import pandas as pd
import json
import re

def plot_venn_5_sets(sets_dict):
    venn(sets_dict)
    plt.title("4 Algo Venn Diagram")
    plt.show()

def f_ext_30(filename):
    results = {}
    current_query = None

    # Define regex patterns for matching query lines and top match lines
    query_pattern = re.compile(r'^Query:\s+>(\S+)')
    match_pattern = re.compile(r'^Top\d+:\s+(\d\.\d{3})\s+>(\S+)')
    
    with open(filename, 'r') as file:
        for line in file:
            # Check for a query line
            query_match = query_pattern.search(line)
            if query_match:
                current_query = query_match.group(1)
                continue  # Move to the next line

            # Check for a top match line
            match_match = match_pattern.search(line)
            if match_match and current_query:
                score = float(match_match.group(1))
                match_id = match_match.group(2)
                if current_query in results.keys():
                    results[current_query].append({
                        'Match ID': match_id,
                        'Score': score })
                else:
                    results[current_query] = [{
                        'Match ID': match_id,
                        'Score': score}]
                    
                
        main_l = []
        multi_match_list = {}
        for q in results:
            que = q
            t1 = results[q][0]['Match ID']
            t2 = results[q][1]['Match ID']
            t3 = results[q][2]['Match ID']

            score1 = results[q][0]['Score']
            score2 = results[q][1]['Score']
            score3 = results[q][2]['Score']

            score_set = [score1, score2, score3]

            if score_set.count(score1)>1:
                temp = [t1,t2,t3][0:score_set.count(score1)]
                temp.sort()
                multi_match_list[que] = temp

            main_l.append({'Query' : que, 'T1' : t1, 'Score' : score1})

    return pd.DataFrame(main_l), multi_match_list

with open("baker2fission.e1.txt") as fileobj: #file path here
    content = fileobj.read()
    lines = content.split('\n')

linewise = []
queries = []
for line in lines[:-1]:
    ph = line.split('\t')
    queries.append(ph[0])
    linewise.append({'Query' : ph[0], 'Ref' : ph[1], 'Score' : eval(ph[-2])})

df = pd.DataFrame(linewise)
set_dict = []
for query in set(queries):
    temp = df[df['Query'] == query].sort_values('Score')
    temp = temp.iloc[0][['Query', 'Ref', 'Score']].to_dict()
    set_dict.append(temp)

e_df = pd.DataFrame(set_dict)

faiss_json = "s_cerevisiae2s_pombe.faiss.k30.txt" #file path here
graph_json = "s_cerevisiae2s_pombe.graph.k30.txt" #file path here
bts_json = "s_cerevisiae2s_pombe.bts.k30.txt" #file path here
bfs_json = "s_cerevisiae2s_pombe.bfs.k30.txt" #file path here

faiss_df, multi_match_list_faiss = f_ext_30(faiss_json)
graph_df, multi_match_list_graph = f_ext_30(graph_json)
bts_df, multi_match_list_bts = f_ext_30(bts_json)
bfs_df, multi_match_list_bfs = f_ext_30(bfs_json)

faiss_q = [*multi_match_list_faiss]
graph_q = [*multi_match_list_graph]
bts_q = [*multi_match_list_bts]
bfs_q = [*multi_match_list_bfs]

dupe_set = [*set(faiss_q + graph_q + bts_q + bfs_q)]
for q in dupe_set:
    try:
        faiss_m = set(multi_match_list_faiss[q])
        graph_m = set(multi_match_list_graph[q])
        bts_m = set(multi_match_list_bts[q])
        bfs_m = set(multi_match_list_bfs[q])

        inter = set.intersection(faiss_m, graph_m, bts_m, bfs_m)
        temp = [*inter]
        temp.sort()
        plh = ''
        
        for t in temp:
            plh+=t

        ind = faiss_df[faiss_df['Query'] == q].index[0]
        faiss_df.loc[ind, "T1"] = plh

        ind = graph_df[graph_df['Query'] == q].index[0]
        graph_df.loc[ind, "T1"] = plh

        ind = bts_df[bts_df['Query'] == q].index[0]
        bts_df.loc[ind, "T1"] = plh

        ind = bfs_df[bfs_df['Query'] == q].index[0]
        bfs_df.loc[ind, "T1"] = plh

    except KeyError:
        print("Query without multiple similar top matches in at least one method detected: ",q)
        top = graph_df[graph_df['Query'] == q]['T1'].values[0]
        ind = faiss_df[faiss_df['Query'] == q].index[0]
        faiss_df.loc[ind, "T1"] = top

        ind = bts_df[bts_df['Query'] == q].index[0]
        bts_df.loc[ind, "T1"] = top

        ind = bfs_df[bfs_df['Query'] == q].index[0]
        bfs_df.loc[ind, "T1"] = top


bts_set = set([i['Query'] + " " + i['T1'] for i in bts_df[bts_df['Score'] <= 1].iloc()])
bfs_set = set([i['Query'] + " " + i['T1'] for i in bfs_df[bfs_df['Score'] <= 1].iloc()])
graph_set = set([i['Query'] + " " + i['T1'] for i in graph_df[graph_df['Score'] <= 1].iloc()])
faiss_set = set([i['Query'] + " " + i['T1'] for i in faiss_df[faiss_df['Score'] <= 1].iloc()])
blast_set = set([i['Query'] + " " + i['Ref'] for i in e_df.iloc()])


plot_set = {
'bts': bts_set,
'bfs': bfs_set,
'graph': graph_set,
'faiss': faiss_set,
'BLAST' : blast_set
}

plot_venn_5_sets(plot_set)
