import os
import re
import matplotlib.pyplot as plt
from venn import venn

def gen_init():
    results_path = "/Users/sakshmenon/Desktop/results"
    os.chdir(results_path)
    files = [file.name for file in os.scandir()]
    
    file_map = {'pc': {'sp': {'bfs': '','bts': '','graph': '', 'faiss': ''}, 'sc': {'bfs': '','bts': '','graph': '', 'faiss': ''}, 'td': {'bfs': '','bts': '','graph': '', 'faiss': ''}, 'td_sp_sc': {'bfs': '','bts': '','graph': '', 'faiss': ''}},
                'pj': {'sp': {'bfs': '','bts': '','graph': '', 'faiss': ''}, 'sc': {'bfs': '','bts': '','graph': '', 'faiss': ''}, 'td': {'bfs': '','bts': '','graph': '', 'faiss': ''}, 'td_sp_sc': {'bfs': '','bts': '','graph': '', 'faiss': ''}},
                'pm': {'sp': {'bfs': '','bts': '','graph': '', 'faiss': ''}, 'sc': {'bfs': '','bts': '','graph': '', 'faiss': ''}, 'td': {'bfs': '','bts': '','graph': '', 'faiss': ''}, 'td_sp_sc': {'bfs': '','bts': '','graph': '', 'faiss': ''}}}

    for file in files:
        segments = file.split('.')
        organisms = segments[0]
        algorithm = segments[1]
        mode = segments[2]

        organisms = organisms.split('2')
        org1 = organisms[0]
        org2 = organisms[1]

        org1 = org1[:3:2]

        if mode == 'k3':
            org2 = org2[0] + org2[org2.find('_') + 1]

        file_map[org1][org2][algorithm] = file

    return file_map

def extract_queries(file_map, org1, org2):
    pattern = r"Query:\s+>([\w.]+)\s.*?\nTop1:\s+[\d.]+\s+>([\w.]+)"
    algo_dict = file_map[org1][org2]
    for algo in algos:
        file = algo_dict[algo]
        with open(file) as obj:
            queries = {}
            file_content = obj.read()
            query_top1_pairs = re.findall(pattern, file_content)
            query_top1_dict = {query: top1 for query, top1 in query_top1_pairs}
            for id in query_top1_dict:
                if id not in queries:
                    queries.update({id: query_top1_dict[id]})
                else:
                    queries[id].append(query_top1_dict[id])
                algo_map[algo].append(query_top1_dict[id])
            algo_query_map[algo] = queries

    return queries, algo_query_map

def plot_venn_4_sets(sets_dict, org1, org2):
    venn(sets_dict)
    plt.title(str(org1) + ' - ' + str(org2) + " Venn Diagram")
    plt.show()

org2s = ['sc', 'sp', 'td', 'td_sp_sc']
algos = ['bfs', 'bts', 'graph', 'faiss']

file_map = gen_init()

algo_map = {'bfs': [],'bts': [],'graph': [], 'faiss': []}
algo_query_map = {'bfs': [],'bts': [],'graph': [], 'faiss': []}

org1 = 'pm'
org2 = 'sc'
queries, aqm = extract_queries(file_map, org1, org2)

plot_set = {
'bts': {str(key) + ' ' + str(value) for key, value in aqm['bts'].items()},
'bfs': {str(key) + ' ' + str(value) for key, value in aqm['bfs'].items()},
'graph': {str(key) + ' ' + str(value) for key, value in aqm['graph'].items()},
'faiss': {str(key) + ' ' + str(value) for key, value in aqm['faiss'].items()}
}

plot_venn_4_sets(plot_set, org1, org2)
