import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde
from venn import venn
import pandas as pd
import json

def plot_density_curves(e_fl, *lists):
    # Prepare the data lists
    data_lists = lists
    colors = ['blue', 'green', 'red', 'purple']
    labels = ['BTS', 'BFS', 'GRAPH', 'FAISS']
    
    # Set up the plot
    fig, axes = plt.subplots(1, 1, figsize=(12, 6))
    
    # Density Plot
    for i, data in enumerate(data_lists):
        density = gaussian_kde(data)
        x_vals = np.linspace(min(data) - 0.5, max(data) + 0.5, 1000)
        y_vals = density(x_vals)
        
        # Plot density curves
    
    # Histogram with Density Curves Overlay
    for i, data in enumerate(data_lists):
        # Plot histogram for each list with transparency
        if e_fl:
            axes.hist(data, bins=30, density=True, color=colors[i], alpha=0.3, label='E Value Histogram')

        else:
            axes.hist(data, bins=30, density=True, color=colors[i], alpha=0.3, label=f'{labels[i]} Histogram')
        
        # Overlay density curve on histogram
        density = gaussian_kde(data)
        y_vals = density(x_vals)
        if e_fl:
            axes.plot(x_vals, y_vals, color=colors[i], label='Density')

        else:
            axes.plot(x_vals, y_vals, color=colors[i], label=f'{labels[i]} Density')
    
    axes.set_title('Histogram with Density Curves')
    axes.set_xlabel('Data Values')
    axes.set_ylabel('Density')
    axes.legend()
    axes.grid(True)
    
    # Adjust layout and display the plots
    plt.tight_layout()
    plt.show()

def plot_venn_5_sets(sets_dict):
    venn(sets_dict)
    plt.title("4 Algo Venn Diagram")
    plt.show()

def f_ext(file):
    with open(file, 'r') as f:
        obj = json.load(f)
    main_l = []
    for q in obj:
        que = q['Query']['ID'].split(' ')[0][1:]
        t1 = q['Query']['Matches'][0]['ID'].split(' ')[0][1:]
        conf = q['Query']['Confidence']
        score = q['Query']['Matches'][0]['Score']
        main_l.append({'Query' : que, 'T1' : t1, 'Score' : score,'Conf' : conf})
    return pd.DataFrame(main_l)

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

faiss_json = "faiss_top10.json" #file path here
graph_json = "graph_top10.json" #file path here
bts_json = "bts_top10.json" #file path here
bfs_json = "bfs_top10.json" #file path here

faiss_df = f_ext(faiss_json)
graph_df = f_ext(graph_json)
bts_df = f_ext(bts_json)
bfs_df = f_ext(bfs_json)

bts_set = set([i['Query'] + i['T1'] for i in bts_df[bts_df['Score'] <= 1].iloc()])
bfs_set = set([i['Query'] + i['T1'] for i in bfs_df[bfs_df['Score'] <= 1].iloc()])
graph_set = set([i['Query'] + i['T1'] for i in graph_df[graph_df['Score'] <= 1].iloc()])
faiss_set = set([i['Query'] + i['T1'] for i in faiss_df[faiss_df['Score'] <= 1].iloc()])

plot_set = {
'bts': bts_set,
'bfs': bts_set,
'graph': graph_set,
'faiss': faiss_set,
}

import math
ref_set = set(bts_df['Query'][bts_df['Score'] <= 1].to_list() + bfs_df['Query'][bfs_df['Score'] <= 1].to_list() + graph_df['Query'][graph_df['Score'] <= 1].to_list() + faiss_df['Query'][faiss_df['Score'] <= 1].to_list())
e_val = []
raw_e_val = []
for q in ref_set:
    if q in e_df['Query'].to_list():
        try:
            raw_e_val.append(float(e_df[e_df['Query'] == q]['Score']))
            e_val.append(-math.log10(float(e_df[e_df['Query'] == q]['Score'])))
        except:
            pass

plot_venn_5_sets(plot_set)

plot_density_curves(0, bts_df[bts_df['Score'] <= 1]['Score'], bfs_df[bfs_df['Score'] <= 1]['Score'], graph_df[graph_df['Score'] <= 1]['Score'], faiss_df[faiss_df['Score'] <= 1]['Score'])

plot_density_curves(1, e_val)

