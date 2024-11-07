import json
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

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


def plot_multiple_score_distributions(cutoff, *data_lists):
    # Define cutoff values from 0.1 to 1.5 with increments of 0.1
    cutoff_values = np.arange(0.1, cutoff, 0.1)
    labels = ['BTS', 'BFS', 'GRAPH', 'FAISS']
    
    # Colors for each dataset
    colors = ["skyblue", "salmon", "lightgreen", "gold"]
    
    # Plot histogram for each dataset
    plt.figure(figsize=(10, 6))
    for i, data in enumerate(data_lists):
        # Calculate the number of items below each cutoff for the current list
        counts = [np.sum(np.array(data) < cutoff) for cutoff in cutoff_values]
        
        # Plot the histogram for this list
        plt.bar(cutoff_values + i * 0.02, counts, width=0.02, color=colors[i % len(colors)],
                edgecolor="black", align="center", label=labels[i])
        
    # Set plot labels and title
    plt.xlabel("Cutoff Score")
    plt.ylabel("Number of Items Below Cutoff")
    plt.title("Number of Items Below Cutoff Values for Multiple Datasets")
    plt.xticks(cutoff_values)
    plt.legend()
    plt.show()

faiss_json = "" #json file path here
graph_json = "" #json file path here
bts_json = "" #json file path here
bfs_json = "" #json file path here

faiss_df = f_ext(faiss_json)
graph_df = f_ext(graph_json)
bts_df = f_ext(bts_json)
bfs_df = f_ext(bfs_json)

plot_multiple_score_distributions(1.5, bts_df['Score'], bfs_df['Score'], graph_df['Score'], faiss_df['Score'])
plot_multiple_score_distributions(1.1, bts_df['Conf'], bfs_df['Conf'], graph_df['Conf'], faiss_df['Conf'])