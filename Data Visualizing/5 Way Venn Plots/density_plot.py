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


def plot_multiple_score_distributions(cutoff, conf_flag, *data_lists):
    # Define cutoff values from 0.1 to 1.5 with increments of 0.1
    cutoff_values = np.arange(0.1, cutoff, 0.1)
    labels = ['BTS', 'BFS', 'GRAPH', 'FAISS']
    
    # Colors for each dataset
    colors = ["skyblue", "salmon", "lightgreen", "gold"]
    
    # Plot histogram for each dataset
    plt.figure(figsize=(10, 6))
    for i, data in enumerate(data_lists):
        # Calculate the number of items below each cutoff for the current list
        if conf_flag:
            counts = [np.sum(np.array(data) >= cutoff) for cutoff in cutoff_values]
        else:
            counts = [np.sum(np.array(data) <= cutoff) for cutoff in cutoff_values]
        
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

def plot_exponential_cutoff_distribution(scores):
    # Define exponential cutoff values from 0.00001 to 1
    cutoff_values = np.logspace(-5, 0, num=6)

    # Calculate the number of items below each cutoff
    counts = [np.sum(np.array(scores) < cutoff) for cutoff in cutoff_values]

    # Plot the histogram
    plt.figure(figsize=(10, 6))
    plt.bar(range(len(cutoff_values)), counts, width=0.2, color='skyblue', edgecolor='black')

    # Customize plot
    plt.xticks(range(len(cutoff_values)), [f"{val:.5f}" for val in cutoff_values], rotation=45)
    plt.xlabel("Cutoff Value")
    plt.ylabel("Number of Scores Below Cutoff")
    plt.title("Distribution of Scores Below Exponential Cutoff Values")

    # Annotate each bar with the count value
    for i, count in enumerate(counts):
        plt.text(i, count + 0.5, str(count), ha='center', va='bottom', fontsize=8)

    plt.tight_layout()
    plt.show()


with open("baker2fission.e1.txt") as fileobj:  #file path here
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

plot_multiple_score_distributions(1.5, 0, bts_df['Score'], bfs_df['Score'], graph_df['Score'], faiss_df['Score'])
plot_multiple_score_distributions(1.1, 1, bts_df['Conf'], bfs_df['Conf'], graph_df['Conf'], faiss_df['Conf'])
plot_exponential_cutoff_distribution(e_df['Score'])