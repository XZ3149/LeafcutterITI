import numpy as np
import pandas as pd
import sys
import networkx as nx
import re






def gen_exon_set(exon_count, threshold = 10):
    
    """
    This function return a dic that give out exon that exist
    """
    
    df = pd.read_csv(exon_count, sep = ' ', index_col = 0)


    if "Gene" in df.columns.tolist(): # deal with the result from different version
        df['sum'] = df.iloc[:, 4:].sum(axis=1)
        
    else:
        df['sum'] = df.iloc[:, 3:].sum(axis=1)
    # eliminate exon that have 
    df = df[df['sum'] < threshold]
    result = set(df.index)

    return result





def gen_intron_exon_dic(intron_exon_connectivity, exon_set):
    
    """
    This function generate a dic that provide the possible exon connection for introns
    """
    df = pd.read_csv(intron_exon_connectivity, sep = ' ')
    df['intron'] = df['intron'].str.replace('-', ':', regex=False)
    df['near_exons'] = df['near_exons'].str.replace('-', ':', regex=False)
    result = {intron: set(exons.split(',')) for intron, exons in zip(df['intron'], df['near_exons'])}
    result = {intron: exons.intersection(exon_set) for intron, exons in result.items()} # this will give out only valid exons
    


    return result






def build_intron(row, intron_to_exon_dic):
    """
    a helper function for process_cluster function
    This function is building intron in format: {intron:}

    Parameters
    ----------
    
    
    index: the index of the df row, should be the name of intron
    row : a df row in format (Cluster Chr Start End Gene samples)
    samples: a list of sample names
    exon_count : df
    intron_to_exon : df

    Returns
    -------
    a dict represent the intron

    """
    intron = dict()
    
    intron['chr'] = row['chr']
    intron['start'] = row['start']
    intron['end'] = row['end']
    intron['splice_sites'] = set([row['start'], row['end']])
    intron['name'] = row['name']
    intron['value'] = row['sum'] # reserve for possible update

    exons = intron_to_exon_dic[row['name']]
    intron['exons'] = exons
    
    return intron




def overlaps(A,B):
    '''
    Checks if intron a overlap with intron b 
    '''

    if A['end'] < B['start'] or B['end'] < A['start']:
        return False
    else: return True


def check_shared_splice_site(A,B):
    '''
    Checks if intron a have shared splice site with intron b 
    '''
    
    return A['splice_sites'].intersection(B['splice_sites'])


def check_overlap_exons(A,B):
    return A['exons'].intersection(B['exons'])



def process_cluster(df, intron_exon_dic, cluster_def):
    """
    This is the function process each clusters and return a dict
    cluster_def: 1: overlap, 2: overlap+share_intron_splice_site,
        3: overlap+share_intron_splice_site+shared_exon_splice_site
    """
    
    intron_list = df.apply(lambda row: build_intron(row, intron_exon_dic), axis=1).tolist()
    
    
    G = nx.Graph() # generate a empty graph

    for i in range(len(intron_list)):
        for j in range(i,len(intron_list)):
            intron1 = intron_list[i]
            intron2 = intron_list[j]
            # weight = intron1['sum'] + intron2['sum'] # update to add weight for edges
            if cluster_def == 1:
                if overlaps(intron1, intron2):
                    G.add_edge(intron1['name'], intron2['name'])
            elif cluster_def == 2:
                if check_shared_splice_site(intron1, intron2):
                    G.add_edge(intron1['name'], intron2['name'])
    
            elif cluster_def == 3:
                if check_shared_splice_site(intron1, intron2) or check_overlap_exons(intron1, intron2):
                    G.add_edge(intron1['name'], intron2['name'])
            else:
                sys.exit("Error: invalid cluster definition...\n")


    edges = list(G.edges())
    edges = [tuple(sorted(e)) for e in edges] # sort the edges to ensure a certain order
    return edges

def extract_number(cluster):
    match = re.search(r'clu_(\d+)_', cluster)
    if match:
        return int(match.group(1))
    return None



def gen_cluster_edges(cluster_file, exon_count, intron_exon_connectivity, exon_treshold = 10, cluster_def = 3, out_prefix = ''):
    """
    This is a function that return a df that contain the edges for clusters

    """

    df = pd.read_csv(cluster_file,sep = ' ', index_col = 0)
    df['sum'] = df.sum(axis=1)
    df[['chr', 'start', 'end', 'cluster']] = df.index.to_series().str.split(':', expand=True)
    # df['cluster'] = df['cluster'].str[:-2] # get rid of the strand sign
    df['name'] = df.index.to_series().apply(lambda x: ':'.join(x.split(':')[:-1]))
    
    

    valid_exons = gen_exon_set(exon_count, exon_treshold)
    intron_exon_dic = gen_intron_exon_dic(intron_exon_connectivity, valid_exons)

    columns = ['cluster_name', 'chr', 'start', 'end', 'cluster_edges'] # possible update to add weight for each edge
    result_df = pd.DataFrame(columns=columns)

    clusters = df.groupby('cluster')
    # Iterate over each group
    for cluster_name, cluster_df in clusters:
        cluster_chr = cluster_df['chr'].unique()[0]
        cluster_start = cluster_df['start'].min()
        cluster_end = cluster_df['end'].max()
        cluster_edges = process_cluster(cluster_df, intron_exon_dic, cluster_def)
        
        cluster_edges = ['-'.join(t) for t in cluster_edges]
        cluster_edges = ",".join(cluster_edges) # a string in format intron1-intron2,intron3-intron4......


        result_df.loc[len(result_df)] = [cluster_name, cluster_chr, cluster_start, cluster_end, cluster_edges]
        
    result_df['number'] = result_df['cluster_name'].apply(extract_number)
    result_df = result_df.sort_values(by='number')

    # Drop the helper column if not needed
    result_df= result_df.drop(columns='number')
    result_df = result_df.reset_index(drop=True)
    result_df.to_csv(f'{out_prefix}clusters_edges.tsv', sep = '\t', index = False)












