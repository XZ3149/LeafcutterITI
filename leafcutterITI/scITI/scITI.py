import pyranges as pr
import numpy as np
import pandas as pd
import sys
import warnings
import leafcutterITI.utils
from leafcutterITI.utils import timing_decorator,write_options_to_file
import scipy
import scipy.sparse
from scipy.sparse import csr_matrix, save_npz, load_npz, vstack
from optparse import OptionParser
from utils import timing_decorator
from utils import write_options_to_file
import random
from pathlib import Path
warnings.simplefilter(action='ignore', category=pd.errors.DtypeWarning)
from shared_functions import build_init_cluster, process_clusters, compute_ratio







def barcode_group_print(barcode_group_dic, out_prefix = ''):
    """
    print a file with barcodes and its corresponding group
    """
    with open(f'{out_prefix}barcodes_to_pseudobulk.txt', 'w') as output:
        for group in barcode_group_dic:
            for barcode in barcode_group_dic[group]:
                print(f'{barcode},{group}', file = output)



def metacell_generation(type_dic, n, out_prefix = ''):
    """
    this is the metacell generate function
    type_dic: a dict that contain type/cluster:[barcodes]
    n: number of cell merge to a metacell
    out_prefix: the prefix of the output file
    """
    
    metacell_dic = {}
    with open(f'{out_prefix}meta_cells_record.txt', 'w') as output, open(f'{out_prefix}meta_group.txt', 'w') as group_file:
        print('type/cluster num', file=output)
        for key in type_dic:
            barcodes = type_dic[key]
            random.shuffle(barcodes) # randomize the barcodes
            num_meta_cell = len(barcodes)//n  #integer division
            # if avalible cell < n, then 0 metacell will be generated
    
            print(f'{key} {num_meta_cell}', file=output)
            for i in range(num_meta_cell):
                print(f'{key}_{i} {key}', file = group_file)
                metacell_dic[f'{key}_{i}'] = barcodes[n*i: n*(i+1)]
    
    return metacell_dic
    

def bootstrapping(type_dic, n, k, out_prefix = ''):
    """
    this is the bootstrapping generate function
    type_dic: a dict that contain type/cluster:[barcodes]
    n: number of cell to merge
    k: round of bootstrapping for each cell type or cluster
    out_prefix: the prefix of the output file
    """
    bootstrap_dic = {}
    with open(f'{out_prefix}bootstrapping_record.txt', 'w') as output, open(f'{out_prefix}bootstrapping_group.txt', 'w') as group_file:
        print('type/cluster num', file=output)
        for key in type_dic:
            barcodes = type_dic[key]
            
            
            if len(barcodes) < n: # if avalible cell < n, then skip this cell type
                print(f'{key} 0', file=output)
                
            else:
                
                for i in range(k):
                    print(f'{key}_{i} {key}', file = group_file)
                    bootstrap_dic[f'{key}_{i}'] = random.choices(barcodes , k=k)

    return bootstrap_dic
    
    
    



def pseudo_group_generation(barcodes_type_file, n, k = 30, group_method = 'metacells', out_prefix = ''):
    '''
    barcodes_type_file: a csv file that contain two columns barcodes cell_type/cluster
    group_method: can be metacells or bootstrapping, otherwise retrun error message
    n: number of barcodes to be merged into a single pseudobulk sample
    k: number of bootstrapping round when generating pseudobulk samples, unused when group_method == metacells
    Returns
    -------
    None.

    '''
    df = pd.read_csv(barcodes_type_file, index_col = 0, header = None, names = ['type/cluster'])
    df = df.sort_values(by='type/cluster')
    # ensure cluster number will be interpreted as str
    df['type/cluster'] = df['type/cluster'].astype(str) 
    # replace the possible space in the cluster or cell type name. Avoid problem in Leafcutter
    # some downstream analysis in leafcutter use space to sepearte values
    df['type/cluster'] = df['type/cluster'].str.replace(" ", "_") 
    

    types = df['type/cluster'].unique()
    type_dic = {}
    for name in types:
        type_dic[name] = list(df[df['type/cluster'] == name].index)




    if group_method == 'metacells':
        merged_dic = metacell_generation(type_dic, n)
        barcode_group_print(merged_dic, f'{out_prefix}meta_')
    else: 
        merged_dic = bootstrapping(type_dic, n, k)
        barcode_group_print(merged_dic, f'{out_prefix}bootstrapping_')














@timing_decorator
def isoform_intron_sparse_generation(isoform_intron_map_file, out_prefix = ''):
    """
    

    Returns
    -------
    None.

    """


    df = pd.read_csv(isoform_intron_map_file, sep = ' ')


    isoform_to_introns = dict(zip(df['Transcript'], df['support_introns']))
    isoform_to_exons = dict(zip(df['Transcript'], df['support_exons']))



    for isoform, introns in isoform_to_introns.items():
    # Check if introns is not NaN (using pandas isna() function)
        if pd.isna(introns):
            isoform_to_introns[isoform] = []
        else:
            isoform_to_introns[isoform] = introns.split(',')


    for isoform, exons in isoform_to_exons.items():
    # Check if introns is not NaN (using pandas isna() function)
        if pd.isna(exons):
            isoform_to_exons[isoform] = []
        else:
            isoform_to_exons[isoform] = exons.split(',')



    # Step 1: Extract unique isoforms and introns
    all_isoforms = list(isoform_to_introns.keys())
    
    all_introns = set()
    for introns in isoform_to_introns.values():
        all_introns.update(introns)
    all_introns = list(all_introns)    
    
    all_exons = set()
    for exons in isoform_to_exons.values():
        all_exons.update(exons)
    all_exons = list(all_exons)



    # Step 2: Create mappings to indices
    isoform_to_index = {isoform: i for i, isoform in enumerate(all_isoforms)}
    intron_to_index = {intron: i for i, intron in enumerate(all_introns)}
    exon_to_index = {exon: i for i, exon in enumerate(all_exons)}




    # Step 3: Populate the matrix
    row_indices = []
    col_indices = []
    for isoform, introns in isoform_to_introns.items():
        row_index = isoform_to_index[isoform]
        for intron in introns:
            col_index = intron_to_index[intron]
            row_indices.append(row_index)
            col_indices.append(col_index)

    num_rows = len(all_isoforms)
    num_intron_cols = len(all_introns)


    isoform_intron_matrix = scipy.sparse.coo_matrix(
        ( [1]*len(row_indices), (row_indices, col_indices) ),
        shape=(num_rows, num_intron_cols)
        )


    row_indices = []
    col_indices = []
    for isoform, exons in isoform_to_exons.items():
        row_index = isoform_to_index[isoform]
        for exon in exons:
            col_index = exon_to_index[exon]
            row_indices.append(row_index)
            col_indices.append(col_index)

    num_exon_cols = len(all_exons)
    
    isoform_exon_matrix = scipy.sparse.coo_matrix(
        ( [1]*len(row_indices), (row_indices, col_indices) ),
        shape=(num_rows, num_exon_cols)
        )




    save_npz(f'{out_prefix}isoform_intron_matrix.npz', isoform_intron_matrix)
    save_npz(f'{out_prefix}isoform_exon_matrix.npz', isoform_exon_matrix)


    with open(f'{out_prefix}isoform_rows.txt', 'w') as output:
        for isoform in all_isoforms:
            print(isoform, file= output)

    with open(f'{out_prefix}intron_cols.txt', 'w') as output:
        for intron in all_introns:
            print(intron, file= output)

    with open(f'{out_prefix}exon_cols.txt', 'w') as output:
        for exon in all_exons:
            print(exon, file= output)






def sample_dic_generation(barcode_pseudo):
    """
    This function return a dict that map pseudo sample to a list of barcodes in that sample

    """
    df = pd.read_csv(barcode_pseudo, header = None, sep = ',', index_col =0,  names = ['sample'])
    samples = df['sample'].unique()
    sample_dic = {}
    for name in samples:
        sample_dic[name] = list(df[df['sample'] == name].index)
    return sample_dic
def sample_index_dic_generation(sample_dic, barcodes):
    """
    This function generate a index dict that all barcodes in the experiment
    then find the indexes of barcodes in a sample among all barcodes

    """
    sample_index_dic = {}
    barcodes_index_dic = {}
    for index, barcode in enumerate(barcodes):
        barcodes_index_dic[barcode] = index
    
    for key in sample_dic:
        sample_index_dic[key] = [barcodes_index_dic.get(barcode) for barcode in sample_dic[key]]
    return sample_index_dic

def process_batch(sparse_matrix, batch_row_groups):
    """
    This function used to process sample groups in batch
    it take a sparse_matrix for all barcodes, and using indexes for each sample
    it generate a pseudobulk row for each sample and stack them to a matrix
    """
    
    
    result_matrix =csr_matrix((0, sparse_matrix.shape[1]))
    for group_indices in batch_row_groups:
        # Sum rows efficiently and ensure the result is a 2-D matrix
        summed_row = sparse_matrix[group_indices].sum(axis=0)
        if summed_row.ndim == 1:
            summed_row = summed_row.reshape(1, -1)  # Reshape to 2-D if necessary
        result_matrix = vstack([result_matrix, summed_row]) 

    return result_matrix






def pseudo_eq_conversion(salmon_dir, ):
    """
    

    Parameters
    ----------
    salmon_dir : TYPE
        DESCRIPTION.
     : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    """















def sc_intron_count(reference_directory, pseudo_matrix_file, pseudo_cols_file, pseudo_rows_file, in_prefix = '' , out_prefix = '', threshold = 10):
    """
    

    Parameters
    ----------
    reference_directory : Str
        The directory that contain the reference
    pseudo_matrix_file : Str
        
    pseudo_cols_file : TYPE
        DESCRIPTION.
    pseudo_rows_file : TYPE
        DESCRIPTION.
    in_prefix : TYPE, optional
        DESCRIPTION. The default is ''.
    out_prefix : TYPE, optional
        DESCRIPTION. The default is ''.

    Returns
    -------
    None.

    """
    
    #loading reference matrix
    intron_isoform_matrix = load_npz(f'{reference_directory}/{in_prefix}isoform_intron_matrix.npz').tocsr()
    exon_isoform_matrix = load_npz(f'{reference_directory}/{in_prefix}isoform_exon_matrix.npz').tocsr()
    
    isoform_rows = []
    with open(f'{reference_directory}/{in_prefix}isoform_rows.txt', 'r') as input_file:
        for line in input_file:
            isoform_rows.append(line.split('.')[0]) # get rid of the version number
            
            
    intron_cols = []
    with open(f'{reference_directory}/{in_prefix}intron_cols.txt', 'r') as input_file:
        for line in input_file:
            intron_cols.append(line[:-1]) # get rid of the \n


    exon_cols = []
    with open(f'{reference_directory}/{in_prefix}exon_cols.txt', 'r') as input_file:
        for line in input_file:
            exon_cols.append(line[:-1]) 



    pseudo_matrix = load_npz(pseudo_matrix_file)
    
    
    pseudo_cols = [] #isoforms
    with open(pseudo_cols_file, 'r') as input_file:
        for line in input_file:
            pseudo_cols.append(line.split('.')[0]) 

    pseudo_rows = [] #barcodes
    with open(pseudo_rows_file, 'r') as input_file:
        for line in input_file:
            pseudo_rows.append(line[:-1]) 





    common_isoforms = set(isoform_rows).intersection(pseudo_cols) # get the common isoforoms
    print(len(common_isoforms))

    # get the index for the isoforms 
    reference_index = {isoform: idx for idx, isoform in enumerate(isoform_rows)}




    # Step 3: Rearrange Matrix 1
    # Filter and reorder columns to align with Matrix 2
    filtered_indices = [i for i, isoform in enumerate(pseudo_cols) if isoform in common_isoforms]
    filtered_psudo_cols = [isoform for i, isoform in enumerate(pseudo_cols) if isoform in common_isoforms]
    filtered_pseudo_matrix = pseudo_matrix[:, filtered_indices]


    # Create a new empty matrix with the same number of rows and the correct number of columns
    new_pseudo_matrix = scipy.sparse.lil_matrix((pseudo_matrix.shape[0], len(reference_index)))


    for idx, isoform in enumerate(filtered_psudo_cols):
        if isoform in reference_index:
            col_new_index = reference_index[isoform]
            new_pseudo_matrix[:, col_new_index] = filtered_pseudo_matrix[:, idx]



    intron_count_matrix = new_pseudo_matrix @ intron_isoform_matrix 
    exon_count_matrix = new_pseudo_matrix @ exon_isoform_matrix 



    filtered_columns = intron_count_matrix.sum(axis=0) > threshold  # Filter out lowerly expressed intron
    # Step 2: Filter the matrix to keep only non-empty columns
    filtered_intron_count_matrix = intron_count_matrix[:, filtered_columns.A.ravel()]
    # Step 3: Filter the column names list
    filtered_intron_cols = [name for name, keep in zip(intron_cols, filtered_columns.A.ravel()) if keep]



    filtered_columns = exon_count_matrix.sum(axis=0) > threshold  # Filter out lowerly expressed exon
    # Step 2: Filter the matrix to keep only non-empty columns
    filtered_exon_count_matrix = exon_count_matrix[:, filtered_columns.A.ravel()]
    # Step 3: Filter the column names list
    filtered_exon_cols = [name for name, keep in zip(exon_cols, filtered_columns.A.ravel()) if keep]



    dense_intron = filtered_intron_count_matrix.toarray()
    df_intron = pd.DataFrame(dense_intron).T
    df_intron.index = filtered_intron_cols
    df_intron.columns = pseudo_rows
    
    

    dense_exon = filtered_exon_count_matrix.toarray()
    df_exon = pd.DataFrame(dense_exon).T
    df_exon.index = filtered_exon_cols
    df_exon.columns = pseudo_rows




    df_intron[['Chr', 'Start', 'End']] = df_intron.index.to_series().str.extract(r'(chr\w+):(\d+)-(\d+)')
    new_order = list(df_intron.columns)[-3:] + list(df_intron.columns)[:-3] 
    df_intron = df_intron[new_order]
    
    
    df_exon[['Chr', 'Start', 'End']] = df_exon.index.to_series().str.extract(r'(chr\w+):(\d+)-(\d+)')
    
    new_order = list(df_exon.columns)[-3:] + list(df_exon.columns)[:-3] 
    df_exon = df_exon[new_order]
    

    df_intron.sort_values(by=['Chr', 'Start', 'End'], inplace=True)
    df_exon.sort_values(by=['Chr', 'Start', 'End'], inplace=True)

    df_intron.to_csv(f'{out_prefix}count_intron', sep = ' ')
    df_exon.to_csv(f'{out_prefix}count_exon', sep = ' ')    





"""

if __name__ == "__main__":

    from optparse import OptionParser

    parser = OptionParser()


    parser.add_option("--map_directory", dest="map_directory", default = None,
                  help="the directory of isoform to intron and exon matrices (default: None)")
    
    parser.add_option("--input_prefix", dest="input_prefix", default = '',
                  help="isoform to intron map file (default: None)")
    
    
    parser.add_option("--alevin-fry_dir", dest="alevin-fry_directory", default = None,
                  help="The directory for alevin-fry results, should contain the eq matrix and other files")
    
    


    parser.add_option("--salmon_dir", dest="salmon_dic", default = None,
                  help="the directory for salmon reference, the pipeline will search for spliceu or splicei")





    parser.add_option("--map", dest="map", default = None,
                  help="isoform to intron map file (default: None)")
    
    
    parser.add_option("--count_files", dest="count_files", default = None,
                  help="the input samples files, should be placed in a txt file (default: None)")
    
    parser.add_option("--connect_file", dest="connect_file", default = None,
                  help="Intron and exon connectivity file (default:None)")
    
    parser.add_option("-a","--annot", dest="annot", default = None,
                  help="transcriptome annotation, required if normalization == True (default: None)")

    
    parser.add_option("--cluster_def", dest="cluster_def", default = 3,
                  help="three def available, 1: overlap, 2: overlap+share_intron_splice_site, \
                      3: overlap+share_intron_splice_site+shared_exon_splice_site")
    
    parser.add_option("-o", "--outprefix", dest="outprefix", default = 'leafcutterITI_',
                  help="output prefix , should include the diretory address if not\
                  in the same dic (default leafcutterITI_)")    
    
                  

    parser.add_option("-n", "--normalization", dest="normalization", default = True,
                  help="whether to performance normalization, if not use TPM directly (default: True)")
    
    parser.add_option("--preprocessed", dest="preprocessed", default = False,
                  help="Whether normalization already been done for the count_files (default: False)")
    
        
    parser.add_option("--samplecutoff", dest="samplecutoff", default = 0,
                  help="minimum count for an intron in a sample to count as exist(default: 0)")

    parser.add_option("--introncutoff", dest="introncutoff", default = 5,
                  help="minimum count for an intron to count as exist(default 5)")
    
    parser.add_option("-m", "--minclucounts", dest="minclucounts", default = 30,
                  help="minimum count in a cluster (default 30 normalized count)")
    
    parser.add_option("-r", "--mincluratio", dest="mincluratio", default = 0.01,
                  help="minimum fraction of reads in a cluster that support a junction (default 0.01)")


    
    (options, args) = parser.parse_args()



    if options.count_files is None:
        sys.exit("Error: no count_files files provided...\n")
        
    if options.map is None:
        sys.exit("Error: no isoform intron map provided...\n")
        
    if options.connect_file is None:
        sys.exit("Error: Intron exon connectivity file provided...\n")
    

    if options.normalization == True and options.annot == None:
        sys.exit("Error: no annotation file provided, this is required if normalization == True...\n")
"""      

   



































