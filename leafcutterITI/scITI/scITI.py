import pyranges as pr
import numpy as np
import pandas as pd
import sys
import warnings
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
import calcutta_functions
import utils
from utils import timing_decorator,write_options_to_file


"""

import leafcutterITI.utils
from leafcutterITI.utils import timing_decorator,write_options_to_file
import calcutta_functions
from leafcutterITI.shared_functions import build_init_cluster, process_clusters, compute_ratio

"""

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
    n = int(n)
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
    n = int(n)
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






def pseudo_dic_generation(barcode_pseudo_df):
    """
    This function return a dict that map pseudo sample to a list of barcodes in that sample

    """
    
    samples = barcode_pseudo_df['sample'].unique()
    sample_dic = {}
    for name in samples:
        sample_dic[name] = list(barcode_pseudo_df[barcode_pseudo_df['sample'] == name].index)
    return sample_dic

def pseudo_index_dic_generation(sample_dic, barcodes):
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


def check_barcodes_exsitent(barcode_pseudo_df, valid_barcodes, out_prefix):
    '''
    This function check whether the barcodes in barcode_pseudo_df is in the valid barcodes list 
    Will write a file that contain the final number of barcodes for each pseudo sample and a list 
    of eliminated samples
    will return a new df that only contain the valid barcodes

    '''

    filtered_index_list = [x for x in barcode_pseudo_df.index if x in valid_barcodes]
    filter_out_index_list = [x for x in barcode_pseudo_df.index if x not in valid_barcodes]
    
    
    filtered_df = barcode_pseudo_df.loc[filtered_index_list]


    with open(f'{out_prefix}eliminated_barcodes.txt', 'w') as file:
        for barcode in filter_out_index_list:
            file.write(f"{barcode}\n")
            
            
    pseudo_count = filtered_df['sample'].value_counts().sort_index()

    with open(f'{out_prefix}pseudo_barcodes_counts.txt', 'w') as file:
        for pseudo, count in pseudo_count.items():
            file.write(f"{pseudo} {count}\n")


    sys.stderr.write(f"There are {len(filter_out_index_list)} barcodes get eliminated \n")
    
    return  filtered_df



def pseudo_eq_conversion(alevin_dir, salmon_ref, barcode_pseudo_file, min_EC = 5, min_transcript = 0, out_prefix= ''):
    """
    

    Parameters
    ----------
    salmon_dir : the dir for where the salmon barcdeos eq matrix and eq transcript matrix are
    salmon_ref: the reference .fa file that salmon used for aligenment

    Returns
    -------
    None.

    """

    
    # step 1: data loading
    transcript_lengths_dic = calcutta_functions.get_transcript_lengths(Path(salmon_ref))    
    
    
    
    barcodes_pseudo = pd.read_csv(barcode_pseudo_file, header = None, sep = ',', index_col =0,  names = ['sample'])
    alevin_dir = Path(alevin_dir)
    
    map_cache = alevin_dir / "gene_eqclass.npz"
    if map_cache.is_file(): 
        ec_transcript_mat = scipy.sparse.load_npz(map_cache)
    else:
        num_genes, num_ec, ecs = calcutta_functions.read_alevin_ec(alevin_dir/ "gene_eqclass.txt.gz") 
        ecs_list = [ ecs[k] for k in range(len(ecs)) ] # this works because indexing is dense
        ec_transcript_mat = calcutta_functions.to_coo(ecs_list)
        scipy.sparse.save_npz(map_cache, ec_transcript_mat)


    # Load cells x EC counts
    npz_cache = alevin_dir  / "geqc_counts.npz"

    if npz_cache.is_file(): 
        cell_ec_sparse = scipy.sparse.load_npz(npz_cache)

    else:
        cell_ec_sparse = scipy.io.mmread(alevin_dir / "geqc_counts.mtx") 
        scipy.sparse.save_npz(npz_cache, cell_ec_sparse)  
    cell_ec_sparse = cell_ec_sparse.tocsr()  # convert coo to csr format, necessary for row operation
    
    
    with open(alevin_dir / "quants_mat_cols.txt", 'r') as file:
        features = np.array([line.strip() for line in file.readlines()])
    with open(alevin_dir / "quants_mat_rows.txt", 'r') as file:
        barcodes = np.array([line.strip() for line in file.readlines()])
    
    
    # step 2: first round of filtering, to facilitate computation speed
    #aggregate all barcodes to give a single pseudobulk sample
    pseudobulk = calcutta_functions.sparse_sum(cell_ec_sparse,0) 
    
    transcript_count = pseudobulk @ ec_transcript_mat
    ec_sizes = calcutta_functions.sparse_sum(ec_transcript_mat,1)
    
    
    # first round of filtering, this will speed up the pseudobulk matrix generation
    ECs_to_keep = pseudobulk >= min_EC
    ec_transcript_filt = ec_transcript_mat.tocsr()[ECs_to_keep,:]
    cell_ec_filt = cell_ec_sparse[:,ECs_to_keep]
    features_to_keep = transcript_count > min_transcript
    features_filt = features[features_to_keep]
    ec_transcript_ff = ec_transcript_filt[:,features_to_keep]
    
    
    barcodes_pseudo = check_barcodes_exsitent(barcodes_pseudo,barcodes, out_prefix= f'{out_prefix}')
    pseudo_dict = pseudo_dic_generation(barcodes_pseudo)
    # get a dict that contain information about pseudo_sample: [barcode1,barcode2,.......]
    pseudo_index_dic = pseudo_index_dic_generation(pseudo_dict, barcodes)
    
    



    # step 3: generate the pseudo eq matrix by merging barcodes
    batch_size = 200 # process in batch to avoid crash
    row_group_values = list(pseudo_index_dic.values())
    row_names = list(pseudo_index_dic.keys())
    batches = [row_group_values[i:i + batch_size] for i in range(0, len(row_group_values), batch_size)]

    # Process each batch and incrementally build the new sparse matrix
    pseudo_ec_sparse = csr_matrix((0, cell_ec_filt.shape[1]))  # Initialize an empty sparse matrix
    for batch in batches:
        batch_result = process_batch(cell_ec_filt, batch)
        pseudo_ec_sparse = vstack([pseudo_ec_sparse, batch_result])

    pseudo_ec_sparse = pseudo_ec_sparse.tocsr() 
    
    # step 4: additional filtering to further reduce the matrix size
    aggragate_cell_ec = calcutta_functions.sparse_sum(pseudo_ec_sparse,0)
    transcript_count = aggragate_cell_ec @ ec_transcript_ff
    

    ECs_to_keep = aggragate_cell_ec >= min_EC
    ec_transcript_fff = ec_transcript_ff[ECs_to_keep,:]
    cell_ec_sparse_pseudo_filt = pseudo_ec_sparse.tocsr()[:,ECs_to_keep]
    
    features_to_keep = transcript_count > min_transcript
    features_ff = features_filt[features_to_keep]
    ec_transcript_ffff = ec_transcript_fff[:,features_to_keep]


    feature_lengths, w = calcutta_functions.get_feature_weights(features_ff, transcript_lengths_dic)


    # step 5: compute pseudo transcript matrix:
    ec_transcript_input = ec_transcript_ffff.tocoo()
    pseudo_TPM = csr_matrix((0, ec_transcript_input.shape[1]))  # Initialize an empty sparse matrix
    pseudo_count = csr_matrix((0, ec_transcript_input.shape[1]))  # Initialize an empty sparse matrix
    for i in range(cell_ec_sparse_pseudo_filt.shape[0]):
    
        temp_sample = cell_ec_sparse_pseudo_filt.getrow(i)
        temp_sample = temp_sample.toarray().ravel() # collapse the row into a vector
        
        total_count = temp_sample.sum()
        temp_sample += 0.0000001 # add 1e-6 to avoid null value
        
        alpha = calcutta_functions.EM(temp_sample, ec_transcript_input, w)
        TPM = alpha * 1000000
        count = alpha * total_count 
        # we don't need the actual count, we just need a scale that represent the support levsd   el of the TPM value
        # TPM * total_count will give a count-like value that retain the TPM ratio of transcripts
    
        """
        count = TPM/w
        count = (count/count.sum()) * total_count
        """
        
        pseudo_TPM = vstack([pseudo_TPM, TPM])
        pseudo_count = vstack([pseudo_count, count])

    # step 6: store matrices and eliminate unspliced isoform as they doesn't excised intron    
    pseudo_names = row_names
    
    spliced_indices = [i for i, s in enumerate(features_ff) if "ENSMUSG" not in s]    
    spliced_feature = [s for i,s in enumerate(features_ff) if "ENSMUSG" not in s]
    unspliced_indices = [i for i, s in enumerate(features_ff) if "ENSMUSG" in s]

    spliced_pseudo_TPM = pseudo_TPM.tocsr()[:,spliced_indices]
    spliced_pseudo_count = pseudo_count.tocsr()[:,spliced_indices]


    with open(f'{out_prefix}pseudo_cols.txt', 'w') as output:
        for trans in features_ff:
            print(trans, file=output)
    with open(f'{out_prefix}pseudo_rows.txt', 'w') as output:
        for pseudo in pseudo_names:
            print(pseudo, file=output)
    with open(f'{out_prefix}pseudo_spliced_cols.txt', 'w') as output:
        for trans in spliced_feature:
            print(trans, file=output)
            
    

    save_npz(f'{out_prefix}pseudo_TPM.npz', pseudo_TPM)
    save_npz(f'{out_prefix}pseudo_count.npz', pseudo_count)
    save_npz(f'{out_prefix}pseudo_spliced_TPM.npz', spliced_pseudo_TPM)
    save_npz(f'{out_prefix}pseudo_spliced_count.npz', spliced_pseudo_count)





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



def LeafcutterITI_scITI(options):
    """
    This is the main function for LeafcutterITI_scITI

    """
    #  step 1: generate or read the barcodes to pseudobulk samples file
    out_prefix = options.outprefix
    if options.pseudobulk_samples == None:
    
        pseudo_group_generation(options.barcodes_cluster, options.num_cell, options.num_bootstrapping, \
                                options.group_method, options.outprefix)
        if options.group_method == 'metacells':
            pseudobulk_name = f'{out_prefix}meta_barcodes_to_pseudobulk.txt'
        else:
            pseudobulk_name = f'{out_prefix}bootstrapping_barcodes_to_pseudobulk.txt'
            
        sys.stderr.write("Barcodes to pseudobulk samples generated\n")

    else:
        pseudobulk_name = options.pseudobulk_samples 
        sys.stderr.write("Barcodes to pseudobulk samples read in\n")
    


    # step 2: compute the pseudobulk isoform matrix 
    
    
    pseudo_eq_conversion(options.alevin_dir, options.salmon_ref, pseudobulk_name, min_EC = options.min_eq, min_transcript = 0, out_prefix= out_prefix)
    
    sys.stderr.write("pseudobulk eq matrix were computed\n")
    
    # step 3: counting intron and exon
    ref_prefix = f'{options.ref_dir}/{options.ref_prefix}'
        

    if options.normalization == True:
        pseudo_matrix_file = f'{out_prefix}pseudo_spliced_count.npz'
    else:
        pseudo_matrix_file = f'{out_prefix}pseudo_spliced_TPM.npz'


    pseudo_cols_file = f'{out_prefix}pseudo_spliced_cols.txt'
    pseudo_rows_file = f'{out_prefix}pseudo_rows.txt'

    sc_intron_count(options.ref_dir, pseudo_matrix_file,pseudo_cols_file, pseudo_rows_file,  in_prefix= options.ref_prefix ,\
                    out_prefix = out_prefix, threshold = options.introncutoff)
    sys.stderr.write("Intron and exon were quantified\n")
    # a csv file that is compatible with the rest of leafcutterITI is obtained
    
    # step 4: clustering, using shared function with leafuctterITI


    build_init_cluster(f'{out_prefix}count_intron')
    
    sys.stderr.write("Finished Initial Clustering\n")


    if options.with_virtual == False:
        connect_file = f'{ref_prefix}intron_exon_connectivity.tsv'
    else:
        connect_file = f'{ref_prefix}intron_exon_connectivity_with_virtual.tsv'

    process_clusters(f'{out_prefix}count_intron', f'{out_prefix}count_exon', \
                     connect_file, \
                     out_prefix = out_prefix,\
                     cutoff = options.introncutoff, \
                     percent_cutoff = options.mincluratio,\
                     min_cluster_val = options.minclucounts)

    
    
    sys.stderr.write("Finished Cluster refinement\n")
    
    compute_ratio(f'{out_prefix}refined_cluster', out_prefix)

    sys.stderr.write("Finished PSI calculation\n")

    record = f'{options.outprefix}clustering_parameters.txt'
    sys.stderr.write(f'saving paramters to {record}\n')
    write_options_to_file(options, record)
    
    sys.stderr.write('CLustering Fnished\n')
    







if __name__ == "__main__":

    from optparse import OptionParser

    parser = OptionParser()


    
    # necessary parameter
    parser.add_option("--alevin_dir", dest="alevin_dir", default = None,
                  help="The directory for alevin results, should contain the eq matrix and other files (default: None)")
    
    
    parser.add_option("--salmon_ref", dest="salmon_ref", default = None,
                  help="The reference used for salmon index, The salmon reference,  maybe spliceu or splicei (default: None)")

    parser.add_option("--ref_dir", dest="ref_dir", default = None,
                  help="leafcutterITI reference direcotry (default: None)")

    parser.add_option("--barcodes_cluster", dest="barcodes_cluster", default = None,
                  help="the file that record which barcodes is belonging to which cluster/cell type in format 'barcode,cluster' \
                  this file will be used to generate pseudobulk samples (default: None)")
                  
    parser.add_option('--pseudobulk_samples', dest="pseudobulk_samples", default = None,
                      help="a txt file with barcodes to pseudobulk sample are expected in format 'barcode pseudobulk_ample', if \
                      this option!= None, then it will overwrite the input to --barcodes_cluster, and use the file in this option \
                      for computation (default: None)")
                      
    parser.add_option("--ref_prefix", dest="ref_prefix", default = '',
               help="the reference prefix that is used to generate isoform to intron map using\
               leacutterITI_map_gen (default: '')")            
               



    # optional parameter


    parser.add_option("-n",'--num_cell', dest="num_cell", default = 100,
                  help="the number of cell/barcode that would like to include in a pseudobulk sample, cluster/cell type that have fewer cell than\
                      this number will not included in the computation (default: 100)")
    
    parser.add_option("-k",'--num_bootstrapping', dest="num_bootstrapping", default = 30,
                  help="the number of bootstrapping samples generated for each cluster/cell type (default: 30)")
    parser.add_option("--min_eq", dest="min_eq", default = 5,
                  help="minimum count for each eq class for it to be included in the EM (default: 5)")


    parser.add_option('--group_method', dest="group_method", default = 'metacells',
                  help="the pseudobulk sample generate method, could be metacell or bootstrapping (default: metacells)")

               
    parser.add_option("--normalization", dest="normalization", default = True,
                  help="whether to performance normalization, if not use TPM directly (default: True)")
    
    parser.add_option("--preprocessed", dest="preprocessed", default = False,
                  help="Whether pseudobulk generation and EM were done, if true, \
                      then the pipeline start from counting intron (default: False)")

    parser.add_option("-v","--with_virtual", dest="with_virtual", default = False,
                  help=" whehter the map contain virtual intron to capture AFE and ALE(default: False)")
    
    
    parser.add_option("--cluster_def", dest="cluster_def", default = 3,
                  help="three def available, 1: overlap, 2: overlap+share_intron_splice_site, \
                      3: overlap+share_intron_splice_site+shared_exon_splice_site")
    
    parser.add_option("-o", "--outprefix", dest="outprefix", default = 'leafcutterITI_',
                  help="output prefix , should include the diretory address if not\
                  in the same dic (default leafcutterITI_)")    

        
    parser.add_option("--samplecutoff", dest="samplecutoff", default = 0,
                  help="minimum count for an intron in a sample to count as exist(default: 0)")

    parser.add_option("--introncutoff", dest="introncutoff", default = 5,
                  help="minimum count for an intron to count as exist(default 5)")
    
    parser.add_option("-m", "--minclucounts", dest="minclucounts", default = 30,
                  help="minimum count in a cluster (default 30 normalized count)")
    
    parser.add_option("-r", "--mincluratio", dest="mincluratio", default = 0.01,
                  help="minimum fraction of reads in a cluster that support a junction (default 0.01)")


    
    (options, args) = parser.parse_args()



    if options.salmon_ref is None:
        sys.exit("Error: no salmon reference is provided...\n")
    
    if options.ref_dir is None:
        sys.exit("Error: no leafcutterITI reference is provided...\n")
    
        
    """
    if map_cache.is_file() == None: 
        ec_transcript_mat = scipy.sparse.load_npz(map_cache)


    if options.count_files is None:
        sys.exit("Error: no count_files files provided...\n")
        
    if options.map is None:
        sys.exit("Error: no isoform intron map provided...\n")
        
    if options.connect_file is None:
        sys.exit("Error: Intron exon connectivity file provided...\n")
    

    if options.normalization == True and options.annot == None:
        sys.exit("Error: no annotation file provided, this is required if normalization == True...\n")      
    """
   

    LeafcutterITI_scITI(options)

































