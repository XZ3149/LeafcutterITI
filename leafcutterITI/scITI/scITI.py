import pyranges as pr
import numpy as np
import pandas as pd
import sys
import warnings
import leafcutterITI.utils
from leafcutterITI.utils import timing_decorator,write_options_to_file
import scipy
import scipy.sparse
from scipy.sparse import csr_matrix, save_npz, load_npz
from optparse import OptionParser
from utils import timing_decorator
from utils import write_options_to_file
warnings.simplefilter(action='ignore', category=pd.errors.DtypeWarning)


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


    # Convert to CSR format if needed
    matrix_csr = isoform_intron_matrix.tocsr()
    return isoform_intron_matrix.tocsr()




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

    df_intron.to_csv(f'{out_prefix}count_intron')
    df_exon.to_csv(f'{out_prefix}count_exon')    






































