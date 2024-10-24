import numpy as np
import pandas as pd
import sys
import warnings
from optparse import OptionParser
warnings.simplefilter(action='ignore', category=pd.errors.DtypeWarning)







def add_intron_to_cluster(cluster_file, intron_file, out_name = 'merge_clusters'):
    """
    
    Parameters
    ----------
    cluster_file : TYPE
        DESCRIPTION.
    intron_file : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    """
    
    intron_df = pd.read_csv(intron_file, sep = ' ', index_col = 0)
    
    cluster_df = pd.read_csv(cluster_file, sep = ' ', index_col = 0 )
    

    intron_df.index = intron_df.index.str.replace("-", ":")
    intron_df = intron_df.iloc[:, 4:] # drop un necessary information

    valid_introns = cluster_df.index.str.split(':').str[:3].str.join(':').to_list()
    


    intron_df= intron_df.reindex(valid_introns, fill_value=0)
    
    
    result_df = pd.concat([cluster_df.reset_index(drop=True), intron_df.reset_index(drop=True)], axis=1)
    result_df.index = cluster_df.index


    result_df.save_csv(out_name, sep = ' ')
    
    with open(out_name, 'r') as file:
        lines = file.readlines()
      
    lines[0] = lines[0].lstrip()

    # Write the modified content back to the file
    with open(out_name, 'w') as file:
        file.writelines(lines)





















