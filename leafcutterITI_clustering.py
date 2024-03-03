import pyranges as pr
import numpy as np
import pandas as pd
import sys
import warnings
import time
from functools import wraps
from support_functions import timing_decorator
from support_functions import write_options_to_file
warnings.simplefilter(action='ignore', category=pd.errors.DtypeWarning)


chr_dic = {'chr1': 0,'chr2': 1,'chr3': 2,'chr4': 3,'chr5': 4,'chr6': 5,'chr7': 6,'chr8': 7,'chr9': 8,'chr10': 9,
               'chr11': 10,'chr12': 11,'chr13': 12, 'chr14': 13,'chr15': 14,'chr16': 15,'chr17': 16,'chr18': 17,
               'chr19': 18,'chr20': 19,'chr21': 20, 'chr22': 21, 'chrX':22, 'chrY':23}




def transcript_to_intron_counts(row, trans_int_map, sample_name, introns_dic, exons_dic,normalize_count, unmapped_transcript):
    """
    
    Parameters
    ----------

    row: a df row from quant.sf
    trans_int_map : a dataframe that contain map for transcript to introns, name of transcript 
                    will be index
                    row in format: transcript_id chr gene_id support_introns support_exons
                    
    
    sample: str, sample name

    introns_dic : a dictionary in format {intron:[gene,{sample_name:value}]}
    
    exons_dic : a dictionary in format {exon:[gene,{sample_name:value}]}
    normalize_count: a boolean value, whether use normalize count or TPM
    unmapped_transcript: a dic use to record unmapped_transcript that appear 
    
    
    Returns
    -------
    None.
    The dic will be updated, no return is needed
    """

  
    transcript_name = row['Name']
    if normalize_count:
        value = row['normalized_count']  
    else: 
        value = row['TPM']

  
    if (transcript_name in trans_int_map.index) == False:
        unmapped_transcript.add(transcript_name)
        return None
    
    transcript = trans_int_map.loc[transcript_name]
    chromosome = transcript['Chr']
    gene = transcript['Gene']

    
  
    if transcript['support_introns'] != None:    # if there have no support intron, skip to exon section
        
        introns = transcript['support_introns'].split(',') # list of introns
        for intron in introns:
            if intron in introns_dic:
                if sample_name in introns_dic[intron][1]: #{sample_name:TPM}
                    introns_dic[intron][1][sample_name] += value
                else: 
                    introns_dic[intron][1][sample_name] = value
            else:
                introns_dic[intron] = [gene, {sample_name:value}]

        
    if transcript['support_exons'] != None:   # if there have no support intron, skip to next section
        exons = transcript['support_exons'].split(',') # list of exons

        for exon in exons:
            if exon in exons_dic:
                if sample_name in exons_dic[exon][1]: #{sample_name:TPM}
                    exons_dic[exon][1][sample_name] += value
                else: 
                    exons_dic[exon][1][sample_name] = value
            else:
                exons_dic[exon] = [gene, {sample_name:value}]
            
    return None
    



def dic_to_csv(dic, output_name, samples):
    """

    Parameters
    ----------
    dic : an intron or exon dic in format {intron:[gene,{sample_name:TPM}]}
    
    output_name : name of output 
    
    samples: name of samples 

    Returns
    -------
    None
    wirte to disk

    """

    output = open(output_name, 'w')
    
    out_str = 'Name Chr Start End Gene'
    for i in samples:
        out_str += f' {i}'
    
    
    
    print(out_str, file = output)
    
    for key in dic: # key = intron/exon e.g. chr1:10000-12000
    
        tmp = key.split(':')
        chromosome = tmp[0] 
        start, end = tmp[1].split('-')
        
        gene_name = dic[key][0]
        val_dic = dic[key][1]  # {sample:value} 

        tmp_str = f'{key} {chromosome} {start} {end} {gene_name}'
        
        for sample in samples:
            
            
            if sample in val_dic:
                tmp_str += f' {val_dic[sample]}'
            else:
                tmp_str += ' 0'
        
        print(tmp_str, file = output)
                

    output.close()







@timing_decorator
def count_introns(samples, trans_int_map_file, out_prefix = '', threshold = 0, normalize_count = False):
    """
    Parameters
    ----------
    samples : list, names of sample files
          
    trans_int_map_file : str, location of transcript_to_intron map file
        
    out_prefix : str, optional, prefix of output file
    
    threshold: the threshold use for filtering out introns, should decide based on TPM or normalize_count
        
    example: count_introns（['liver_quant.sf', 'heart_quant.sf'],'hg38_transcript_intron_map.tsv' ）
    Returns
    -------
    df of introns count infor

    """
    
    
    try:
        trans_map = pd.read_csv(trans_int_map_file, sep = ' ')
    except:
        sys.exit("%s does not exist... check your count files.\n"%trans_map) 
    
    trans_map = pd.read_csv(trans_int_map_file, sep = ' ')
    trans_map = trans_map.replace({np.nan: None}) 
    trans_map.index = trans_map.Transcript    # use transcript id as index to facilitate search
    introns_dic = {}
    exons_dic = {}                            
    unmapped_transcript = set()
    
    
    base_names = []
    for name in samples:
        try:
            df = pd.read_csv(name, sep = '\t')
            if normalize_count:
                df = df[df.normalized_count > threshold].reset_index(drop = True)   # drop trancript expressed too low
            else:
                df = df[df.TPM > threshold].reset_index(drop = True)
        except:
            sys.exit("%s does not exist... check your count files.\n"%name) 
          
        base_name = name.split("/")[-1]
        base_names.append(base_name)

        _ = df.apply(transcript_to_intron_counts, axis = 1, args=(trans_map,base_name, introns_dic, exons_dic, normalize_count, unmapped_transcript))

    dic_to_csv(introns_dic, f'{out_prefix}count_intron', base_names)  # save the dic in csv format
    dic_to_csv(exons_dic, f'{out_prefix}count_exon', base_names)
    not_mapped_transcript = open(f'{out_prefix}dropped_trancripts', 'w')
    for droped in unmapped_transcript:
        print(droped, file = not_mapped_transcript)
    not_mapped_transcript.close()
    
    
    
    if len(unmapped_transcript) == 0:
        print('All transcript mapped to intron')
    
    elif len(unmapped_transcript) <= 500:
        print(f'There are {len(unmapped_transcript)} transcript not found in Transcript to intron map, like due to transcript id difference or filtering step')
    else:
        print(f'There are {len(unmapped_transcript)} transcript not found in Transcript to intron map, please check the version of annotation or may due to filtering step')
    
    return pd.read_csv(f'{out_prefix}count_intron', sep = ' ')




@timing_decorator
def build_init_cluster(intron_count_file):
    """
    Parameters
    ----------
    intron_count_file : str
        Path to a CSV file containing information about introns and their counts.
        Sample row format: Chr Start End Gene sample_1 sample_2........

    Adds a 'Cluster' column to the CSV with cluster assignments for introns.
    """
    
    # Read in the file
    df = pd.read_csv(intron_count_file, sep=' ')

    # Sort the DataFrame as required
    df.sort_values(by=['Chr', 'Start', 'End'], inplace=True)

    # Initialize variables for tracking clusters
    cluster_num = 1
    cluster_gene = df.iloc[0]['Gene'] 
    cluster_end = df.iloc[0]['End']
    cluster_chr = df.iloc[0]['Chr']

    # Function to apply to each row to determine cluster
    def check_new_cluster(row):
        nonlocal cluster_num
        nonlocal cluster_chr
        nonlocal cluster_end
        nonlocal cluster_gene #used in previous version, plan to removed in further version
        

        # Start a new cluster if the gene changes or there's no overlap with the current cluster
        # no need to check for sttart as we sort value based on Start
        if row['Chr'] != cluster_chr or row['Start'] > cluster_end:
            cluster_num += 1
            cluster_gene = row['Gene']
            cluster_end = row['End']
            cluster_chr = row['Chr']
            
        else:
            cluster_end = max(cluster_end, row['End'])
        return f'cluster_{cluster_num}'

    # Apply the function to each row in the DataFrame
    # df['Cluster'] = df.apply(check_new_cluster, axis=1)
    df.insert(1, 'Cluster', df.apply(check_new_cluster, axis=1))

    # Save the modified DataFrame back to a file
    df.to_csv(intron_count_file, sep=' ', index=False)






@timing_decorator
def process_clusters(init_clus_file, exon_count_file, intron_to_exon_file, out_prefix = '', mode = 3,
                     cutoff = 0.1, percent_cutoff = 0.01, min_cluster_val = 1):
    """
    
    Parameters
    ----------
    init_clus_file : Str, the file contain the init cluster result
    exon_count_file : Str, the file contain exon counts
    intron_to_exon_file : Str, the file contain information about          

    out_prefix : prefix of output file
    mode: mode that use to refine clusters, 1 for overlap, 2 for overlap+share_intron_splice_site, \
        3 for overlap+share_intron_splice_site+shared_exon_splice_site

    Returns
    -------
    None.

    """
    # test run code: process_clusters('count_intron', 'count_exon', 'intron_to_exon.tsv')
    
    # format of a cluster: {start: ; end: ; splicesites: ; exon_splicesites: ; introns: [[start,end,count,name, [exon_slicesites]]}}

    
    exon_count = pd.read_csv(exon_count_file, sep = ' ', index_col = 0 )
    initial_clus = pd.read_csv(init_clus_file, sep = ' ', index_col = 0)
    intron_to_exon = pd.read_csv(intron_to_exon_file, sep = ' ', index_col = 0)
    
    num_cluters = 0
    samples = list(initial_clus.columns[5:])
    output = open(f'{out_prefix}refined_cluster', 'w')
    
    out_str = ''
    for sample in samples:
        out_str += f'{sample} '
    out_str = out_str[:-1]
    print(out_str, file = output)

    # Sum the values for the samples. It will create a columns at -1 position 
    initial_clus['Sum'] = initial_clus[samples].sum(axis=1)
    clu_start = 0
    num_rows = len(initial_clus)
    
    
    
    for idx in range(num_rows + 1): # +1 to deal with -1 row have a diff clu than -2 row
    
        if idx != num_rows and initial_clus.iloc[clu_start]['Cluster'] == initial_clus.iloc[idx]['Cluster']:
            
            continue
        
        else: #
            if idx == num_rows:
                cluster_df = initial_clus.iloc[clu_start:]
            
            else: # there are more rows left
                cluster_df = initial_clus.iloc[clu_start:idx]
        
            clu_start = idx # update idx
            
            if len(cluster_df) < 2: # skip the clu that only contain one intron
                continue
            
            introns = []
            
            # build intron from initial cluster df
            for index, row in cluster_df.iterrows():
                intron = build_intron(index, row, samples, exon_count, intron_to_exon)
                # in format [start,end,count,name, {exon_slicesites}]
                introns.append(intron)
            
            
            cluster_dic = {}
            
            # perform filtering and reclu
            if mode == 1: #overlap only
                process_clu(introns, cluster_dic, cutoff, percent_cutoff, mode)
            else: 
                new_clus = []
                
                if mode == 2: #shared splice site
                    new_clus = refine_links(introns, exon_connection = False)
                else: 
                    new_clus = refine_links(introns, exon_connection = True)
            
                for clu in new_clus:
                    process_clu(clu, cluster_dic, cutoff, percent_cutoff, mode)
            
            

            # we should have a complete cluster_dic at this point
            for cluster_num in cluster_dic:
                clu_val = 0
                out_str = ''
                
                
                for clu_intron in cluster_dic[cluster_num]:
                    
                    infor = initial_clus.loc[clu_intron[3]]
                    strand = str(intron_to_exon.loc[clu_intron[3]]['strand'])
                    
                    out_str += f'{infor.Chr}:{infor.Start}:{infor.End}:clu_{num_cluters}_{strand}'
                    for value in list(infor)[5:-1]: # start of sample values, and not include the sum columns
                        out_str += f' {value}'
                        clu_val += float(value)
                    
                    out_str += '\n'
            
                if clu_val <= min_cluster_val:
                    continue
                out_str = out_str[:-1] # get rid of the last \n
                print(out_str, file = output)
                    
          
                num_cluters += 1

    output.close()









def build_intron(index, row, samples, exon_count,intron_to_exon):
    """
    a helper function for process_clusters function

    Parameters
    ----------
    index: the index of the df row, should be the name of intron
    row : a df row in format (Cluster Chr Start End Gene samples)
    samples: a list of sample names
    exon_count : df
    intron_to_exon : df

    Returns
    -------
    None.

    """
    intron = [row['Start'], row['End']]
    
    total_val = row['Sum']

    intron.append(total_val)
    intron.append(index)
    
    
    exon_sites = set()
    
    if index in intron_to_exon.index and pd.notna(intron_to_exon.at[index,'near_exons']): 
        #add this check because virtual intron 0 and -1 will not have connected exon
        exons = intron_to_exon.at[index,'near_exons'].split(',')
    else:
        exons = []
    
    
    for exon in exons:
        if exon in exon_count.index: # check whether this exon exist:
            exon_infor = exon_count.loc[exon]
            exon_sites.add(exon_infor['Start'])
            exon_sites.add(exon_infor['End'])
            
    intron.append(exon_sites)
    
    return intron
    




def process_clu(clu, output_dic, cutoff, percent_cutoff, mode):
    """ 
    this is the helper function of process_clusters 
    clu: a list that contain intron in format [[start,end,count,name, {exon_slicesites}]]
    output_dic: a dic that contain final cluster information
    cutoff: float or int, cutoff value for an intron
    percent_cutoff: float, percent cutoff value
    mode: str, the processing mode of the cluster
    
    """
    
    introns, reclu = filter_introns(clu, cutoff, percent_cutoff)
    
    if reclu == True and len(introns) > 1:
        
        clus = cluster_intervals(introns)
        
        
        for cluster in clus:
            
            if mode == 'overlap':
                process_clu(cluster,output_dic, cutoff, percent_cutoff, mode)
            
            elif mode == 'exon_connection' : 
                temp_clus = refine_links(cluster, True)
            
                for i in temp_clus:
                    process_clu(i,output_dic, cutoff, percent_cutoff, mode)
            else:
                temp_clus = refine_links(cluster, False)
            
                for i in temp_clus:
                    process_clu(i,output_dic, cutoff, percent_cutoff, mode)
            
    else: 
        if len(introns) > 1: # skip cluster that is less than 2 intron
            num_clus = len(output_dic)
            output_dic[num_clus] = introns
    
    
    
    
    


def filter_introns(introns, cutoff, percent_cutoff):
    """ 
    this is the helper function of process_clusters, will return a new introns list
    and True if reclu needed, flase if not needed
    introns:a list in format [[start,end,count,name, {exon_slicesites}]]
    cutoff: int or float
    percent_cutoff: float
    """
    
    total_TPM = 0
    for intron in introns: #loop over to get total TPM
        total_TPM += intron[2]
    
    new_introns = []
    reclu = False
    
    for intron in introns: 
        
        count = intron[2]
        if (count/total_TPM >= percent_cutoff and count >= cutoff):
            new_introns.append(intron)
        else:
            reclu = True
    

    return new_introns, reclu



def cluster_intervals(introns):
    """ 
    this is the helper function of process_clusters, it will generate a list of introns clu
    """
    introns.sort()
    
    clusters = []
    


    tmp_clu = [introns[0]] # init with the first intron
    current_intron = introns[0]    
    
    tmp_start = current_intron[0]
    tmp_end = current_intron[1]
    
    for i in range(1,len(introns)):
        
        current_intron = introns[i]
        
        if overlaps([tmp_start, tmp_end], [current_intron[0], current_intron[1]]):
            tmp_clu.append(current_intron)
        else:
            clusters.append(tmp_clu)
            tmp_clu = [current_intron]
        tmp_start = current_intron[0]
        tmp_end = max(tmp_end, current_intron[1])
            

    if len(tmp_clu) > 0:
        clusters.append(tmp_clu)



    return clusters




def overlaps(A,B):
    '''
    Checks if range A overlap with range B
    range in format [start,end]
    '''

    if A[1] < B[0] or B[1] < A[0]:
        return False
    else: return True




def refine_links(introns, exon_connection = True):

    """ 
    this is the helper function of process_clusters, it will generate a list of Intron_cluster object
    input should be a list of introns in format [[start,end,count,name, {exon_slicesites}]]
    """
    
    unassigned = introns[1:]
    current = [introns[0]]
    

    splicesites = set([current[0][0],current[0][1]])
    exon_splicesites = set(current[0][4])
    newClusters = []
    
    
    while len(unassigned) > 0:
        finished = False

        while not finished:
            finished = True
            torm = []
            for intron in unassigned:
                count = intron[2]
                start = intron[0]
                end = intron[1]
                
                
                
                if start in splicesites or end in splicesites:
                    current.append(intron)
                    splicesites.add(start)
                    splicesites.add(end)
                    finished = False
                    torm.append(intron)
                    exon_splicesites = exon_splicesites | intron[4] # union of exon_sites
                    
                    
                    
                elif exon_connection == True: 
                    if len(intron[4] & exon_splicesites) > 0:
                        current.append(intron)
                        splicesites.add(start)
                        splicesites.add(end)
                        finished = False
                        torm.append(intron)
                        exon_splicesites = exon_splicesites | intron[4] # union of exon_sites
               
                    
               
            for intron in torm:
                unassigned.remove(intron)

                
                
        newClusters.append(current)
        current = []
        if len(unassigned) > 0:
            current = [unassigned[0]]
            splicesites = set([current[0][0],current[0][1]])
            exon_splicesites = set(current[0][4])
            unassigned = unassigned[1:]
    
    if len(current) > 0:
        newClusters.append(current) #this line is not necessary, it is use to append the single intron cluster tha come last
    

    return newClusters   




@timing_decorator
def compute_ratio(cluster_file, output_prefix=''):
    """

    Parameters
    ----------
    cluster_file : the refined cluster file from process_clusters
    output_prefix: str

    Returns
    -------
    None
    save to disk, result in a file called ratio_count

    """
    df = pd.read_csv(cluster_file, sep=' ', index_col=0)
    samples = df.columns
    
    # Split the index to get the cluster information
    df['Cluster'] = df.index.to_series().apply(lambda x: x.split(':')[3])
    
    # Calculate the sum of each sample within each cluster
    cluster_sums = df.groupby('Cluster').sum()
    
    # Join the cluster sums back to the original dataframe
    df = df.join(cluster_sums, on='Cluster', rsuffix='_sum')
    
    
    # Calculate ratios in a vectorized way
    for sample in samples:  
        df[sample] = df[sample] / df[f'{sample}_sum']
    
    # Prepare the output dataframe, dropping sum columns and the Cluster column
    df.fillna(0, inplace=True)
    df = df.round(6)
    output_df = df[[col for col in df.columns if '_sum' not in col and col != 'Cluster']]
    # Exclude the last two columns which are 'Cluster' and sums
    
    # Write the result to a file
    output_df.to_csv(f'{output_prefix}ratio_count', sep=' ')

    #Remove the leading space of the first line, this step is specifically to fit the leafcutter_df format
    with open(f'{output_prefix}ratio_count', 'r') as file:
        lines = file.readlines()
      
    lines[0] = lines[0].lstrip()

    # Write the modified content back to the file
    with open(f'{output_prefix}ratio_count', 'w') as file:
        file.writelines(lines)



def input_file_processing(input_file):
    """
 
    Parameters
    ----------
    input_file : a string that direct to the file that contain the name of samples

    Returns
    -------
    a list that contain the name of samples

    """
    
    result = []
    with open(input_file, 'r') as file:
        for line in file:
            result += [line[:-1]] # [:-1] to remove the \n 
    
    return result




def count_TPM_normalization(quant_file, gft_file, out_prefix= None):
    # Your function implementation
    df = pd.read_csv(quant_file, sep='\t')
    df = df[df['TPM'] > 0]
    
    transcript_to_gene_dic = extract_transcript_to_gene_map(gft_file)
    
    gene_dic = {}
    df['gene'] = df['Name'].map(transcript_to_gene_dic)
    
    for index, row in df.iterrows():
        TPM = row['TPM']
        count = row['NumReads']
        gene = row['gene']
        
        gene_dic.setdefault(gene, {'TPM': 0, 'count': 0})
        gene_dic[gene]['TPM'] += TPM
        gene_dic[gene]['count'] += count

    df['normalized_count'] = df.apply(
        lambda x: (x['TPM'] / gene_dic[x['gene']]['TPM']) * gene_dic[x['gene']]['count'], axis=1)
    
    # Use the out_prefix to determine the output file name
    output_file = f'{out_prefix}_normalized.sf' if out_prefix else 'normalized.sf'
    df.to_csv(output_file, sep='\t', index=False)




def extract_transcript_to_gene_map(gtf_file):
    """
    generate a map that map transcript_id to gene_id

    Parameters
    ----------
    gtf_file : a string that contain the name for the gtf file

    Returns
    -------
    map_transcript_to_gene : a dic that map transcript_id to gene_id

    """
    
    map_transcript_to_gene = {}

    with open(gtf_file, 'r') as file:
        for line in file:
            if line.startswith('#'):
                continue  # Skip header lines

            fields = line.strip().split('\t')
            
            if fields[2] == 'transcript':
                attributes = fields[8]
                attribute_parts = attributes.split(';')

                transcript_id = ''
                gene_id = ''

                for attribute in attribute_parts:
                    if 'transcript_id' in attribute:
                        transcript_id = attribute.split('"')[1]
                    elif 'gene_id' in attribute:
                        gene_id = attribute.split('"')[1]

                if transcript_id and gene_id:
                    map_transcript_to_gene[transcript_id] = gene_id

    return map_transcript_to_gene



@timing_decorator
def LeafcutterITI_clustering(options):
    """
    This is the main function for LeafcutterITI
    """

    samples = input_file_processing(options.count_files)
    
    # perform normalization
    if options.normalization == True and options.preprocessed == False:
        new_samples = []
        for sample in samples:
            prefix = sample.split(".")[0] # get rid of the .sf
            count_TPM_normalization(sample, options.annot, prefix)
            new_samples += [f'{prefix}_normalized.sf']
               
        samples = new_samples
        # saved the normalized sample for furture usage
        with open(f'{options.count_files.split(".")[0]}_normalized.txt', 'w') as file:
            for sample in samples:
                file.write(f"{sample}\n")
            
        sys.stderr.write("Finished Normalization\n")
        

    out_prefix = options.outprefix
    
    # count_introns
    count_introns(samples, \
                  options.map, \
                  out_prefix= out_prefix,\
                  threshold = options.samplecutoff,\
                  normalize_count = options.normalization)
    sys.stderr.write("Finished Introns Counting\n")
    
    
    # build initial cluster
    build_init_cluster(f'{out_prefix}count_intron')
    sys.stderr.write("Finished Initial Clustering\n")


    process_clusters(f'{out_prefix}count_intron', f'{out_prefix}count_exon', \
                     options.connect_file, \
                     out_prefix = options.outprefix,\
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
        

   
    LeafcutterITI_clustering(options)









# not in use
class Intron_cluster:
    # a class for introns cluster to improve readability of code
    def __init__(self, intron):
        # intron in format [start,end,count,name, [exon_slicesites]]
        
        
        self.start = intron[0]
        self.end = intron[1]
        self.splicesites = set()
        self.splicesites.add(intron[0])
        self.splicesites.add(intron[1])
        
        self.exon_splicesites = set()
        for site in intron[4]:
            self.exon_splicesites.add(site)
        
        self.introns = [intron]
        
    def change_end(self,new_end):
        
        if self.end < new_end:
            self.end = new_end

    def change_start(self, new_start):
        
        if self.start > new_start:
            self.start = new_start








































        
    




