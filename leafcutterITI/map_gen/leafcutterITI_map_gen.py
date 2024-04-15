import pyranges as pr
import numpy as np
import pandas as pd
import sys
import warnings
import leafcutterITI.utils
from leafcutterITI.utils import timing_decorator,write_options_to_file

from optparse import OptionParser

warnings.simplefilter(action='ignore', category=pd.errors.DtypeWarning)



chr_dic = {'chr1': 0,'chr2': 1,'chr3': 2,'chr4': 3,'chr5': 4,'chr6': 5,'chr7': 6,'chr8': 7,'chr9': 8,'chr10': 9,
               'chr11': 10,'chr12': 11,'chr13': 12, 'chr14': 13,'chr15': 14,'chr16': 15,'chr17': 16,'chr18': 17,
               'chr19': 18,'chr20': 19,'chr21': 20, 'chr22': 21, 'chrX':22, 'chrY':23}


##################################isoform to intron and exon map generation
def compute_transcript_intron_map(annot, out_prefix = '', annot_type = 'gencode', min_length = 50, max_length = 5000000, quality_control = True):
    """

    Parameters
    ----------
    annot : a gtf format file that contain the annotation for the genome
    out_prefix : the output prefix, example: "mm10_"
    annot_type : support gencode, Stringtie (in progress, will add support to ref_seq, ensembl)
    min_length: int, the minimun intron length for intron to be consider 
    max_length: int, the maximun intron length to be consider
    quality_control: whether to get rid of pseudogene, decay transcript

    Returns
    -------
    None.
    Result will be saved to disk

    """

    num_gene = 0
    num_transcript = 0
    
    
    try:
        open(annot)
    except:
        sys.exit("%s does not exist... check your annotation files.\n"%annot)
        
    
    gr = pr.read_gtf(annot)
    df = gr.df
    df = df.replace({np.nan: None})
    df = df[df['Chromosome'].str.contains('chr')].reset_index(drop=True) #eliminate gene and transcript not mapped into a chr
    df = df[~ df['gene_type'].str.contains('artifact', na = False)]
    df = df[df['Chromosome'] != ('chrM')].reset_index(drop=True) #eliminate mito genes as they are not likely undergo alternative splicing
    # deal with output from stringtie
    if annot_type == 'Stringtie':
        gene_dic = {}
        name_dic = {}


        filtered_df = df[df['ref_gene_id'].notna()]

        # Use `set` to get unique gene_ids and update dictionaries
        for gene_id, gene_name, ref_gene_id in zip(filtered_df['gene_id'], filtered_df['gene_name'], filtered_df['ref_gene_id']):
            if gene_id not in gene_dic or gene_id not in name_dic:
                name_dic[gene_id] = gene_name
                gene_dic[gene_id] = ref_gene_id

        df['ref_gene_id'] = df['gene_id'].map(gene_dic).fillna(df['ref_gene_id'])
        df['gene_name'] = df['gene_id'].map(name_dic).fillna(df['gene_name'])   
    

        reference_list = open(f'{out_prefix}stringtie_transcriptome_reference.tsv', 'w') 
        print('Stringtie_id reference_id reference_name', file=reference_list)
        for key in gene_dic:
            print(f'{key} {gene_dic[key]} {name_dic[key]}', file = reference_list)
        reference_list.close()

    elif quality_control == True:
        df = df[~ df['gene_type'].str.contains('pseudogene', na = False)]
        df = df[~ df['transcript_type'].str.contains('decay', na = False)]
        # df = df[~ df['transcript_type'].str.contains('retained_intron', na = False)]
        # maybe useful to account for retain_intron isoforms


    genes = list(df['gene_id'].unique()) #get all gene_id    

    gene_start_index = 0
    num_rows = len(df) - 1


    output = open(f'{out_prefix}isoform_intron_map.tsv', 'w')

    if annot_type != 'Stringtie':
        print('Chr Gene Transcript support_introns support_exons Transcript_type', file = output)
    else: 
        print('Chr Gene Transcript support_introns support_exons', file = output)

    
    near_exon_dic = {}  # the dic that keep infor about the nearby exon of a intron
    # should in format {intron: {exon: counts}} counts not use 
    intron_str_dic = {}

    for idx in range(len(df)): # don't need to deal with the case that -1 row have a different gene than -2 row, since a gene must have at lease one transcript
    
        if df.iloc[gene_start_index]['gene_id'] == df.iloc[idx]['gene_id'] and idx != num_rows:
            continue
        else: #
            if idx == num_rows:
                gene_df = df.iloc[gene_start_index:]
            
            else: # there are more rows left
                gene_df = df.iloc[gene_start_index:idx]
        
        
        
            num_gene += 1
            gene = gene_df['gene_id'].unique()[0]
            chromosome = gene_df['Chromosome'].unique()[0]
            transcripts = list(gene_df['transcript_id'].unique()) #get all gene_id, will include None 
        
        
            for trans in transcripts:
                if trans == None: # skip the None value 
                    continue 
                num_transcript += 1
                exons = []
                trans_df = gene_df[(gene_df['transcript_id'] == trans) & (gene_df['Feature'] == 'exon')] # only care about exon for the transcript
                if annot_type != 'Stringtie': 
                    trans_type = trans_df['transcript_type'].unique()[0]
        
                for index,row in trans_df.iterrows():
                    exon = [int(row['exon_number']), row['Chromosome'], row['Start'], row['End']] # exon number is the order of exon for this transcript
                    
                    exons += [exon]
        
                # computer introns:
                exons = sorted(exons)  # guranttee that exon will follow the order
                introns = []
                for i in range(len(exons) - 1):  # -1 as only len(exons) - 1 introns exist 
                    current_exon = exons[i] 
                    next_exon = exons[i + 1]
                    current_exon_name = f'{current_exon[1]}:{current_exon[2]}-{current_exon[3]}'
                    next_exon_name = f'{next_exon[1]}:{next_exon[2]}-{next_exon[3]}'
                    
            
                    if current_exon[3] < next_exon[2]: # current_exon end < next_exon start, foward direction
                        introns += [[i+1,current_exon[1], (current_exon[3]), next_exon[2]]]
                        
                        # the precise location could be problematic, but it is defined th
                        temp_name = f'{current_exon[1]}:{current_exon[3]}-{next_exon[2]}'
                        
                        add_near_exon_dic(near_exon_dic, temp_name, current_exon_name, next_exon_name) # add value to near_exon_dic
                        # use helper function to improve readbility
                        
                        intron_str_dic[temp_name] = '+'
                        
                
                    else: # current_exon start > next_exon end, backward direction
                        introns += [[i+1,current_exon[1], (next_exon[3]), current_exon[2]]]
                        
                        temp_name = f'{current_exon[1]}:{next_exon[3]}-{current_exon[2]}'
                
                        add_near_exon_dic(near_exon_dic, temp_name, current_exon_name, next_exon_name)
                        
                        
                        intron_str_dic[temp_name] = '-'
        
        
                # the order of intron doesn't really matter, but will keep in case have usage in future 
            
                out_introns = ''
                out_exons = ''
                for i in introns: 
                    intron_len = abs(i[2]-i[3])
                    if intron_len > min_length and intron_len < max_length: # limit the length of intron
                        out_introns += f'{chromosome}:{i[2]}-{i[3]},'
                
            
                for i in exons:
                    out_exons += f'{chromosome}:{i[2]}-{i[3]},'
            
                out_introns = out_introns[:-1]  
                out_exons = out_exons[:-1] # eliminate the last ","
                if annot_type != 'Stringtie':
                    print(f'{chromosome} {gene} {trans} {out_introns} {out_exons} {trans_type}', file = output)
                else: 
                    print(f'{chromosome} {gene} {trans} {out_introns} {out_exons}', file = output)
            gene_start_index = idx
        
    output.close()
    print_near_exon_dic(near_exon_dic, intron_str_dic, out_prefix)
    
    sys.stderr.write("\nFinished process %d genes and %d isofroms\n"% (num_gene, num_transcript))



def add_near_exon_dic(near_exon_dic, intron_name, current_exon_name, next_exon_name):
    """ 
    this is the helper function to build near_exon_dic
    """
    
    if intron_name not in near_exon_dic:
        near_exon_dic[intron_name] = {current_exon_name:1}
        near_exon_dic[intron_name][next_exon_name] = 1
                            
    else:
        if current_exon_name not in near_exon_dic[intron_name]:
            near_exon_dic[intron_name][current_exon_name] = 1
        else: 
            near_exon_dic[intron_name][current_exon_name] += 1
                        
        if next_exon_name not in near_exon_dic[intron_name]:
            near_exon_dic[intron_name][next_exon_name] = 1
        else: 
            near_exon_dic[intron_name][next_exon_name] += 1
                        

def print_near_exon_dic(near_exon_dic, intron_str_dic, out_prefix):
    """ 
    this is the ulepre function that print out the 
    """
    
    output = open(f'{out_prefix}intron_exon_connectivity.tsv', 'w')

    print('intron near_exons strand', file = output)

    for intron in near_exon_dic:
        exons = ''
        for exon in near_exon_dic[intron]:
            exons += f'{exon},'
        exons = exons[:-1]
        print(f'{intron} {exons} {intron_str_dic[intron]}',file = output)
    
    output.close()
    
##############################################################################################

#these function is some additional features, not fully tested yet

def add_virtual_first_last_introns(trancript_to_intron_map, intron_to_exon_map, out_prefix = '', include_exon = False):
    """
    This is the function that use to add virtual intron 0 and -1 to capture 
    AFE and ALE usage. Will use a location before or after all transcript to generate 
    a virtual exon and compute virtual intron between it and the first exon of the 
    transcript

    """
    df = pd.read_csv(trancript_to_intron_map, sep = ' ')
    df = df.where(pd.notna(df), None)
    gene_dic = {}
    intron_exon_df = pd.read_csv(intron_to_exon_map, sep = ' ')
    
    
    # this round to get the start and end for a gene
    for index, row in df.iterrows():
        current_gene = row['Gene']
        current_exons= row['support_exons'].split(',')
        
        
        
        # exon in format "chr1:3143475-3144545"
        first_exon = current_exons[0].split(':')[1].split('-')
        last_exon = current_exons[-1].split(':')[1].split('-')
        

        # tmp_start always smaller than tmp_end
        if int(first_exon[0]) <= int(last_exon[0]): # forward direction
            tmp_start = int(first_exon[0])
            tmp_end = int(last_exon[1]) 
            tmp_direction = '+'
        else: # reverse direction
            tmp_start = int(last_exon[0]) 
            tmp_end = int(first_exon[1])
            tmp_direction = '-'
        
        
        if current_gene not in gene_dic:
            gene_dic[current_gene] = {'start':tmp_start, 'end': tmp_end, 'direction': tmp_direction}
            
        else:
            if tmp_start < gene_dic[current_gene]['start']:
                gene_dic[current_gene]['start'] = tmp_start
            if tmp_end > gene_dic[current_gene]['end']:
                gene_dic[current_gene]['end'] = tmp_end
            
            if tmp_direction != gene_dic[current_gene]['direction'] and len(current_exons) >1:
                #correct the possible direction mistake when single exon appear
                gene_dic[current_gene]['direction'] = tmp_direction
        
    list_virtual = open(f'{out_prefix}virtual.tsv', 'w')
    print('intron type', file = list_virtual)
    
    #not sure whether to update the intron_to_exon map as virtual intron may cause
    # arbitray large cluster depend on reference structure, but without it hard to
    # define hybrid exon
    intron_exon_dic = {}
    intron_dic = {}
    
    for index, row in df.iterrows():
        current_chr = row['Chr']
        current_gene = row['Gene']
        current_exons= row['support_exons'].split(',')
        # exon in format "chr1:3143475-3144545"
        first_exon = current_exons[0].split(':')[1].split('-')
        last_exon = current_exons[-1].split(':')[1].split('-')
        
        
        
        direction = gene_dic[current_gene]['direction']
        
        if direction == '+':
            intron_0 = f'{current_chr}:{gene_dic[current_gene]["start"] - 500}-{int(first_exon[0])}'
            intron_minus1 = f'{current_chr}:{int(last_exon[1])}-{gene_dic[current_gene]["end"] + 500}'
            
            
            
        else:
            intron_0 = f'{current_chr}:{int(first_exon[1])}-{gene_dic[current_gene]["end"] + 500}'
            intron_minus1 = f'{current_chr}:{gene_dic[current_gene]["start"] - 500}-{int(last_exon[0])}'
            
    
        
        if include_exon == True:
            if intron_0 not in intron_exon_dic:
                intron_exon_dic[intron_0] = {'exons': set([current_exons[0]]), 'strand': direction}
            else: 
                intron_exon_dic[intron_0]['exons'].append(current_exons[0])
                
            if intron_minus1 not in intron_exon_dic:
                intron_exon_dic[intron_minus1] = {'exons': set([current_exons[-1]]), 'strand': direction}
            else: 
                intron_exon_dic[intron_minus1]['exons'].append(current_exons[-1])  
        else:
            if intron_0 not in intron_exon_dic:
                intron_exon_dic[intron_0] = {'exons': 'not_avalible', 'strand': direction}
            if intron_minus1 not in intron_exon_dic:
                intron_exon_dic[intron_minus1] = {'exons': 'not_avalible', 'strand': direction}
            # can use N/A or None to replace 'not_avalible'
            
            
        
        if df.loc[index]['support_introns'] == None:
            df.loc[index]['support_introns'] = f'{intron_0},{intron_minus1}'
        else:
            df.loc[index]['support_introns'] = f'{intron_0},{df.loc[index]["support_introns"]},{intron_minus1}'
        
        
        if intron_0 not in intron_dic:
            intron_dic[intron_0] = 'intron_0'
        if intron_minus1 not in intron_dic:
            intron_dic[intron_minus1] = 'intron_-1'
        
        
    intron_exon_list = []
    for intron in intron_exon_dic:
        if include_exon == True:
            tmp_str = ''
            for element in intron_exon_dic[intron]['exons']:
                tmp_str += f'{element},'
            tmp_str = tmp_str[:-1]
            # get rid of last , 
            intron_exon_list.append([intron,tmp_str,intron_exon_dic[intron]['strand']])
        else:
            intron_exon_list.append([intron,intron_exon_dic[intron]['exons'],intron_exon_dic[intron]['strand']])
    
    for intron in intron_dic:
        print(f'{intron} {intron_dic[intron]}', file=list_virtual)
    

    list_virtual.close()
    df.to_csv(f'{out_prefix}isoform_intron_map_with_virtual.tsv', sep = ' ', index=False)
    intron_exon_df = intron_exon_df.append(pd.DataFrame(intron_exon_list, columns = ['intron', 'near_exons', 'strand']), ignore_index=True)
    intron_exon_df.to_csv(f'{out_prefix}intron_exon_connectivity_with_virtual.tsv', sep = ' ', index=False)    




def intron_source_generation(transcript_intron_map, out_prefix = ''):
    df = pd.read_csv(transcript_intron_map, sep = ' ')
    df = df[df['support_introns'].notna()]
    df['support_introns'] = df['support_introns'].str.split(',')
    df = df.explode('support_introns')
    df = df.rename(columns={'support_introns': 'intron', 'Transcript_type': 'source_type'})
    df = df[['intron', 'source_type']]
    # we now get a intron source type map, but there are duplucation
    df = df.groupby('intron')['source_type'].agg(lambda x: ','.join(set(x))).reset_index()
    # get rid of duplication of intron and merge all source_type that can mapped to an intron
    
    df.to_csv(f'{out_prefix}intron_source_map.tsv')



#######################################################################################################################################
@timing_decorator
def LeafcutterITI_map_generation(options):
    """
    This this the main warper function for the LeafcutterITI map generation

    """
    
    
    
    compute_transcript_intron_map(options.annot, out_prefix= options.outprefix,\
                                  annot_type = options.annot_source,\
                                  min_length = options.minintronlen,\
                                  max_length= options.maxintronlen, \
                                  quality_control= options.quality_control)
        
    out_prefix = options.outprefix    
    if options.annot_source == 'gencode':
        intron_source_generation(f'{out_prefix}isoform_intron_map.tsv',out_prefix)
        
    if options.virtual_intron == True:
        
        f'{out_prefix}intron_exon_connectivity.tsv'
        f'{out_prefix}isoform_intron_map.tsv'

        add_virtual_first_last_introns(f'{out_prefix}isoform_intron_map.tsv', f'{out_prefix}intron_exon_connectivity.tsv',\
                                        out_prefix)





if __name__ == "__main__":

    parser = OptionParser()

    parser.add_option('-a',"--annot", dest="annot", default = None,
              help="annotation source, support gencode, Stringtie \
                  (in progress, will add support to ref_seq, ensembl)")    
                  
    parser.add_option("--annot_source", dest="annot_source", default = 'gencode',
              help="annotation source, support gencode, Stringtie \
                  (in progress, will add support to ref_seq, ensembl)")    
                 
    parser.add_option("-o", "--outprefix", dest="outprefix", default = 'leafcutter_',
                  help="output prefix (default leafcutter_), should include the diretory address if not\
                  in the same dic")    

    parser.add_option("--maxintronlen", dest="maxintronlen", default = 5000000,
                  help="maximum intron length in bp (default 500,000bp)")
        
    parser.add_option("--minintronlen", dest="minintronlen", default = 50,
                  help="minimum intron length in bp (default 50bp)")

    parser.add_option("--quality_control", dest="quality_control", default = True,
                 help="whether to remove pseudogene, and decay transcript")

    parser.add_option("-v", "--virtual_intron", dest="virtual_intron", default = False,
                 help="whether to use virtual intron to capture alternative first and last exon usage")
   
    (options, args) = parser.parse_args()



    if options.annot is None:
        sys.exit("Error: no annotation file provided...\n")


    sys.stderr.write(f"Start Processing transcriptome annotation {options.annot}")
    sys.stderr.write(f"annot source: {options.annot_source}")
    sys.stderr.write(f"outprefix: {options.outprefix}")
    sys.stderr.write(f"Max intron length: {options.maxintronlen}")
    sys.stderr.write(f"Min intron length: {options.minintronlen}")
    sys.stderr.write(f"remove pseudogene and decay transcript: {options.quality_control}")
    sys.stderr.write(f"virtual intron: {options.virtual_intron}")
    
    record = f'{options.outprefix}map_parameters.txt'
    sys.stderr.write(f'saving record to {record}')
    
    write_options_to_file(options, record)

    LeafcutterITI_map_generation(options)
    
    sys.stderr.write('Finish building Isofrom to intron map')











