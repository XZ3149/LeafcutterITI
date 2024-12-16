import subprocess as sp
import sys, re, copy, os, codecs, gzip
from argparse import ArgumentParser, Action as ArgParseAction
from collections import OrderedDict
import pandas as pd
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

__version__ = "1.1.5"


def get_version():
        """Return version information."""
        prog = 'ggsashimi'
        version = '{} v{}'.format(prog, __version__)
        return version

def gen_intron_exon_stand_dic(intron_exon_connect_file):
    df = pd.read_csv(intron_exon_connect_file, sep = ' ')
    
    # Initialize an empty dictionary
    strand_dict = {}

    # Iterate over each row in the DataFrame
    for index, row in df.iterrows():
        intron = row['intron']
        near_exons = row['near_exons']
        strand = row['strand']
        
        # Map the intron to its strand
    
        strand_dict[intron] = strand

        # Split the 'near_exons' string by ',' and map each exon to its strand
        exons = near_exons.split(',')
        for exon in exons:
            if exon not in strand_dict:
                strand_dict[exon] = strand

    return strand_dict

def process_txt_file(file_name):
    """
    This function is used for processing args.overlay and args.color
    Each file should be txt and contain 

    """
    
    mapping = {}

    with open(file_name, "r") as file:
        for line in file:
            # Strip whitespace (e.g., newline at the end)
            line = line.strip()
            # Split the line into two parts: sample and color
            sample, value = line.split()
            # Store in the dictionary
            mapping[sample] = value
    return mapping
    
    
def get_celltype(colname, sample_cell_type_dic, normalized, data_type = "bulk", ):
    parts = colname.split('_')
    if data_type == 'sc':
        cell_type =  '_'.join(parts[:-1]) # ASC_0
    elif normalized == "True":
        cell_type = '_'.join(parts[:-3]) #astrocyte_0_quant_normalized.sf
    elif normalized == "False":
        cell_type = '_'.join(parts[:-2]) #astrocyte_0_quant.sf
    else:
        print("ERROR: invalid value in normalized")
        exit(1)
        
    sample_cell_type_dic[colname] = cell_type
    return cell_type


def modify_index(idx):
    parts = idx.split(":")
    parts[2] = str(int(parts[2]) - 1)  # Subtract 1 from the "end" part
    return ":".join(parts)



def aggregate_cell_types(df, data_type = "bulk", normalized = 'True', aggregation = 'True'):
    
    if "Cluster" in df.columns:
        df.drop(columns = ['Cluster'], inplace = True)
    
    
    if 'Gene' in df.columns:
        df_info = df.iloc[:, :5]
        # Step 2: Process the remaining columns (cell type data)
        df_cells = df.iloc[:, 5:]
    else: 
        df_info = df.iloc[:, :4]
        df_cells = df.iloc[:, 4:]
    
    sample_cell_type_dic = {}
    # Extract cell type names for each column in df_cells
    celltype_names = [get_celltype(col, sample_cell_type_dic, normalized, data_type) for col in df_cells.columns]
    
    if aggregation == "True":

        # Assign a MultiIndex to the columns (cell type, original column name)
        df_cells.columns = pd.MultiIndex.from_arrays([celltype_names, df_cells.columns])
        # **Compute the mean for each cell type**
        df_celltype_mean = df_cells.groupby(level=0, axis=1).mean()

        # Step 3: Concatenate the first four columns back with the mean cell type data
        df_result = pd.concat([df_info, df_celltype_mean], axis=1)
        
        sample_cell_type_dic = {}
        celltypes = set(celltype_names)
        for name in celltypes:
            sample_cell_type_dic[name] = name
        
    elif aggregation == "False":
        
        df_result = df
    else:
        print("ERROR: invalid value in aggregation")
        exit(1)
    
    return df_result, sample_cell_type_dic




def parse_coordinates(c):
    
        chr = c.split(":")[0]
        start, end = c.split(":")[1].split("-")
        # tealeaf is 1-based
        start, end = int(start), int(end)
        return chr, start, end


def read_intron_exon_input(df_intron, df_exon, sample_cell_type_dic):
    """
    Overlay, color, label should be a dict that contain the value each celltype should be

    """
    if 'Gene' in df_intron.columns:
        samples = list(df_intron.columns[5:])
    else:
        samples = list(df_intron.columns[4:])
        
    for sample in samples:
        name = sample
        cell_type = sample_cell_type_dic[sample]
        
        
        overlay_level = cell_type
        color_level = cell_type
        label_text = cell_type

        
        
        yield name, cell_type, overlay_level, color_level, label_text

        
def read_intron_exon(df_intron, df_exon, cell_type, strand_dic, c, s = None):
    """
    This function take two df based on intron counts and exon count, also the cell type name
    to build the coverage and junction for a specific cell type
    """
    
    
    chrom, start, end = parse_coordinates(c)

    # Initialize coverage array and junction dict
    a = {"+" : [0] * (end - start) }
    junctions = {"+": OrderedDict()}
    if s != "NONE":
        a["-"] = [0] * (end - start)
        junctions["-"] = OrderedDict()
        
    
    
    # Filter the DataFrame for overlapping intervals
    overlapping_rows = df_intron[
        (df_intron['Chr'] == chrom) &
        (df_intron['Start'] <= end) &
        (df_intron['Start'] >= start) & 
        (df_intron['End'] - 1<= end) &
        (df_intron['End'] - 1 >= start) 
    ] # we need the junction to be fullly included in the region
      # -1 as tealeaf intron have + 1 to align with leafcutter

    if len(overlapping_rows) > 0:
        for index, row in overlapping_rows.iterrows():

            
            if s == "NONE": 
                strand = "+"
            else:
                strand = strand_dic[row['Name']]
            count = round(row[cell_type], 1)
            junctions[strand][(row['Start'], row['End'] - 1)] = count # -1 as tealeaf intron have + 1 to align with leafcutter
            
            
            
    overlapping_rows = df_exon[
        (df_exon['Chr'] == chrom) &
        (df_exon['Start'] <= end) &
        (df_exon['End'] >= start)
    ]
        
    if len(overlapping_rows) > 0:
        for index, row in overlapping_rows.iterrows():
            
            if s == "NONE": 
                strand = "+"
            else:
                strand = strand_dic[row['Name']]
            count = round(row[cell_type], 1)
            
            for i in range(row['Start'], row['End']):
                if i < start or i >= end:
                    continue
                ind = i - start
                a[strand][ind] += count
        
    return a, junctions


def prepare_for_R(a, junctions, c, m):

        _, start, _ = parse_coordinates(c)

        # Convert the array index to genomic coordinates
        x = list(i+start for i in range(len(a)))
        y = a

        # Arrays for R
        dons, accs, yd, ya, counts = [], [], [], [], []

        # Prepare arrays for junctions (which will be the arcs)
        for (don, acc), n in junctions.items():

                # Do not add junctions with less than defined coverage
                if n < m:
                        continue

                dons.append(don)
                accs.append(acc)
                counts.append(n)

                yd.append( a[ don - start -1 ])
                ya.append( a[ acc - start +1 ])

        return x, y, dons, accs, yd, ya, counts




def intersect_introns(data):
        data = sorted(data)
        it = iter(data)
        a, b = next(it)
        for c, d in it:
                if b > c:  # Use `if b > c` if you want (1,2), (2,3) not to be
                                # treated as intersection.
                        b = min(b, d)
                        a = max(a, c)
                else:
                        yield a, b
                        a, b = c, d
        yield a, b



def shrink_density(x, y, introns):
        new_x, new_y = [], []
        shift = 0
        start = 0
        # introns are already sorted by coordinates
        for a,b in introns:
                end = x.index(a)+1
                new_x += [int(i-shift) for i in x[start:end]]
                new_y += y[start:end]
                start = x.index(b)
                l = (b-a)
                shift += l-l**0.7
        new_x += [int(i-shift) for i in x[start:]]
        new_y += y[start:]
        return new_x, new_y

def shrink_junctions(dons, accs, introns):
        new_dons, new_accs = [0]*len(dons), [0]*len(accs)
        real_introns = dict()
        shift_acc = 0
        shift_don = 0
        s = set()
        junctions = list(zip(dons, accs))
        for a,b in introns:
                l = b - a
                shift_acc += l-int(l**0.7)
                real_introns[a - shift_don] = a
                real_introns[b - shift_acc] = b
                for i, (don, acc) in enumerate(junctions):
                        if a >= don and b <= acc:
                                if (don,acc) not in s:
                                        new_dons[i] = don - shift_don
                                        new_accs[i] = acc - shift_acc
                                else:
                                        new_accs[i] = acc - shift_acc
                                s.add((don,acc))
                shift_don = shift_acc
        return real_introns, new_dons, new_accs

def read_palette(f):
        palette = "#ff0000", "#00ff00", "#0000ff", "#000000"
        if f:
                with open(f) as openf:
                        palette = list(line.split("\t")[0].strip() for line in openf)
        return palette


def read_gtf(f, c):
        exons = OrderedDict()
        transcripts = OrderedDict()
        chr, start, end = parse_coordinates(c)
        end = end -1
        with gzip.open(f, 'rt') if f.endswith(".gz") else open(f) as openf:
                for line in openf:
                        if line.startswith("#"):
                                continue
                        el_chr, _, el, el_start, el_end, _, strand, _, tags = line.strip().split("\t")
                        if el_chr != chr:
                                continue
                        if el not in ("transcript", "exon"):
                                continue
                        try:
                                transcript_id = re.findall('transcript_id ("[^"]+")', tags)[0]
                        except KeyError:
                                print("ERROR: 'transcript_id' attribute is missing in the GTF file.")
                                exit(1)
                        el_start, el_end = int(el_start) -1, int(el_end)
                        strand = '"' + strand + '"'
                        if el == "transcript":
                                if (el_end > start and el_start < end):
                                        transcripts[transcript_id] = max(start, el_start), min(end, el_end), strand
                                continue
                        if el == "exon":
                                if (start < el_start < end or start < el_end < end):
                                        exons.setdefault(transcript_id, []).append((max(el_start, start), min(end, el_end), strand))

        return transcripts, exons


def make_introns(transcripts, exons, intersected_introns=None):
        new_transcripts = copy.deepcopy(transcripts)
        new_exons = copy.deepcopy(exons)
        introns = OrderedDict()
        if intersected_introns:
                for tx, (tx_start,tx_end,strand) in new_transcripts.items():
                        total_shift = 0
                        for a,b in intersected_introns:
                                l = b - a
                                shift = l - int(l**0.7)
                                total_shift += shift
                                for i, (exon_start,exon_end,strand) in enumerate(exons.get(tx,[])):
                                        new_exon_start, new_exon_end = new_exons[tx][i][:2]
                                        if a < exon_start:
                                                if b > exon_end:
                                                        if i ==  len(exons[tx])-1:
                                                                total_shift = total_shift - shift + (exon_start - a)*(1-int(l**-0.3))
                                                        shift = (exon_start - a)*(1-int(l**-0.3))
                                                        new_exon_end = new_exons[tx][i][1] - shift
                                                new_exon_start = new_exons[tx][i][0] - shift
                                        if b <= exon_end:
                                                new_exon_end = new_exons[tx][i][1] - shift
                                        new_exons[tx][i] = (new_exon_start,new_exon_end,strand)
                        tx_start = min(tx_start, sorted(new_exons.get(tx, [[sys.maxsize]]))[0][0])
                        new_transcripts[tx] = (tx_start, tx_end - total_shift, strand)

        for tx, (tx_start,tx_end,strand) in new_transcripts.items():
                intron_start = tx_start
                ex_end = 0
                for ex_start, ex_end, strand in sorted(new_exons.get(tx, [])):
                        intron_end = ex_start
                        if tx_start < ex_start:
                                introns.setdefault(tx, []).append((intron_start, intron_end, strand))
                        intron_start = ex_end
                if tx_end > ex_end:
                        introns.setdefault(tx, []).append((intron_start, tx_end, strand))
        d = {'transcripts': new_transcripts, 'exons': new_exons, 'introns': introns}
        return d


def gtf_for_ggplot(annotation, start, end, arrow_bins):
        arrow_space = int((end - start)/arrow_bins)
        s = """

        # data table with exons
        ann_list = list(
                "exons" = data.table(),
                "introns" = data.table()
        )
        """

        if annotation["exons"]:

                s += """
                ann_list[['exons']] = data.table(
                        tx = rep(c(%(tx_exons)s), c(%(n_exons)s)),
                        start = c(%(exon_start)s),
                        end = c(%(exon_end)s),
                        strand = c(%(strand)s)
                )
                """ %({
                "tx_exons": ",".join(annotation["exons"].keys()),
                "n_exons": ",".join(map(str, map(len, annotation["exons"].values()))),
                "exon_start" : ",".join(map(str, (v[0] for vs in annotation["exons"].values() for v in vs))),
                "exon_end" : ",".join(map(str, (v[1] for vs in annotation["exons"].values() for v in vs))),
                "strand" : ",".join(map(str, (v[2] for vs in annotation["exons"].values() for v in vs))),
                })

        if annotation["introns"]:

                s += """
                ann_list[['introns']] = data.table(
                        tx = rep(c(%(tx_introns)s), c(%(n_introns)s)),
                        start = c(%(intron_start)s),
                        end = c(%(intron_end)s),
                        strand = c(%(strand)s)
                )
                # Create data table for strand arrows
                txarrows = data.table()
                introns = ann_list[['introns']]
                # Add right-pointing arrows for plus strand
                if ("+" %%in%% introns$strand && nrow(introns[strand=="+" & end-start>5, ]) > 0) {
                        txarrows = rbind(
                                txarrows,
                                introns[strand=="+" & end-start>5, list(
                                        seq(start+4,end,by=%(arrow_space)s)-1,
                                        seq(start+4,end,by=%(arrow_space)s)
                                        ), by=.(tx,start,end)
                                ]
                        )
                }
                # Add left-pointing arrows for minus strand
                if ("-" %%in%% introns$strand && nrow(introns[strand=="-" & end-start>5, ]) > 0) {
                        txarrows = rbind(
                                txarrows,
                                introns[strand=="-" & end-start>5, list(
                                        seq(start,max(start+1, end-4), by=%(arrow_space)s),
                                        seq(start,max(start+1, end-4), by=%(arrow_space)s)-1
                                        ), by=.(tx,start,end)
                                ]
                        )
                }
                """ %({
                        "tx_introns": ",".join(annotation["introns"].keys()),
                        "n_introns": ",".join(map(str, map(len, annotation["introns"].values()))),
                        "intron_start" : ",".join(map(str, (v[0] for vs in annotation["introns"].values() for v in vs))),
                        "intron_end" : ",".join(map(str, (v[1] for vs in annotation["introns"].values() for v in vs))),
                        "strand" : ",".join(map(str, (v[2] for vs in annotation["introns"].values() for v in vs))),
                        "arrow_space" : arrow_space,
                })

        s += """

        gtfp = ggplot()
        if (length(ann_list[['introns']]) > 0) {
                gtfp = gtfp + geom_segment(data=ann_list[['introns']], aes(x=start, xend=end, y=tx, yend=tx), size=0.3)
                gtfp = gtfp + geom_segment(data=txarrows, aes(x=V1,xend=V2,y=tx,yend=tx), arrow=arrow(length=unit(0.02,"npc")))
        }
        if (length(ann_list[['exons']]) > 0) {
                gtfp = gtfp + geom_segment(data=ann_list[['exons']], aes(x=start, xend=end, y=tx, yend=tx), size=5, alpha=1)
        }
        gtfp = gtfp + scale_y_discrete(expand=c(0,0.5))
        gtfp = gtfp + scale_x_continuous(expand=c(0,0.25))
        gtfp = gtfp + coord_cartesian(xlim = c(%s,%s))
        gtfp = gtfp + labs(y=NULL)
        gtfp = gtfp + theme(axis.line = element_blank(), axis.text.x = element_blank(), axis.ticks = element_blank())
        """ %(start, end)

        return s


def setup_R_script(h, w, b, label_dict):
        s = """
        library(ggplot2)
        library(grid)
        library(gridExtra)
        library(data.table)
        library(gtable)

        scale_lwd = function(r) {
                lmin = 0.1
                lmax = 4
                return( r*(lmax-lmin)+lmin )
        }

        base_size = %(b)s
        height = ( %(h)s + base_size*0.352777778/67 ) * 1.02
        width = %(w)s
        theme_set(theme_bw(base_size=base_size))
        theme_update(
                plot.margin = unit(c(15,15,15,15), "pt"),
                panel.grid = element_blank(),
                panel.border = element_blank(),
                axis.line = element_line(size=0.5),
                axis.title.x = element_blank(),
                axis.title.y = element_text(angle=0, vjust=0.5)
        )

        labels = list(%(labels)s)

        density_list = list()
        junction_list = list()

        """ %({
                'h': h,
                'w': w,
                'b': b,
                'labels': ",".join(('"%s"="%s"' %(id,lab) for id,lab in label_dict.items())),
        })
        return s

def median(lst):
    quotient, remainder = divmod(len(lst), 2)
    if remainder:
        return sorted(lst)[quotient]
    return sum(sorted(lst)[quotient - 1:quotient + 1]) / 2.

def mean(lst):
        return sum(lst)/len(lst)


def make_R_lists(id_list, d, overlay_dict, aggr, intersected_introns):
        s = ""
        aggr_f = {
                "mean": mean,
                "median": median,
        }
        id_list = id_list if not overlay_dict else overlay_dict.keys()
        # Iterate over ids to get bam signal and junctions
        shrinked_introns = dict()
        for k in id_list:
                shrinked_introns_k, shrinked_intronsid = dict(), dict()
                x, y, dons, accs, yd, ya, counts = [], [], [], [], [], [], []
                if not overlay_dict:
                        x, y, dons, accs, yd, ya, counts = d[k]
                        if intersected_introns:
                                x, y = shrink_density(x, y, intersected_introns)
                                shrinked_introns_k, dons, accs = shrink_junctions(dons, accs, intersected_introns)
                                shrinked_introns.update(shrinked_introns_k)
                else:
                        for id in overlay_dict[k]:
                                xid, yid, donsid, accsid, ydid, yaid, countsid = d[id]
                                if intersected_introns:
                                        xid, yid = shrink_density(xid, yid, intersected_introns)
                                        shrinked_intronsid, donsid, accsid = shrink_junctions(donsid, accsid, intersected_introns)
                                        shrinked_introns.update(shrinked_intronsid)
                                x += xid
                                y += yid
                                dons += donsid
                                accs += accsid
                                yd += ydid
                                ya += yaid
                                counts += countsid
                        
                        
                        if aggr and "_j" not in aggr:
                                x = d[overlay_dict[k][0]][0]
                                y = list(map(aggr_f[aggr], zip(*(d[id][1] for id in overlay_dict[k]))))
                                if intersected_introns:
                                        x, y = shrink_density(x, y, intersected_introns)
                        #dons, accs, yd, ya, counts = [], [], [], [], []
                s += """
                density_list[["%(id)s"]] = data.frame(x=c(%(x)s), y=c(%(y)s))
                junction_list[["%(id)s"]] = data.frame(x=c(%(dons)s), xend=c(%(accs)s), y=c(%(yd)s), yend=c(%(ya)s), count=c(%(counts)s))
                """ %({
                        'id': k,
                        'x' : ",".join(map(str, x)),
                        'y' : ",".join(map(str, y)),
                        'dons' : ",".join(map(str, dons)),
                        'accs' : ",".join(map(str, accs)),
                        'yd' : ",".join(map(str, yd)),
                        'ya' : ",".join(map(str, ya)),
                        'counts' : ",".join(map(str, counts))
                })
        if intersected_introns:
                s+= """
                coord_dict = data.frame(shrinked=c(%(shrinked_introns_keys)s), real=c(%(shrinked_introns_values)s))
                intersected_introns = data.frame(real_x=c(%(intersected_introns_x)s), real_xend=c(%(intersected_introns_xend)s))
                """ %({
                        'shrinked_introns_keys': ','.join(map(str, shrinked_introns.keys())),
                        'shrinked_introns_values': ','.join(map(str, shrinked_introns.values())),
                        'intersected_introns_x': ','.join([str(coord[0]) for coord in intersected_introns]),
                        'intersected_introns_xend': ','.join([str(coord[1]) for coord in intersected_introns])
                })
        return s


def plot(R_script):
        p = sp.Popen("R --vanilla --slave", shell=True, stdin=sp.PIPE)
        p.communicate(input=R_script.encode('utf-8'))
        p.stdin.close()
        p.wait()
        return


def colorize(d, p):
        levels = list(OrderedDict.fromkeys(d.values()).keys())
        n = len(levels)
        if n > len(p):
                p = (p*n)[:n]

        s = "color_list = list(%s)\n" %( ",".join('"%s"="%s"' %(k, p[levels.index(v)]) for k,v in d.items()) )

        return s

def get_debug_info():
        """
        Return useful debug information:
        - OS info
        - Linux distribution
        - Python version
        - program version
        - output of R `sessionInfo()`
        """
        # get system info
        import platform
        system = platform.system()
        info = OrderedDict()
        info["OS"] = "{}-{}".format(system, platform.machine())
        if system == "Linux":
                release = sp.check_output(["lsb_release", "-ds"])
                info["Distro"] = release.strip().decode('utf-8')
        info["Python"] = platform.python_version()
        print('')
        maxlen = max(map(len, info.keys()))
        for k,v in info.items():
                print("{:{width}}: {:>}".format(k ,v, width=maxlen))
        print('')

        # get R session info
        r_exec = ';'.join([
                "library(ggplot2)",
                "library(grid)",
                "library(gridExtra)",
                "library(data.table)",
                "library(gtable)",
                "sessionInfo()",
        ])

        r_command = "R --vanilla --slave -e '{}'".format(r_exec)
        r_info = sp.check_output(r_command, shell=True, stderr=sp.STDOUT)
        print(r_info.strip().decode('utf-8'))




def define_options():

        class DebugInfoAction(ArgParseAction):

                def __init__(self, option_strings, dest, **kwargs):
                        super(DebugInfoAction, self).__init__(option_strings, dest, nargs=0, **kwargs)

                def __call__(self, parser, namespace, values, option_string=None):
                        returncode = 0
                        try:
                                get_debug_info()
                        except sp.CalledProcessError as CPE:
                                print("ERROR: {}".format(CPE.output.strip().decode('utf-8')))
                                returncode = CPE.returncode
                        except Exception as e:
                                print("ERROR: {}".format(e))
                                returncode = 1
                        finally:
                                parser.exit(returncode)

        """
        # Argument parsing
        parser = ArgumentParser(description='Create sashimi plot for a given genomic region')
        # parser.register('action', 'debuginfo', DebugInfoAction)
        parser.add_argument("--intron", type=str, required=True,
                help="
                The intron count file by tealeaf
                ")
        
        parser.add_argument("--exon", type=str, required=True,
                help="
                The exon count file by tealeaf
                ")
        
        parser.add_argument("--strand_info", type=str, required=True,
                help="
                The intron exon connectivity file that contain information for strandness
                ")
        parser.add_argument("-c", "--coordinates", type=str, required=True,
                help="Genomic region. Format: chr:start-end. Remember that bam coordinates are 0-based")
        
        """
        

        # Argument parsing
        parser = ArgumentParser(description='Create sashimi plot for a given genomic region')
        # parser.register('action', 'debuginfo', DebugInfoAction)
        parser.add_argument("--intron", dest = "intron", default = None,
                help="""The intron count file by tealeaf""")
        
        parser.add_argument("--exon", dest= "exon",  default = None,
                help="""
                The exon count file by tealeaf
                """)
        parser.add_argument("--strand_info", dest = "strand_info", type=str, default = "/gpfs/commons/home/xzhang/reference/vM32_leafcutterITI/vM32_intron_exon_connectivity.tsv",
                help="""
                The intron exon connectivity file that contain information for strandness
                """)
  
        parser.add_argument("--data_type", dest = "data_type", type=str, default = "bulk",
                help="The data type that input into tealeaf, could be bulk or sc, use to determine the sample name (default: bulk)")
                    
        parser.add_argument("--normalized", dest = "normalized",type=str, default = "True",
                help="")
        
        parser.add_argument("--aggregation", dest = "aggregation",type=str, default = "True",
                help="Whether to aggregate the value from tealeaf by mean, should use either this parameter or -A paramether for ggsashimi")
                   
        parser.add_argument("-c", "--coordinates", type=str,dest = "coordinates", required = True,
                help="Genomic region. Format: chr:start-end.")
        parser.add_argument("-g", "--gtf", dest = "gtf", required= True,
                help="Gtf file with annotation (only exons is enough)")
        
        
        
        
        
        
        parser.add_argument("-o", "--out-prefix", type=str, dest="out_prefix", default="sashimi",
                help="Prefix for plot file name [default=%(default)s]")
        parser.add_argument("-S", "--out-strand", type=str, dest="out_strand", default="both",
                help="Only for --strand other than 'NONE'. Choose which signal strand to plot: <both> <plus> <minus> [default=%(default)s]")
        parser.add_argument("-M", "--min-coverage", type=int, default=1, dest="min_coverage",
                help="Minimum number of reads supporting a junction to be drawn [default=1]")
        parser.add_argument("-j", "--junctions-bed", type=str, dest = "junctions_bed", default="",
                help="Junction BED file name [default=no junction file]")
        parser.add_argument("-s", "--strand", default="NONE", type=str,
                help="Strand specificity: <NONE> <SENSE> <ANTISENSE> <MATE1_SENSE> <MATE2_SENSE> [default=%(default)s]")
        parser.add_argument("--shrink", action="store_true",
                help="Shrink the junctions by a factor for nicer display [default=%(default)s]")
        
        
        
        parser.add_argument("-A", "--aggr", type=str, default="",
                help="""Aggregate function for overlay: <mean> <median> <mean_j> <median_j>.
                        Use mean_j | median_j to keep density overlay but aggregate junction counts [default=no aggregation]""")
                        
                        
        parser.add_argument("--alpha", type=float, default=0.5,
                help="Transparency level for density histogram [default=%(default)s]")
        
        parser.add_argument("-P", "--palette", type=str,
                help="Color palette file. tsv file with >=1 columns, where the color is the first column. Both R color names and hexadecimal values are valid")
        

        
        parser.add_argument("--fix-y-scale", default=False, action="store_true", dest = "fix_y_scale",
                help="Fix y-scale across individual signal plots [default=%(default)s]")
        parser.add_argument("--height", type=float, default=2,
                help="Height of the individual signal plot in inches [default=%(default)s]")
        parser.add_argument("--ann-height", type=float, default=1.5, dest="ann_height",
                help="Height of annotation plot in inches [default=%(default)s]")
        parser.add_argument("--width", type=float, default=10,
                help="Width of the plot in inches [default=%(default)s]")
        parser.add_argument("--base-size", type=float, default=14, dest="base_size",
                help="Base font size of the plot in pch [default=%(default)s]")
        parser.add_argument("-F", "--out-format", type=str, default="pdf", dest="out_format",
                help="Output file format: <pdf> <svg> <png> <jpeg> <tiff> [default=%(default)s]")
        parser.add_argument("-R", "--out-resolution", type=int, default=300, dest="out_resolution",
                help="Output file resolution in PPI (pixels per inch). Applies only to raster output formats [default=%(default)s]")
        parser.add_argument("--debug-info", action=DebugInfoAction,
                help="Show several system information useful for debugging purposes [default=%(default)s]")
        #parser.add_argument("-s", "--smooth", action="store_true", default=False, help="Smooth the signal histogram")
        return parser


        
        
        
        
if __name__ == "__main__":

        strand_dict = {"plus": "+", "minus": "-"}

        parser = define_options()
        """
        if len(sys.argv)==1:
            parser.print_help()
            sys.exit(1)
        """
        args = parser.parse_args()

#       args.coordinates = "chrX:9609491-9612406"
#       args.coordinates = "chrX:9609491-9610000"
#       args.bam = "/nfs/no_backup/rg/epalumbo/projects/tg/work/8b/8b0ac8705f37fd772a06ab7db89f6b/2A_m4_n10_toGenome.bam"
                
        
        if args.exon == None or args.intron == None:
            print("ERROR: No intron or exon file is provided.")
            exit(1)
                
                
                

        palette = read_palette(args.palette)

        sample_dict, overlay_dict, color_dict, id_list, label_dict = {"+":OrderedDict()}, OrderedDict(), OrderedDict(), [], OrderedDict()
        if args.strand != "NONE": sample_dict["-"] = OrderedDict()
        if args.junctions_bed != "": junctions_list = []

        
        df_intron = pd.read_csv(args.intron , sep = ' ')
        df_exon = pd.read_csv(args.exon, sep = ' ')
        strand_dic = gen_intron_exon_stand_dic(args.strand_info)
        
        df_intron, sample_cell_type_dic = aggregate_cell_types(df_intron, args.data_type, args.normalized,args.aggregation)
        df_exon, sample_cell_type_dic = aggregate_cell_types(df_exon, args.data_type, args.normalized, args.aggregation)
        
        
        
        
                               
        for id, sample, overlay_level, color_level, label_text in read_intron_exon_input(df_intron, df_exon, sample_cell_type_dic):
                # overlay, color, label need to be fix 
 
                a, junctions = read_intron_exon(df_intron, df_exon, id, strand_dic, args.coordinates, args.strand)
                
                print(f'finish a and junction for {sample}')
                
                if a.keys() == ["+"] and all(map(lambda x: x==0, list(a.values()[0]))):
                        print("WARN: Sample {} has no reads in the specified area.".format(id))
                        continue
                id_list.append(id)
                label_dict[id] = label_text
                for strand in a:
                        # Store junction information
                        if args.strand == "NONE" or args.out_strand == 'both' or strand == strand_dict[args.out_strand]:
                                if args.junctions_bed:
                                        for k,v in zip(junctions[strand].keys(), junctions[strand].values()):
                                                if v >= args.min_coverage:
                                                        junctions_list.append('\t'.join([args.coordinates.split(':')[0], str(k[0]), str(k[1]), id, str(v), strand]))
                        sample_dict[strand][id] = prepare_for_R(a[strand], junctions[strand], args.coordinates, args.min_coverage)
                        
                
                
                
                if color_level is None:
                        color_dict.setdefault(id, id)
                if overlay_level is not None:
                        overlay_dict.setdefault(overlay_level, []).append(id)
                        label_dict[overlay_level] = overlay_level
                        color_dict.setdefault(overlay_level, overlay_level)
                if overlay_level is None:
                        color_dict.setdefault(id, color_level)




        # No intron and exon files
        if not sample_dict["+"]:
                print("ERROR: No available intron and exon files.")
                exit(1)

        # Write junctions to BED
        if args.junctions_bed:
                if not args.junctions_bed.endswith('.bed'):
                        args.junctions_bed = args.junctions_bed + '.bed'
                jbed = open(args.junctions_bed, 'w')
                jbed.write('\n'.join(sorted(junctions_list)))
                jbed.close()

        if args.gtf:
                transcripts, exons = read_gtf(args.gtf, args.coordinates)

        if args.out_format not in ('pdf', 'png', 'svg', 'tiff', 'jpeg'):
                print("ERROR: Provided output format '%s' is not available. Please select among 'pdf', 'png', 'svg', 'tiff' or 'jpeg'" % args.out_format)
                exit(1)

        # Iterate for plus and minus strand
        for strand in sample_dict:

                # Output file name (allow tiff/tif and jpeg/jpg extensions)
                if args.out_prefix.endswith(('.pdf', '.png', '.svg', '.tiff', '.tif', '.jpeg', '.jpg')):
                        out_split = os.path.splitext(args.out_prefix)
                        if (args.out_format == out_split[1][1:] or
                        args.out_format == 'tiff' and out_split[1] in ('.tiff','.tif') or
                        args.out_format == 'jpeg' and out_split[1] in ('.jpeg','.jpg')):
                                args.out_prefix = out_split[0]
                                out_suffix = out_split[1][1:]
                        else:
                                out_suffix = args.out_format
                else:
                        out_suffix = args.out_format
                out_prefix = args.out_prefix + "_" + strand
                if args.strand == "NONE":
                        out_prefix = args.out_prefix
                else:
                        if args.out_strand != "both" and strand != strand_dict[args.out_strand]:
                                continue

                # Find set of junctions to perform shrink
                intersected_introns = None
                if args.shrink:
                        introns = (v for vs in sample_dict[strand].values() for v in zip(vs[2], vs[3]))
                        intersected_introns = list(intersect_introns(introns))


                # *** PLOT *** Define plot height
                bam_height = args.height * len(id_list)

                bam_height = args.height * len(overlay_dict)
                if args.gtf:
                        bam_height += args.ann_height

                # *** PLOT *** Start R script by loading libraries, initializing variables, etc...
                R_script = setup_R_script(bam_height, args.width, args.base_size, label_dict)

                R_script += colorize(color_dict, palette)

                # *** PLOT *** Prepare annotation plot only for the first bam file
                arrow_bins = 50
                if args.gtf:
                        # Make introns from annotation (they are shrunk if required)
                        annotation = make_introns(transcripts, exons, intersected_introns)
                        x = list(sample_dict[strand].values())[0][0]
                        if args.shrink:
                                x, _ = shrink_density(x, x, intersected_introns)
                        R_script += gtf_for_ggplot(annotation, x[0], x[-1], arrow_bins)

                R_script += make_R_lists(id_list, sample_dict[strand], overlay_dict, args.aggr, intersected_introns)

                R_script += """

                pdf(NULL) # just to remove the blank pdf produced by ggplotGrob

                if(packageVersion('ggplot2') >= '3.0.0'){  # fix problems with ggplot2 vs >3.0.0
                        vs = 1
                } else {
                        vs = 0
                }

                if(%(fix_y_scale)s) {

                        maxheight = max(unlist(lapply(density_list, function(df){max(df$y)})))
                        breaks_y = labeling::extended(0, maxheight, m = 4)
                        maxheight_j = 0
                        for (bam_index in 1:length(density_list)) { # Emulate code below to obtain max top arch height

                                id = names(density_list)[bam_index]
                                d = data.table(density_list[[id]])
                                junctions = data.table(junction_list[[id]])

                                row_i = c()
                                if (nrow(junctions) >0 ) {
                                        junctions$jlabel = as.character(junctions$count)
                                        junctions = setNames(junctions[,.(max(y), max(yend),round(%(args.aggr)s(count)),paste(jlabel,collapse=",")), keyby=.(x,xend)], names(junctions))
                                        if ("%(args.aggr)s" != "") {
                                                junctions = setNames(junctions[,.(max(y), max(yend),round(%(args.aggr)s(count)),round(%(args.aggr)s(count))), keyby=.(x,xend)], names(junctions))
                                        }
                                        row_i = 1:nrow(junctions)
                                }
                                for (i in row_i) {
                                        j = as.numeric(junctions[i,1:5])
                                        if ("%(args.aggr)s" != "") {
                                                j[3] = ifelse(length(d[x==j[1]-1,y])==0, 0, max(as.numeric(d[x==j[1]-1,y])))
                                                j[4] = ifelse(length(d[x==j[2]+1,y])==0, 0, max(as.numeric(d[x==j[2]+1,y])))
                                        }
                                        if (i%%%%2 != 0) { #top
                                                set.seed(mean(j[3:4]))
                                                maxheight_j = max(maxheight_j, max(j[3:4]) * runif(1, 1.2, 1.5))
                                        }
                                }
                        }
                }

                if(exists('coord_dict')){
                        all_pos_shrinked = do.call(rbind, density_list)$x
                        s2r = merge(intersected_introns, coord_dict, by.x = 'real_xend', by.y = 'real')
                        s2r = merge(s2r, coord_dict, by.x = 'real_x', by.y = 'real', suffixes = c('_xend', '_x'))
                        breaks_x_shrinked = labeling::extended(min(all_pos_shrinked), max(all_pos_shrinked), m = 5)
                        breaks_x = c()
                        out_range = c()
                        for (b in breaks_x_shrinked){
                                iintron = FALSE
                                for (j in 1:nrow(s2r)){
                                        l = s2r[j, ]
                                        if(b >= l$shrinked_x && b <= l$shrinked_xend){
                                                # Intersected intron
                                                p = (b-l$shrinked_x)/(l$shrinked_xend - l$shrinked_x)
                                                realb = round(l$real_x + p*(l$real_xend - l$real_x))
                                                breaks_x = c(breaks_x, realb)
                                                iintron = TRUE
                                                break
                                        }
                                }
                                if (!iintron){
                                        # Exon, upstream/downstream intergenic region or intron (not intersected)
                                        if(b <= min(s2r$shrinked_x)) {
                                                l <- s2r[which.min(s2r$shrinked_x), ]
                                                if(any(b == all_pos_shrinked)){
                                                        # Boundary (subtract)
                                                        s = l$shrinked_x - b
                                                        realb = l$real_x - s
                                                        breaks_x = c(breaks_x, realb)
                                                } else {
                                                        out_range <- c(out_range, which(breaks_x_shrinked == b))
                                                }
                                        } else if (b >= max(s2r$shrinked_xend)){
                                                l <- s2r[which.max(s2r$shrinked_xend), ]
                                                if(any(b == all_pos_shrinked)){
                                                        # Boundary (sum)
                                                        s = b - l$shrinked_xend
                                                        realb = l$real_xend + s
                                                        breaks_x = c(breaks_x, realb)
                                                } else {
                                                        out_range <- c(out_range, which(breaks_x_shrinked == b))
                                                }
                                        } else {
                                                delta = b-s2r$shrinked_xend
                                                delta[delta < 0] = Inf
                                                l = s2r[which.min(delta), ]
                                                # Internal (sum)
                                                s = b - l$shrinked_xend
                                                realb = l$real_xend + s
                                                breaks_x = c(breaks_x, realb)
                                        }
                                }
                        }
                        if(length(out_range)) {
                                breaks_x_shrinked = breaks_x_shrinked[-out_range]
                        }
                }

                density_grobs = list();

                for (bam_index in 1:length(density_list)) {

                        id = names(density_list)[bam_index]
                        d = data.table(density_list[[id]])
                        junctions = data.table(junction_list[[id]])

                        # Density plot
                        gp = ggplot(d) + geom_bar(aes(x, y), width=1, position='identity', stat='identity', fill=color_list[[id]], alpha=%(alpha)s)
                        gp = gp + labs(y=labels[[id]])
                        if(exists('coord_dict')) {
                                gp = gp + scale_x_continuous(expand=c(0, 0.25), breaks = breaks_x_shrinked, labels = breaks_x)
                        } else {
                                gp = gp + scale_x_continuous(expand=c(0, 0.25))
                        }

                        if(!%(fix_y_scale)s){
                                maxheight = max(d[['y']])
                                breaks_y = labeling::extended(0, maxheight, m = 4)
                                gp = gp + scale_y_continuous(breaks = breaks_y)
                        } else {
                                gp = gp + scale_y_continuous(breaks = breaks_y, limits = c(NA, max(maxheight, maxheight_j)))
                        }

                        # Aggregate junction counts
                        row_i = c()
                        if (nrow(junctions) >0 ) {

                                junctions$jlabel = as.character(junctions$count)
                                junctions = setNames(junctions[,.(max(y), max(yend),round(%(args.aggr)s(count)),paste(jlabel,collapse=",")), keyby=.(x,xend)], names(junctions))
                                if ("%(args.aggr)s" != "") {
                                        junctions = setNames(junctions[,.(max(y), max(yend),round(%(args.aggr)s(count)),round(%(args.aggr)s(count))), keyby=.(x,xend)], names(junctions))
                                }
                                # The number of rows (unique junctions per bam) has to be calculated after aggregation
                                row_i = 1:nrow(junctions)
                        }


                        for (i in row_i) {

                                j_tot_counts = sum(junctions[['count']])

                                j = as.numeric(junctions[i,1:5])

                                if ("%(args.aggr)s" != "") {
                                        j[3] = ifelse(length(d[x==j[1]-1,y])==0, 0, max(as.numeric(d[x==j[1]-1,y])))
                                        j[4] = ifelse(length(d[x==j[2]+1,y])==0, 0, max(as.numeric(d[x==j[2]+1,y])))
                                }

                                # Find intron midpoint
                                xmid = round(mean(j[1:2]), 1)
                                set.seed(mean(j[3:4]))
                                ymid = max(j[3:4]) * runif(1, 1.2, 1.5)

                                # Thickness of the arch
                                lwd = scale_lwd(j[5]/j_tot_counts)

                                curve_par = gpar(lwd=lwd, col=color_list[[id]])

                                # Arc grobs

                                # Choose position of the arch (top or bottom)
                                nss = i
                                if (nss%%%%2 == 0) {  #bottom
                                        ymid = -runif(1, 0.2, 0.4) * maxheight
                                        # Draw the arcs
                                        # Left
                                        curve = xsplineGrob(x=c(0, 0, 1, 1), y=c(1, 0, 0, 0), shape=1, gp=curve_par)
                                        gp = gp + annotation_custom(grob = curve, j[1], xmid, 0, ymid)
                                        # Right
                                        curve = xsplineGrob(x=c(1, 1, 0, 0), y=c(1, 0, 0, 0), shape=1, gp=curve_par)
                                        gp = gp + annotation_custom(grob = curve, xmid, j[2], 0, ymid)
                                }

                                if (nss%%%%2 != 0) {  #top
                                        # Draw the arcs
                                        # Left
                                        curve = xsplineGrob(x=c(0, 0, 1, 1), y=c(0, 1, 1, 1), shape=1, gp=curve_par)
                                        gp = gp + annotation_custom(grob = curve, j[1], xmid, j[3], ymid)
                                        # Right
                                        curve = xsplineGrob(x=c(1, 1, 0, 0), y=c(0, 1, 1, 1), shape=1, gp=curve_par)
                                        gp = gp + annotation_custom(grob = curve, xmid, j[2], j[4], ymid)
                                }

                                # Add junction labels
                                gp = gp + annotate("label", x = xmid, y = ymid, label = as.character(junctions[i,6]),
                                        vjust=0.5, hjust=0.5, label.padding=unit(0.01, "lines"),
                                        label.size=NA, size=(base_size*0.352777778)*0.6
                                )


                #               gp = gp + annotation_custom(grob = rectGrob(x=0, y=0, gp=gpar(col="red"), just=c("left","bottom")), xmid, j[2], j[4], ymid)
                #               gp = gp + annotation_custom(grob = rectGrob(x=0, y=0, gp=gpar(col="green"), just=c("left","bottom")), j[1], xmid, j[3], ymid)


                        }

                        gpGrob = ggplotGrob(gp);
                        gpGrob$layout$clip[gpGrob$layout$name=="panel"] <- "off"
                        if (bam_index == 1) {
                                maxWidth = gpGrob$widths[2+vs] + gpGrob$widths[3+vs];    # fix problems ggplot2 vs
                                maxYtextWidth = gpGrob$widths[3+vs];                     # fix problems ggplot2 vs
                                # Extract x axis grob (trim=F --> keep empty cells)
                                xaxisGrob <- gtable_filter(gpGrob, "axis-b", trim=F)
                                xaxisGrob$heights[8+vs] = gpGrob$heights[1]              # fix problems ggplot2 vs
                                x.axis.height = gpGrob$heights[7+vs] + gpGrob$heights[1] # fix problems ggplot2 vs
                        }


                        # Remove x axis from all density plots
                        kept_names = gpGrob$layout$name[gpGrob$layout$name != "axis-b"]
                        gpGrob <- gtable_filter(gpGrob, paste(kept_names, sep="", collapse="|"), trim=F)

                        # Find max width of y text and y label and max width of y text
                        maxWidth = grid::unit.pmax(maxWidth, gpGrob$widths[2+vs] + gpGrob$widths[3+vs]); # fix problems ggplot2 vs
                        maxYtextWidth = grid::unit.pmax(maxYtextWidth, gpGrob$widths[3+vs]); # fix problems ggplot2 vs
                        density_grobs[[id]] = gpGrob;
                }

                # Add x axis grob after density grobs BEFORE annotation grob
                density_grobs[["xaxis"]] = xaxisGrob

                # Annotation grob
                if (%(args.gtf)s == 1) {
                        gtfGrob = ggplotGrob(gtfp);
                        maxWidth = grid::unit.pmax(maxWidth, gtfGrob$widths[2+vs] + gtfGrob$widths[3+vs]); # fix problems ggplot2 vs
                        density_grobs[['gtf']] = gtfGrob;
                }

                # Reassign grob widths to align the plots
                for (id in names(density_grobs)) {
                        density_grobs[[id]]$widths[1] <- density_grobs[[id]]$widths[1] + maxWidth - (density_grobs[[id]]$widths[2+vs] + maxYtextWidth); # fix problems ggplot2 vs
                        density_grobs[[id]]$widths[3+vs] <- maxYtextWidth # fix problems ggplot2 vs
                }

                # Heights for density, x axis and annotation
                heights = unit.c(
                        unit(rep(%(signal_height)s, length(density_list)), "in"),
                        x.axis.height,
                        unit(%(ann_height)s*%(args.gtf)s, "in")
                        )

                # Arrange grobs
                argrobs = arrangeGrob(
                        grobs=density_grobs,
                        ncol=1,
                        heights = heights,
                );

                # Save plot to file in the requested format
                if ("%(out_format)s" == "tiff"){
                        # TIFF images will be lzw-compressed
                        ggsave("%(out)s", plot = argrobs, device = "tiff", width = width, height = height, units = "in", dpi = %(out_resolution)s, compression = "lzw", limitsize = FALSE)
                } else {
                        ggsave("%(out)s", plot = argrobs, device = "%(out_format)s", width = width, height = height, units = "in", dpi = %(out_resolution)s, limitsize = FALSE)
                }

                dev.log = dev.off()

                """ %({
                        "out": "%s.%s" % (out_prefix, out_suffix),
                        "out_format": args.out_format,
                        "out_resolution": args.out_resolution,
                        "args.gtf": float(bool(args.gtf)),
                        "args.aggr": args.aggr.rstrip("_j"),
                        "signal_height": args.height,
                        "ann_height": args.ann_height,
                        "alpha": args.alpha,
                        "fix_y_scale": ("TRUE" if args.fix_y_scale else "FALSE")
                        })
                if os.getenv('GGSASHIMI_DEBUG') is not None:
                        with open("R_script", 'w') as r:
                                r.write(R_script)
                else:
                        plot(R_script)
        exit()


