import seqhandlers as sq
import bamnostic as bs
import collections

#cigar tags to exclude to get reference coordinates
ref_exclude = [4, 1, 5, 6, 8]

#Summarizing the data
def create_full_dic(read_info, var, opposite):
    """Put together info of all reads"""
    
    readnames = list(set([read[0] for read in read_info]))
    fullreaddic = {readname : [] for readname in readnames}
    
    for read in read_info:
        ref_snps = sq.snp_coordinates(read, ref_exclude)
        nonsnps = sq.nonsnp_coordinates(var, read, ref_snps)
            
        readdic = sq.opposite_snps(opposite,
                                sq.create_readdic(var, read, ref_snps, nonsnps))
    
        if len(fullreaddic[read[0]]) != 0:
            fullreaddic[read[0]] = fullreaddic[read[0]] | readdic
        else:
            fullreaddic[read[0]] = readdic

    return(fullreaddic)

#Filtering
def filter_readdic(fullreaddic, minstretch, minscore):
    """Filter interesting reads from all"""
    filtereddic = collections.defaultdict(lambda : [])
    
    for read, readdic in fullreaddic.items():
        ps_cigar = sq.read_summary(list(readdic.values()))
        res = sq.rec_possible(ps_cigar, minstretch, minscore)

        if res[0] == 'True':
            filtereddic[read] = [readdic, ps_cigar, res]

    return(filtereddic)

#Data for plotting
def plot_data(fullreaddic):
    """Create tsv that makes plotting easy"""
    data = [[k, a, b]
            for k, v in fullreaddic.items()
           for a, b in v.items()]
    
    return(data)

#Output
def dicprint(filtdic):
    f.open()
    out.write(f'{read}\t' + ''.join(ps_cigar) + '\t' + '\t'.join(res) + '\n')