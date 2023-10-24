import collections
import bamnostic as bs

#Data
#Generic, works for every sam file
#cigar tags to exclude to get reference coordinates
ref_exclude = [4, 1, 5, 6, 8] 
#cigar tags to exclude to get query coordinates
q_exclude = [2, 3, 4, 5, 6, 8]

#Functions

#Chr, pos, ref, alt from vcf
def get_variants(vcf_file):
    """Read in variant coordinates from vcf"""
    
    with open(vcf_file) as vcf:
        #Get chr, id, pos, ref, alt
        variant_info = [line.strip().split('\t')[0:5]
                        for line in vcf
                        if not line.startswith('#')] 
    
        #HAve to add indel-control, or just filter them from vcf
    variant_info = [var for var in variant_info
                    if len(var[3]) == 1]
                
    return(variant_info)


#Figure out the variants
def snp_coordinates(read, exclusion):
    """Function to get sequence coordinates"""
    """Use predefined cigar tag list as exclusion"""
    
    rq_snps = []
    current_coord = read[1] #start
    
    for op in read[4]: #cigar
        if op[0] not in exclusion:
            current_coord += int(op[1])
        elif op[0] == 8:
            for i in range(0,op[1]):
                current_coord += 1
                rq_snps.append(current_coord)
                i += 1
                
    return(rq_snps)

def nonsnp_coordinates(var, read, rq_snps):
    """Get variant positions with no actual variant"""
    non_snps = [va[1] for va in var
                if int(va[1]) < read[2] and
                int(va[1]) > read[1]
                and int(va[1]) not in rq_snps]

    non_snps = non_snps
    return(non_snps)

def create_readdic(var, read, ref_snps, non_snps):
    """Generic sequence creator"""
    readseq = ""
    readcoords = ""
    
    readseq = [readseq + va[3] for va in var
              if int(va[1]) > read[1] #start
              and int(va[1]) < read[2]] #end
    readcoords = [readcoords + va[1] for va in var
                  if int(va[1]) > read[1]
                  and int(va[1]) < read[2]]
    read_dic = dict(zip(readcoords, readseq))

    for pos, va in read_dic.items():  #check if this works as intended
        if pos in non_snps:
            read_dic[pos] = "A"
        elif int(pos) in ref_snps:
            read_dic[pos] = "B"
        
    return(read_dic)

#Pseudocigar for pseudoseq
def read_summary(pseudoseq):
    i = 0
    result = []
    
    pseudoseq.append(object()) # append a dummy item first
    while(i < len(pseudoseq)-1):
        count = 1
        result.append(pseudoseq[i])
        while(pseudoseq[i] == pseudoseq[i+1]):
            count += 1
            i += 1
    
        result.append(str(count))
        i += 1

    pseudoseq.pop()
    return(result)

#Very basic check for recombination possibility
def rec_possible(ps_cigar, minstretch, minscore):
    """Scoring the possibility of recombination"""

    score = 0
    if len(ps_cigar) > 2:
        for i in ps_cigar[1::2]:
            #print(i)
            if int(i) >= minstretch:
                score += 1
            elif int(i) == 1:
                score -= 1
    if score >= minscore:
        return(["True", str(score)])
    else:
        return(["False", str(score)])
    

def opposite_snps(opposite, readdic):
    """Called snps that are actually in the ref type"""
    for pos, base in readdic.items():
        if pos in opposite:
            if readdic[pos] == "A":
                readdic[pos]= "B"
            else:
                readdic[pos] = "A"

    return(readdic)
