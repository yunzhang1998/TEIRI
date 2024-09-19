#!/usr/bin/python
import argparse
import sys
import re
import random
import string
import time
import math
import subprocess
from collections import Counter

parser = argparse.ArgumentParser(description='TEIRI_consolidate')
parser.add_argument('-i', '--gtf_list', help="A text file with a list of GTF files (required)")
parser.add_argument('-r', '--reference_gtf', help="Reference genome annotation in GTF format (required)")
parser.add_argument('--TE_anno', help="TE annotation in BED format (required)")
parser.add_argument('--tss_merge_distance', default=50, type=int, help='The min distance to merge TE-derived TSSs of different GTF files (default: 50)')
parser.add_argument('--min_exon_length', default=20, type=int, help='The min length of first exon (default: 20)')
parser.add_argument('-p', '--prefix', default='TEIRI', help="Prefix for output file (default: TEIRI)")


args = parser.parse_args()


def process_gtf(gtf_path):
    transcript = {}
    print((gtf_path+"_load"))
    with open(gtf_path) as file:
        for eachline in file:
            if not eachline.strip().startswith('#'):
                temp = eachline.strip().split("\t")
                type = temp[2]
                if type == "exon":
                    anno = temp[8].split("\"; ")
                    for transcriptinfo in anno:
                        if "transcript_id" in transcriptinfo:
                            transcriptid = transcriptinfo.split("\"")[1]
                            break
                    id = "\t".join([temp[0], str(temp[6]), transcriptid])
                    if id not in transcript:
                        transcript[id] = []
                    transcript[id].append(int(temp[3]))
                    transcript[id].append(int(temp[4]))
    return transcript

def tss_merge(transcript,tss_merge_distance,min_exon_length):
    merged_transcript = {}
    merged_tss = {}
    for gtf in transcript:
        for transcriptID in transcript[gtf]:
            id = transcriptID.split("\t")
            exons=transcript[gtf][transcriptID]
            num_exons=int(len(exons)/2)
            chr = "\t".join([id[0],id[1]])
            if id[1]=="+":
                exons.sort() 
            elif id[1]=="-":
                exons.sort(reverse=True)
            strand = id[1]
            tss = str(exons[0])
            tes = str(exons[-1])
            exons = list(map(str, exons))
            if num_exons >= 2: 
                junctions = "\t".join(exons[1:-1])
            else:
                junctions = "SE"
            if chr not in merged_transcript:
                merged_transcript[chr] = {}
                merged_tss [chr] = {}
            if  tss not in merged_transcript[chr]:
                merged_transcript[chr][tss]={}
                merged_tss[chr][tss]=[gtf]
            else:
                merged_tss[chr][tss].append(gtf)
            if junctions not in merged_transcript[chr][tss]:
                merged_transcript[chr][tss][junctions] = tes
            else:
                if strand == "+" and int(tes) < int(merged_transcript[chr][tss][junctions]):
                    merged_transcript[chr][tss][junctions] = tes
                elif strand == "-" and int(tes) > int(merged_transcript[chr][tss][junctions]):
                    merged_transcript[chr][tss][junctions] = tes

    for chr in merged_tss:
        strand= chr.split("\t")[1]
        chr_= chr.split("\t")[0]
        tss_list=list(map(int,list(merged_tss[chr])))
        if len(tss_list) >1 :
            if strand == "+" :
                tss_list.sort()
            elif strand == "-":
                tss_list.sort(reverse=True)
            tss_list=list(map(str,tss_list))
            tss_cluster = []
            sublist = [tss_list[0]]
            for i in range(1, len(tss_list)):
                if abs(int(tss_list[i]) - int(sublist[-1])) <= tss_merge_distance:
                    sublist.append(tss_list[i])
                else:
                    tss_cluster.append(sublist)
                    sublist = [tss_list[i]]
            tss_cluster.append(sublist)
            for cluster in tss_cluster:
                if len(cluster) >1 :
                    max_score = 0
                    for temp in cluster:
                        if len(list(set(merged_tss[chr][temp])))> max_score:
                            max_score = len(list(set(merged_tss[chr][temp])))
                            best_tss=temp
                    for temp_ in cluster:
                        if temp_!=best_tss:
                            del_state = 0
                            for junctions in merged_transcript[chr][temp_]:
                                tes=merged_transcript[chr][temp_][junctions]
                                if junctions not in merged_transcript[chr][best_tss] :
                                    if junctions == "SE":
                                        exon1_end = tes
                                    else:
                                        exon1_end = junctions.split("\t")[0]
                                    if strand == "+" and int(exon1_end) > (int(best_tss) + min_exon_length):
                                        merged_transcript[chr][best_tss][junctions] = tes  
                                    elif strand == "-" and (int(exon1_end) + min_exon_length) < int(best_tss) :
                                        merged_transcript[chr][best_tss][junctions] = tes 
                                    else:
                                        del_state = 1
                                else:
                                    if strand == "+" and int(tes) < int(merged_transcript[chr][best_tss][junctions]):
                                        merged_transcript[chr][best_tss][junctions] = tes
                                    elif strand == "-" and int(tes) > int(merged_transcript[chr][best_tss][junctions]):
                                        merged_transcript[chr][best_tss][junctions] = tes
                            merged_tss[chr][best_tss]=merged_tss[chr][best_tss] + merged_tss[chr][temp_]
                            if del_state==0:
                                del(merged_transcript[chr][temp_])      
                                del(merged_tss[chr][temp_])

    return merged_transcript,merged_tss          

def custom_sort_key(s):
    s = s.split("\t")[0]
    match = re.match(r'([a-zA-Z]+)(\d+)', s)
    if match:
        alpha, num = match.groups()
        return alpha, int(num)
    else:
        return s, 0


def write_exon(exons,gene_id,transcriptid,chr_,strand):
    exons = list(map(int, exons))
    num_exons=int(len(exons)/2)
    exons.sort()
    exons = list(map(str, exons))
    stringreturn=[chr_, "TE_RNA","transcript" , exons[0] ,exons[-1],".", strand,".",str("transcript_id \"teRNA_"+ str(gene_id) + "." + str(transcriptid) + "\"; gene_id \"TE_" + str(gene_id) + "\";" ) ]
    stringreturn = "\t".join([str(x) for x in stringreturn])
    yield stringreturn
    for i in range(0,num_exons):
        stringreturn=[chr_, "TE_RNA","exon" , exons[2*i] ,exons[2*i+1],".", strand,".",str("transcript_id \"teRNA_"+ str(gene_id) + "." + str(transcriptid) + "\"; gene_id \"TE_" + str(gene_id) + "\";" ) ]
        stringreturn = "\t".join([str(x) for x in stringreturn])
        yield stringreturn

def write_gtf(merged_transcript,merged_tss,output_file):
    
    gene_id=0

    with open(output_file, "w") as f_out:
        chr_list = list(merged_transcript)
        chr_list = sorted(chr_list, key=custom_sort_key)
        for chr in chr_list:
            tss_list=list(map(int,list(merged_transcript[chr])))
            tss_list.sort()
            tss_list = list(map(str, tss_list))
            strand= chr.split("\t")[1]
            chr_= chr.split("\t")[0]
            for tss in tss_list:
                transcriptid = 0
                gene_id += 1
                for junctions in merged_transcript[chr][tss]:
                    transcriptid += 1
                    tes = merged_transcript[chr][tss][junctions]
                    if junctions =="SE":
                        exons = [tss] + [tes]
                    else:
                        exons = [tss] + junctions.split("\t") + [tes]
                    for stringreturn in write_exon(exons,gene_id,transcriptid,chr_,strand):
                        f_out.write(stringreturn + "\n")


def write_bed(merged_transcript,merged_tss,output_file):
    gene_id=0

    with open(output_file, "w") as f_out:
        chr_list = list(merged_transcript)
        chr_list = sorted(chr_list, key=custom_sort_key)
        for chr in chr_list:
            tss_list=list(map(int,list(merged_transcript[chr])))
            tss_list.sort()
            tss_list = list(map(str, tss_list))
            strand= chr.split("\t")[1]
            chr_= chr.split("\t")[0]
            for tss in tss_list:
                gene_id += 1
                tissue_list = list(set(merged_tss[chr][tss]))
                for tissue in tissue_list:
                    stringreturn=[chr_, int(tss)-1, int(tss) , "TE_"+ str(gene_id),re.sub(r"_TE\.gtf", "", tissue), strand ]
                    stringreturn = "\t".join([str(x) for x in stringreturn])
                    f_out.write(stringreturn + "\n")


def extract_tss(reference_gtf,output):
    with open(output_file, "w") as f_out:
        for transciptID in reference_gtf:
            id = transciptID.split("\t")
            exons = reference_gtf[transciptID]
            exons.sort() 
            num_exons=int(len(exons)/2)
            chr = id[0]
            strand = id[1]
            if id[1]=="+":
                tss = exons[0]  
            elif id[1]=="-":
                tss = exons[-1]
            junctions = ",".join([str(x) for x in exons[1:-1]])
            stringreturn=[chr, int(tss)-1, int(tss) , id[2] ,".", strand ]
            stringreturn = "\t".join([str(x) for x in stringreturn])
            f_out.write(stringreturn + "\n")

def generate_random_string(length=8):
    letters = string.ascii_lowercase
    random_string = ''.join(random.choice(letters) for i in range(length))
    return random_string


def ref_merge(temp_file,TE_anno,reference_gtf,prefix):
    cmd = "bedtools intersect -wo -a "+ temp_file + "_ref_tss.bed -b "+ TE_anno + "  > " + temp_file +  "_TE_anno.bed "
    subprocess.check_call(cmd, shell=True, executable='/bin/bash')

    cmd = "cut -f 4   " + temp_file +  "_TE_anno.bed |sort|uniq > " + temp_file +  "_transcript.txt"
    subprocess.check_call(cmd, shell=True, executable='/bin/bash')

    cmd = "awk -F \"\\\"\" 'BEGIN{while(getline<\"" + temp_file +  "_transcript.txt\") a[$1]=1;} {if(a[$4]!=1) print $0}' "+ reference_gtf + " > " + temp_file +  "_reference_filtered.gtf "
    subprocess.check_call(cmd, shell=True, executable='/bin/bash')

    cmd = "cat " + prefix + ".gtf " + temp_file +  "_reference_filtered.gtf  > " + temp_file +  "_TE_merge.gtf"
    subprocess.check_call(cmd, shell=True, executable='/bin/bash')

    cmd = "gffread -T --sort-alpha -o " + args.prefix + "_sorted.gtf " + temp_file +  "_TE_merge.gtf"
    subprocess.check_call(cmd, shell=True, executable='/bin/bash')



transcripts = {}

gtf_list = []
with open(args.gtf_list) as files:
    for each in files:
        gtf_list.append(each.strip())
files.close()  

for gtf in gtf_list:
    transcripts[gtf] = process_gtf(gtf)

if args.reference_gtf:
    ref_transcript = process_gtf(args.reference_gtf)

merged_transcripts,merged_tss  = tss_merge(transcripts,args.tss_merge_distance,args.min_exon_length)

output_file = args.prefix + ".gtf"
write_gtf(merged_transcripts,merged_tss,output_file)

temp_file = generate_random_string()
output_file = temp_file + "_ref_tss.bed"
extract_tss(ref_transcript,output_file)
ref_merge(temp_file,args.TE_anno,args.reference_gtf,args.prefix)

output_file = args.prefix + "_tss.bed"
write_bed(merged_transcripts,merged_tss,output_file)

cmd = "rm " + temp_file +  "*"
subprocess.check_call(cmd, shell=True, executable='/bin/bash')
