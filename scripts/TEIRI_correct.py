#!/usr/bin/python
import argparse
import sys
import re
import random
import string
import time
import subprocess
import logging

parser = argparse.ArgumentParser(description='TEIRI_correct')
parser.add_argument('-i', '--gtf_list', help="A text file with a list of GTF files (required)")
parser.add_argument('-r', '--reference_gtf', help="Reference genome annotation in GTF format (required)")
parser.add_argument('-l', '--corrected_bed12', help="The flair-corrected bed12 file")
parser.add_argument('-c', '--tss_score', help="A tsv file with the counts of unique mapping reads supporting the tss (CAGE and/or RAMPAGE)")
parser.add_argument('--TE_anno', help="TE annotation in BED format (required)")
parser.add_argument('--max_exon1_length', default=2588, type=int, help='The max length of first exon (default: 2588)')
parser.add_argument('--min_exon1_length', default=20, type=int, help='The min length of first exon (default: 20)')
parser.add_argument( '--tss_window', default=50, type=int, help='The window size used to calculate the weight score (default: 50)')
parser.add_argument('--max_tss', default=2, type=int, help='Maximum number of TSS picked per first exon (default: 2)')
parser.add_argument('--TGS_weight', default=1, type=float, help='The weight of TGS reads supporting the TSS (default: 1)')
parser.add_argument('--min_NGS_ratio', default=0.05, type=float, help='The min NGS ratio supporting the first exon (default: 0.05)')
parser.add_argument('--min_TGS_reads', default=2, type=float, help='The min TGS reads supporting the first exon (default: 2)')
parser.add_argument('--threshold', default=10, type=float, help='The threshold of weight score for TSS (default: 10)')
parser.add_argument('--min_NGS_ratio_SE', default=0.25, type=float, help='The min NGS ratio supporting the single exon (default: 0.25)')
parser.add_argument('--min_TGS_reads_SE', default=10, type=float, help='The min TGS reads supporting the single exon (default: 10)')
parser.add_argument('--SE_threshold', default=20, type=float, help='The threshold of weight score for TSS of the single exon (default: 20)')

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

def process_bed(bed_path):
    tss = {}
    tss_set = {}

    with open(bed_path) as bed:
            for eachline in bed:
                if not eachline.strip().startswith('#') :
                    temp = eachline.strip().split("\t")
                    chr="\t".join([temp[0],temp[3]])
                    tss_temp=str(int(temp[1]))
                    if chr not in tss:
                        tss[chr] = {}
                        tss_set[chr] = set()
                    tss[chr][tss_temp] = int(temp[4])
                    tss_set[chr].add(int(tss_temp))

    return tss,tss_set

def process_bed12(bed12_path,tss_count,tss_set,TGS_weight):
    transcript = {}

    with open(bed12_path) as nanopore:
        for eachline in nanopore:
            if not eachline.strip().startswith('#'):
                temp = eachline.strip().split("\t")
                id = "\t".join([temp[0],str(temp[5]),temp[3]])
                exon_length = temp[10].split(",")
                exon_start = temp[11].split(",")
                transcript[id] = []
                for i in range(0,int(temp[9])):
                    transcript[id].append(int(temp[1]) + int(exon_start[i]) + 1)
                    transcript[id].append(int(temp[1]) + int(exon_start[i]) + int(exon_length[i]))
                if int(temp[9]) >= 2 :
                    chr = "\t".join([temp[0],str(temp[5])])
                    if temp[5]=="+":
                        tss_temp=str(int(temp[1]))
                    elif temp[5]=="-":
                        tss_temp=str(int(temp[1])+ int(exon_start[-2]) + int(exon_length[-2]) - 1)
                    if chr not in tss_count:
                        tss_count[chr] ={}
                        tss_set[chr] = set()
                    if tss_temp not in tss_count[chr]:
                        tss_count[chr][tss_temp] = 0
                        tss_set[chr].add(int(tss_temp))
                    tss_count[chr][tss_temp] += TGS_weight
                    

    return transcript,tss_count,tss_set


def exon_extract(transcipt,min_exon1_length):
    output_exons = {}
    output_exons["first_exon"] = {}
    output_exons["inner_exon"] = {}
    output_exons["final_exon"] = {}
    output_exons["single_exon"] = {}

    for gtf in transcipt:
        for transciptID in transcipt[gtf]:
            id = transciptID.split("\t")
            exons=transcipt[gtf][transciptID]
            exons.sort() 
            num_exons=int(len(exons)/2)
            chr = "\t".join([id[0],id[1]])
            if num_exons>=2 :
                exons[0] = exons[1] - min_exon1_length
                exons[-1] = exons[-2] + min_exon1_length
                for i in range(num_exons):
                    exon_temp="\t".join([str(exons[2*i]-1),str(exons[2*i+1])])
                    exon_type = "inner_exon"
                    if id[1]=="+":
                        if i==0:
                            exon_type = "first_exon"
                        elif i==(num_exons-1):
                            exon_type = "final_exon"
                    elif id[1]=="-":
                        if i==0:
                            exon_type = "final_exon"
                        elif i==(num_exons-1):
                            exon_type = "first_exon"
                    if chr not in output_exons[exon_type]:
                        output_exons[exon_type][chr]={}
                    if exon_temp not in output_exons[exon_type][chr]:
                        output_exons[exon_type][chr][exon_temp] = [[],0,"no_ref"]
                    if gtf == "reference":
                        output_exons[exon_type][chr][exon_temp][2]="ref"
                    elif gtf == "nanopore":
                        output_exons[exon_type][chr][exon_temp][1]+=1
                    else:
                        output_exons[exon_type][chr][exon_temp][0].append(gtf)
            elif num_exons==1 :
                exon_temp="\t".join([str(exons[0]-1),str(exons[1])])
                exon_type = "single_exon"
                chr_ = id[0]
                if chr_ not in output_exons[exon_type]:
                    output_exons[exon_type][chr_]={}
                if exon_temp not in output_exons[exon_type][chr_]:
                    output_exons[exon_type][chr_][exon_temp] = [0,0,0]
                if gtf == "reference":
                    output_exons[exon_type][chr_][exon_temp][2]=1
                elif gtf == "nanopore":
                    output_exons[exon_type][chr_][exon_temp][1]+=1
                else:
                    output_exons[exon_type][chr_][exon_temp][0]+=1
    
    return(output_exons)

def generate_random_string(length=8):
    logging.info("Generating random string for temp file name...")
    letters = string.ascii_lowercase
    random_string = ''.join(random.choice(letters) for i in range(length))
    logging.info("Random string generation completed.")
    return random_string

def write_exon_bed_file(results, output_file,min_samples,min_TGS_reads,single_exon):
    logging.info(f'Writing BED file to: {output_file}')
    with open(output_file, "w") as f_out:
        for chr in results:
            for exon_temp in results[chr]:
                nanopore_reads = results[chr][exon_temp][1]
                ref = results[chr][exon_temp][2]
                if single_exon==1 :
                    sample_counts = results[chr][exon_temp][0]
                    result = "\t".join([chr,exon_temp,str(sample_counts),str(nanopore_reads),"+",str(ref)])
                    f_out.write(result + "\n")
                    if ref==1 and (sample_counts!=0 or nanopore_reads!=0) :
                        result = "\t".join([chr,exon_temp,str(sample_counts),str(nanopore_reads),"+",str(0)])
                        f_out.write(result + "\n")
                else:
                    samples = list(set(results[chr][exon_temp][0]))
                    chr_ = chr.split("\t")
                    result = "\t".join([chr_[0],exon_temp,ref,ref,chr_[1]]) 
                    if len(samples) >= min_samples or nanopore_reads >= min_TGS_reads   :
                        f_out.write(result + "\n")
                    elif ref=="ref" and (len(samples)!=0 or nanopore_reads!=0) :
                        f_out.write(result + "\n")
    logging.info('BED file writing completed.')


def exon_filter(temp_file):
    first = temp_file  + "_first_exon.bed"
    inner = temp_file  + "_inner_exon.bed" 
    final = temp_file  + "_final_exon.bed"
    single = temp_file  + "_single_exon.bed"

    cmd = "cat " + inner + " " + final + " > " + temp_file +  "_inner_final.bed"
    subprocess.check_call(cmd, shell=True, executable='/bin/bash')
    cmd = "awk '$4==\"ref\"' "  + first  + " > " + temp_file +  "_reference_exon1.bed"
    subprocess.check_call(cmd, shell=True, executable='/bin/bash')
    cmd = "bedtools subtract -a " + first + " -b " + temp_file +  "_inner_final.bed  -s -A > " + temp_file +  "_filtered_exon1.bed"
    subprocess.check_call(cmd, shell=True, executable='/bin/bash')
    cmd = "cat " + temp_file +  "_filtered_exon1.bed " + temp_file +  "_reference_exon1.bed|sort -k1,1 -k2,2n|uniq > " + temp_file +  "_final_exon1.bed"
    subprocess.check_call(cmd, shell=True, executable='/bin/bash')
    cmd = "cat " + first   + " " + inner  + " " + final    + " > " + temp_file +  "_total.bed"
    subprocess.check_call(cmd, shell=True, executable='/bin/bash')
    cmd = "bedtools subtract -a " +  single + " -b " + temp_file +  "_total.bed   -A > " + temp_file +  "_filtered_single_exon.bed"
    subprocess.check_call(cmd, shell=True, executable='/bin/bash')

def exon1_deal(filter_exon1,tss_count,tss_set,min_exon1_length,max_exon1_length,max_tss,tss_window,threshold):
    with open(filter_exon1) as file:
        for eachline in file:
            temp = eachline.strip().split("\t")
            chr = "\t".join([temp[0],temp[5]])
            if temp[5] == "+" :
                tss_range=set(range(int(temp[2])-max_exon1_length,int(temp[2])-min_exon1_length+1))
                tes = temp[2]
            elif temp[5] == "-" :
                tss_range=set(range(int(temp[1])+min_exon1_length-1,int(temp[1])+max_exon1_length))
                tes = temp[1]
            outputs=list(tss_correct(chr,temp[5],tss_range,tes,tss_count,tss_set,max_tss,tss_window,threshold))
            if len(outputs) >=1:
                for output in outputs:
                    yield output


def tss_correct(chr,strand,tss_range,tes,tss_count,tss_set,max_tss,tss_window,threshold):
    if chr in tss_set:
        tss_range = list(tss_range.intersection(tss_set[chr]))
        best_tss = {}
        for k in range(0,int(max_tss)):
            if len(tss_range) >= 1 :
                best_tss[str(k)] = [0,0,0]
                weighted_score = dict.fromkeys(tss_range, 0) 
                if strand == "+" :
                    sorted_tss=sorted(tss_range,reverse=True)
                elif strand == "-":
                    sorted_tss=sorted(tss_range)
                for tss_1 in sorted_tss:
                    for tss_2 in sorted_tss:
                        if tss_1 == tss_2:
                            weighted_score[tss_1] += tss_count[chr][str(tss_1)]
                        elif (abs(tss_1 - tss_2) < tss_window) :
                            weighted_score[tss_1] += ((tss_window - abs(tss_1 - tss_2))/float(tss_window * tss_count[chr][str(tss_2)]))
                    if weighted_score[tss_1] > best_tss[str(k)][2] :
                        best_tss[str(k)] = [tss_1,tss_count[chr][str(tss_1)],weighted_score[tss_1] ]
                
                if strand == "+" :
                    output = [chr.split("\t")[0], best_tss[str(k)][0],tes,best_tss[str(k)][1],best_tss[str(k)][2],strand]
                elif strand == "-":
                    output = [chr.split("\t")[0], tes,best_tss[str(k)][0] + 1, best_tss[str(k)][1],best_tss[str(k)][2],strand]
                temp = sorted_tss
                for tss_1 in temp:
                    if abs(tss_1 - best_tss[str(k)][0]) < tss_window:
                        tss_range.remove(tss_1)
                if k==0 and best_tss[str(k)][2] >= threshold:
                    yield output
                elif best_tss[str(k)][2] >= threshold and best_tss[str(k)][2]/best_tss[str(0)][2] > 0.5:
                    yield output


def write_bed_file(corrected_tss, output_file):
    logging.info(f'Writing BED file to: {output_file}')
    with open(output_file, "w") as f_out:
        for tss in corrected_tss:
            f_out.write("\t".join(map(str, tss)) + "\n")
    logging.info('BED file writing completed.')


def te_tss(tss_bed,TE_anno,temp_file):
    cmd = "bedtools intersect -wo -a " +  tss_bed + " -b " + TE_anno +  " > " + temp_file +  "_intersect.bed"
    subprocess.check_call(cmd, shell=True, executable='/bin/bash')
    intersect_file = temp_file +  "_intersect.bed"
    with open(intersect_file) as file:
        for eachline in file:
            temp = eachline.strip().split("\t")
            if temp[5] == "+" and int(temp[1]) < int(temp[8]) and int(temp[1]) > int(temp[7])  :
                yield temp[:6]
            elif temp[5] == "-" and int(temp[2]) < int(temp[8]) and int(temp[2]) > int(temp[7])  :
                yield temp[:6]


def single_exon_deal(single,output_file,TE_anno,temp_file,min_samples,min_TGS_reads,tss_count,tss_set,max_tss,tss_window,threshold):
    filter_single = temp_file +  "_filter2_single_exon.bed"
    merge_file = temp_file +  "_merge.bed"
    final_single = temp_file +  "_final_single_exon.bed"
    intersect_file = temp_file +  "_intersect.bed"

    cmd = "bedtools subtract -a " +  single + " -b " + output_file +" -A | sort -k1,1 -k2,2n > " + filter_single
    subprocess.check_call(cmd, shell=True, executable='/bin/bash')
    cmd = "bedtools merge -i " +  filter_single + " -c 4 -o sum | awk -v OFS='\t' '$4>=" + str(min_samples) +  " ' |cut -f 1-3 >  " + merge_file
    subprocess.check_call(cmd, shell=True, executable='/bin/bash')
    cmd = "bedtools merge -i " +  filter_single + " -c 5 -o sum | awk -v OFS='\t' '$4>=" + str(min_TGS_reads) +  " ' |cut -f 1-3 >>  " + merge_file
    subprocess.check_call(cmd, shell=True, executable='/bin/bash')
    cmd = "bedtools merge -i " +  filter_single + " -c 7 -o sum | awk -v OFS='\t' '$4>=1' |cut -f 1-3 >> "  + merge_file
    subprocess.check_call(cmd, shell=True, executable='/bin/bash')
    cmd = "bedtools intersect -wo -a " +  filter_single + " -b " + merge_file +  " | awk -v OFS='\t' '$7==0' |cut -f 1-6 |sort -k1,1 -k2,2n|uniq > " + final_single
    subprocess.check_call(cmd, shell=True, executable='/bin/bash')
    cmd = "bedtools intersect -wo -a " +  final_single + " -b " + TE_anno +  " > " + intersect_file
    subprocess.check_call(cmd, shell=True, executable='/bin/bash')
    
    single_exon = {}
    t = 0
    old_exon_id = "NA"
    with open(intersect_file) as file:
        for eachline in file:
            temp = eachline.strip().split("\t")
            te_id = ",".join(temp[6:12])
            exon_id = ",".join(temp[:3])
            if t == 0:
                old_exon_id = exon_id
                te_list = [te_id]
                strand = [temp[11]]
                t = t + 1
            elif exon_id == old_exon_id:
                te_list.append(te_id)
                strand.append(temp[11])
            else:
                if len(list(set(strand))) ==1 :
                    chr = "\t".join([old_exon_id.split(",")[0],strand[0]])
                    te_temp = "\t".join(te_list)
                    if chr not in single_exon:
                        single_exon[chr]={}
                    if te_temp not in single_exon[chr]:
                        single_exon[chr][te_temp]=[]
                    single_exon[chr][te_temp].append(old_exon_id)   
                old_exon_id = exon_id
                te_list = [te_id]
                strand = [temp[11]]

    for chr in single_exon:
        for te_temp in single_exon[chr]:
            temp = chr.split("\t")
            start = []
            end = []
            for exon_id in single_exon[chr][te_temp]:
                start.append(int(exon_id.split(",")[1]))
                end.append(int(exon_id.split(",")[2]))
            if temp[1]=="+":
                tss_range = set(range(min(start),max(start)+1))
                tes = max(end)
            elif temp[1]=="-":
                tss_range = set(range(min(end)-1,max(end)))
                tes = min(start)
            outputs=list(tss_correct(chr,temp[1],tss_range,tes,tss_count,tss_set,max_tss,tss_window,threshold))
            if len(outputs) >=1:
                for output in outputs:
                    yield output




    

transcipts = {}
extracted_exons = {}

gtf_list = []
with open(args.gtf_list) as files:
    for each in files:
        gtf_list.append(each.strip())
files.close()  

for gtf in gtf_list:
    transcipts[gtf] = process_gtf(gtf)

if args.reference_gtf:
    transcipts["reference"] = process_gtf(args.reference_gtf)

if args.tss_score:
    tss_count,tss_set = process_bed(args.tss_score)

if args.corrected_bed12:
    transcipts["nanopore"],tss_count,tss_set = process_bed12(args.corrected_bed12,tss_count,tss_set,args.TGS_weight)

output_exons = exon_extract(transcipts,args.min_exon1_length)
temp_file = generate_random_string()

min_samples = args.min_NGS_ratio * len(gtf_list)
min_samples_SE = args.min_NGS_ratio_SE * len(gtf_list)

for exon_type in output_exons:
    output_file = temp_file + "_" + exon_type + ".bed"
    if exon_type == "single_exon":
        write_exon_bed_file(output_exons[exon_type],output_file,min_samples_SE,args.min_TGS_reads_SE,1)
    else:
        write_exon_bed_file(output_exons[exon_type],output_file,min_samples,args.min_TGS_reads,0)
    

exon_filter(temp_file)

filter_exon1 = temp_file  + "_final_exon1.bed"
corrected_tss = exon1_deal(filter_exon1,tss_count,tss_set,args.min_exon1_length,args.max_exon1_length,args.max_tss,args.tss_window,args.threshold)

corrected_tss_file = temp_file + "_tss_corrected.bed"
write_bed_file(corrected_tss,corrected_tss_file)

output_file = args.prefix + "_te_tss.bed"
write_bed_file(te_tss(corrected_tss_file,args.TE_anno,temp_file),output_file)



filtered_single_exon = temp_file +  "_filtered_single_exon.bed"
corrected_tss = single_exon_deal(filtered_single_exon,output_file,args.TE_anno,temp_file,min_samples_SE,args.min_TGS_reads_SE,tss_count,tss_set,1,args.tss_window,args.SE_threshold)

corrected_tss_file = temp_file + "_tss_corrected_SE.bed"
write_bed_file(corrected_tss,corrected_tss_file)

output_file = args.prefix + "_te_tss_SE.bed"
write_bed_file(te_tss(corrected_tss_file,args.TE_anno,temp_file),output_file)

cmd = "rm " +  temp_file + "*"
subprocess.check_call(cmd, shell=True, executable='/bin/bash')
