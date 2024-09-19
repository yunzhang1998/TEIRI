#!/usr/bin/python
import argparse
import sys
import re
import random
import string
import time
import subprocess
import logging

parser = argparse.ArgumentParser(description='TEIRI_merge')
parser.add_argument('-i', '--gtf_list', help="A text file with a list of GTF files (required)")
parser.add_argument('-r', '--reference_gtf', help="Reference genome annotation in GTF format (required)")
parser.add_argument('-l', '--corrected_bed12', help="The flair-corrected bed12 file")
parser.add_argument('--corrected_tss', help="A tsv file with the corrected TE-derived TSSs (TEIRI_correct.py generated)")
parser.add_argument('--corrected_tss_single', help="A tsv file with the corrected TE-derived TSSs for the single-exon transcript (TEIRI_correct.py generated)")
parser.add_argument('--TGS_weight', default=1000, type=float, help='The weight of TGS reads supporting the transcript (default: 1000)')
parser.add_argument( '--illumina_threshold', default=0.05, type=float, help='The min NGS ratio supporting the transcript (default: 0.05)')
parser.add_argument( '--nanopore_threshold', default=2, type=float, help='The min TGS reads supporting the transcript (default: 2)')
parser.add_argument( '--ref_transcript_length', default=2725, type=float, help='The average length of reference transcripts (default: 2725)')
parser.add_argument('--max_transcripts', default=10, type=int, help='The max counts of transcripts for a TE-initiated RNA (default: 10)')
parser.add_argument('--trunctated_exclude', default=True, help='Truncated transcripts were excluded, as they may represent fragments of the full-length transcript (default: True)')
parser.add_argument('--min_transcript_length', default=200, type=float, help='The min transcript length (default: 200)')
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

def process_bed12(bed12_path):
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

    return transcript

def process_bed(bed_path):
    te_tss = {}

    with open(bed_path) as bed:
        for eachline in bed:
            if not eachline.strip().startswith('#') :
                temp = eachline.strip().split("\t")
                chr="\t".join([temp[0],temp[5]])
                if chr not in te_tss:
                    te_tss[chr] = {}
                id=[str(int(temp[1])+1),str(temp[2])]
                if temp[5]=="-":
                    id.reverse()
                if id[1] not in te_tss[chr]:
                    te_tss[chr][id[1]]=[str(id[0])]
                else:
                    te_tss[chr][id[1]].append(str(id[0]))
    
    return te_tss



def splicing_sites_extract(reference_gtf):
    ss1 = {}
    ss2 = {}
    
    for transcriptID in reference_gtf:
        id = transcriptID.split("\t")
        exons=reference_gtf[transcriptID]
        exons.sort() 
        num_exons=int(len(exons)/2)
        chr = "\t".join([id[0],id[1]])
        if  num_exons>=2 :
            if chr not in ss1:
                ss1[chr]={}
                ss2[chr]={}
            for i in range(0,num_exons-1):
                ss1_temp = str(exons[2*i+2])
                ss2_temp = str(exons[2*i+1])
                if ss1_temp not in ss1[chr]:
                    ss1[chr][ss1_temp]=[transcriptID]
                else:
                    ss1[chr][ss1_temp].append(transcriptID)
                if ss2_temp not in ss2[chr]:
                    ss2[chr][ss2_temp]=[transcriptID]
                else:
                    ss2[chr][ss2_temp].append(transcriptID)
    return ss1,ss2

def ref_exons_generate(reference_gtf,ss,i,chr,strand,ss1,ss2):
    reference_transcripts=[]
    if strand == "+":
        if ss in ss2[chr] and i % 2==1:
            reference_transcripts = ss2[chr][ss]
        elif ss in ss1[chr] and i % 2==0:
            reference_transcripts = ss1[chr][ss]
        if len(reference_transcripts)>=1:
            for ID in reference_transcripts:
                ref_exons = reference_gtf[ID]
                ref_exons.sort()
                ref_exons = ref_exons[1:]
                ref_exons = list(map(str, ref_exons))
                ref_exons = ref_exons[ref_exons.index(ss):]
                yield ID,ref_exons
        else:
            yield "None","None"

    elif strand == "-":
        if ss in ss1[chr] and i % 2==1:
            reference_transcripts = ss1[chr][ss]
        elif ss in ss2[chr] and i % 2==0:
            reference_transcripts = ss2[chr][ss]
        if len(reference_transcripts)>=1:
            for ID in reference_transcripts:
                ref_exons = reference_gtf[ID]
                ref_exons.sort(reverse=True)
                ref_exons = ref_exons[1:]
                ref_exons = list(map(str, ref_exons))
                ref_exons = ref_exons[ref_exons.index(ss):]
                yield ID,ref_exons
        else:
            yield "None","None"





def transcript_correct(transcript,reference_gtf,te_tss,ss1,ss2):
    corrected_transcript={}
    for transcriptID in list(transcript.keys()):
        id = transcriptID.split("\t")
        exons=transcript[transcriptID]
        num_exons=int(len(exons)/2)
        chr = "\t".join([id[0],id[1]])
        if id[1]=="+":
            exons.sort() 
        elif id[1]=="-":
            exons.sort(reverse=True)
        if chr in te_tss:
            if  (num_exons>=2) and (str(exons[1]) in te_tss[chr]):
                new_exons = ["\t".join(te_tss[chr][str(exons[1])])]
                correct = 0
                if chr in ss1:
                    for i in range(1,2*num_exons-1):
                        for refID,ref_exons in ref_exons_generate(reference_gtf,str(exons[i]),i,chr,id[1],ss1,ss2):
                            if refID !="None":
                                new_transcriptID = "_".join([transcriptID , refID ])
                                corrected_transcript[new_transcriptID]= new_exons + ref_exons
                                correct = 1
                        if correct == 1:
                            break
                        else:
                            new_exons.append(str(exons[i]))
                if correct != 1 : 
                    new_exons.append(str(exons[-1]))
                    corrected_transcript[transcriptID] = new_exons
    return corrected_transcript

def transcript_merge(corrected_transcript,TGS_weight,ref_transcript_length,illumina_threshold,nanopore_threshold,trunctated_exclude,max_transcripts,min_transcript_length):
    mergeded_transcript = {}
    for gtf in corrected_transcript:
        for transcriptID in list(corrected_transcript[gtf].keys()):
            id = transcriptID.split("\t")
            exons=corrected_transcript[gtf][transcriptID]
            ss_all = "\t".join(exons[1:-1])
            num_exons=int(len(exons)/2)
            chr = "\t".join([id[0],id[1]])
            for tss_ in exons[0].split("\t"):
                if chr not in mergeded_transcript:
                    mergeded_transcript[chr] = {}
                if tss_ not in mergeded_transcript[chr]:
                    mergeded_transcript[chr][tss_]={}
                if ss_all not in mergeded_transcript[chr][tss_]:
                    if gtf == "ref" :
                        mergeded_transcript[chr][tss_][ss_all] = [[],1000,int(tss_),int(exons[-1])]
                    elif gtf == "nanopore" :
                        mergeded_transcript[chr][tss_][ss_all] = [[],1,int(tss_),int(exons[-1])]
                    else:
                        mergeded_transcript[chr][tss_][ss_all] = [[gtf],0,int(tss_),int(exons[-1])]
                else:
                    if gtf == "ref" :
                        mergeded_transcript[chr][tss_][ss_all][1] += 1000
                    elif gtf == "nanopore" :
                        mergeded_transcript[chr][tss_][ss_all][1] += 1
                    else:
                        mergeded_transcript[chr][tss_][ss_all][0].append(gtf)
                    if (id[1]=="+") and (int(exons[-1]) < mergeded_transcript[chr][tss_][ss_all][3]):
                        mergeded_transcript[chr][tss_][ss_all][3] = int(exons[-1])
                    if (id[1]=="-") and (int(exons[-1]) > mergeded_transcript[chr][tss_][ss_all][3]):
                        mergeded_transcript[chr][tss_][ss_all][3] = int(exons[-1])


    for chr in mergeded_transcript:
        for tss_ in list(mergeded_transcript[chr].keys()):
            for ss_all in list(mergeded_transcript[chr][tss_].keys()):
                exons = [mergeded_transcript[chr][tss_][ss_all][2]] + ss_all.split("\t") + [mergeded_transcript[chr][tss_][ss_all][3]]
                
                transcript_length = float( transcript_length_calculate(exons))
                nanopore_reads = mergeded_transcript[chr][tss_][ss_all][1]
                illumina_score = len(list(set(mergeded_transcript[chr][tss_][ss_all][0]))) 
                if args.corrected_bed12:
                    nanopore_score = nanopore_reads*TGS_weight*transcript_length/ref_transcript_length
                    weight_score = illumina_score + nanopore_score
                else:
                    weight_score = illumina_score
                if illumina_score < illumina_threshold and nanopore_reads < nanopore_threshold :
                    del mergeded_transcript[chr][tss_][ss_all]
                elif transcript_length < min_transcript_length :
                    del mergeded_transcript[chr][tss_][ss_all]
                else:
                    mergeded_transcript[chr][tss_][ss_all][0]=float(weight_score)

    if trunctated_exclude:
        for chr in mergeded_transcript:
            for tss_ in list(mergeded_transcript[chr].keys()):
                for ss_all in list(mergeded_transcript[chr][tss_].keys()):
                    for ss_all_ in list(mergeded_transcript[chr][tss_].keys()):
                        if (ss_all in ss_all_) and ss_all.split("\t")[0]==ss_all_.split("\t")[0]  and mergeded_transcript[chr][tss_][ss_all][0] < mergeded_transcript[chr][tss_][ss_all_][0]:
                            del mergeded_transcript[chr][tss_][ss_all]
                            break

    final_transcript = {}
    for chr in mergeded_transcript:
        final_transcript[chr] = []
        for tss_ in mergeded_transcript[chr]:
            score = []
            for ss_all in mergeded_transcript[chr][tss_]:
                score.append(float( mergeded_transcript[chr][tss_][ss_all][0]))
            score.sort(reverse=True)
            max_score = set(score[0:max_transcripts])
            for ss_all in mergeded_transcript[chr][tss_]:
                exons = [tss_] + ss_all.split("\t") + [mergeded_transcript[chr][tss_][ss_all][3]]
                weight_score = mergeded_transcript[chr][tss_][ss_all][0]
                if weight_score in max_score :
                    new_transcript = "\t".join([str(x) for x in exons])
                    final_transcript[chr].append(new_transcript)

    return final_transcript



def transcript_length_calculate(exons_temp):
    exons_temp = list(map(int, exons_temp))
    exons_temp.sort()
    num_exons=int(len(exons_temp)/2)
    transcript_length_temp = 0
    if  num_exons>=2 :
        for i in range(0,num_exons):
            transcript_length_temp += (int(exons_temp[2*i+1])-int(exons_temp[2*i]))
    return (transcript_length_temp)


def SE_merge(single_te_tss):
    transcript_SE={}
    for chr in single_te_tss:
        transcript_SE[chr] = {} 
        strand = chr.split("\t")[1]
        for tes_ in single_te_tss[chr]:
            for tss_ in single_te_tss[chr][tes_]:
                if tss_ not in transcript_SE[chr] :
                    transcript_SE[chr][tss_] = tes_
                elif strand=="+" and int(tes_) > int(transcript_SE[chr][tss_]) :
                    transcript_SE[chr][tss_] = tes_
                elif strand=="-" and int(tes_) < int(transcript_SE[chr][tss_]) :
                    transcript_SE[chr][tss_] = tes_
    
    final_transcript = {}
    for chr in transcript_SE:
        final_transcript[chr] = []
        for tss_ in transcript_SE[chr]:
            exons = [tss_] + [transcript_SE[chr][tss_]]
            new_transcript = "\t".join([str(x) for x in exons])
            final_transcript[chr].append(new_transcript)
    return final_transcript

def write_exon(exons,transcriptid,chr_,strand):
    exons = list(map(int, exons))
    num_exons=int(len(exons)/2)
    exons.sort()
    exons = list(map(str, exons))
    stringreturn=[chr_, "TE_RNA","transcript" , exons[0] ,exons[-1],".", strand,".",str("transcript_id \"teRNA_"+ str(transcriptid) + "\";") ]
    stringreturn = "\t".join([str(x) for x in stringreturn])
    yield stringreturn
    for i in range(0,num_exons):
        stringreturn=[chr_, "TE_RNA","exon" , exons[2*i] ,exons[2*i+1],".", strand,".",str("transcript_id \"teRNA_"+ str(transcriptid) + "\";") ]
        stringreturn = "\t".join([str(x) for x in stringreturn])
        yield stringreturn

def write_gtf(mergeded_transcripts,SE_transcripts,output_file):
    transcriptid = 0

    with open(output_file, "w") as f_out:
        for chr in mergeded_transcripts:
            strand= chr.split("\t")[1]
            chr_= chr.split("\t")[0]
            for id in mergeded_transcripts[chr]:
                transcriptid += 1
                exons = id.split("\t")
                for stringreturn in write_exon(exons,transcriptid,chr_,strand):
                    f_out.write(stringreturn + "\n")
            if chr in SE_transcripts:
                for id in SE_transcripts[chr]:
                    transcriptid += 1
                    exons = id.split("\t")
                    for stringreturn in write_exon(exons,transcriptid,chr_,strand):
                        f_out.write(stringreturn + "\n")


                        

transcripts = {}

gtf_list = []
with open(args.gtf_list) as files:
    for each in files:
        gtf_list.append(each.strip())
files.close()  

for gtf in gtf_list:
    transcripts[gtf] = process_gtf(gtf)

if args.reference_gtf:
    transcripts["ref"] = process_gtf(args.reference_gtf)

if args.corrected_bed12:
    transcripts["nanopore"] = process_bed12(args.corrected_bed12)

te_tss = process_bed(args.corrected_tss)

ss1,ss2 = splicing_sites_extract(transcripts["ref"])

corrected_transcripts = {}

for gtf in transcripts:
    corrected_transcripts[gtf] = transcript_correct(transcripts[gtf],transcripts["ref"],te_tss,ss1,ss2)
    print (gtf+"_corrected")

min_samples = len(gtf_list) * args.illumina_threshold

mergeded_transcripts=transcript_merge(corrected_transcripts,args.TGS_weight,args.ref_transcript_length,min_samples,args.nanopore_threshold,args.trunctated_exclude,args.max_transcripts,args.min_transcript_length)

single_te_tss = process_bed(args.corrected_tss_single)

SE_transcripts = SE_merge(single_te_tss)

output_file = args.prefix + "_TE.gtf"
write_gtf(mergeded_transcripts,SE_transcripts,output_file)
