import sys, csv, glob, re, os, argparse
import numpy as np
from . import config
# SUPPORT -u mode. will parse both form

def parsePredFile(path_to_peptide_file,fin,fout_pass,IC50_cutoff,gene_name,junction,pred_med,length):
	seq_list=[]
	seq_list_seq=[]
	if os.path.exists(path_to_peptide_file):
		path_seq_file=path_to_peptide_file
		for seq in open(path_seq_file):
			if seq.startswith('>'):
				seq_list.append(seq.strip('>').strip())
			else:
				seq_list_seq.append(seq.strip())
		for l in csv.DictReader(open(fin),dialect='excel-tab'):
			ic50=[k for k in l.keys() if k.find('ic50')!=-1]
			p = re.compile('\d+(\.\d+)?')
			ic50_value=[float(l[k]) for k in ic50 if p.match(l[k]) != None] 
			med_ic50=np.median(ic50_value)
			junc_pep=seq_list_seq[int(l['seq_num'])-1]
			if junc_pep.find(l['peptide'])==-1:
				pos= junc_pep.upper().find(l['peptide'])
				if junc_pep.upper()==junc_pep[:pos]+l['peptide']+junc_pep[pos+len(l['peptide']):]:#TODO
					if med_ic50<=IC50_cutoff:
						fout_pass.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(gene_name, junction,seq_list[int(l['seq_num'])-1],l['start'],l['end'],length,pred_med, l['allele'],l['peptide'],junc_pep,med_ic50))

def writePositivePrediction(output_path, IC50_cutoff):
	print '[INFO] Collecting binding predictions:'+output_path.rstrip('/').split('/')[-1]
	pred_fin_list=glob.glob(output_path+'/tmp/pred/*/*')
	fout_pass=open(output_path+'/pred_filtered.score'+str(IC50_cutoff)+'.txt','w')
	tot=len(pred_fin_list)-1
	for i,fin in enumerate(pred_fin_list):
		config.update_progress(i/(0.0+tot))
		pep_path=output_path+'/tmp/prot.compared/'
		file_name=fin.split('/')[-1]
		des=file_name.split('.')
		gene_name=des[1]
		length=des[-3]
		junction=des[-6]
		pred_med=des[-5]
		pep_path_skp=pep_path+'skp/'+'.'.join(file_name.split('.')[:-5])+'.prot.fa'
		pep_path_inc=pep_path+'inc/'+'.'.join(file_name.split('.')[:-5])+'.prot.fa'
		parsePredFile(pep_path_skp,fin,fout_pass,IC50_cutoff,gene_name,junction,pred_med,length)
		parsePredFile(pep_path_inc,fin,fout_pass,IC50_cutoff,gene_name,junction,pred_med,length)
	fout_pass.close()

def loadSampleHLA(sampleHLA_fin):
	sample_HLA={}
	sample_list=[]
	for l in open(sampleHLA_fin):
		ls=l.strip().split()
		for hla in ls[1:]:
			if hla not in sample_HLA:
				sample_HLA[hla]=[]
			sample_HLA[hla].append(ls[0])
		sample_list.append(ls[0])
	return sample_list,sample_HLA

def loadScreening(screening_result):
	screening_result_dict={}
	header_list=[]
	for n,l in enumerate(open(screening_result)):
		ls=l.strip().split('\t')
		if n==0:
			header_list=ls[1:]
			continue
		reformat_name=ls[0].split(':')
		reformat_name=':'.join([reformat_name[0],reformat_name[1].split('.')[0]]+reformat_name[2:])
		screening_result_dict[reformat_name]=ls[1:]
	return screening_result_dict,header_list

def loadGeneExp(gene_exp_matrix_fin):
	Exp={}
	i=0
	for l in open(gene_exp_matrix_fin):
		if i==0:
			i+=1
			continue
		ls=l.strip().split('\t')
		name=ls[0].split('_')[0].split('.')[0]
		exp_list=map(float, ls[1:])
		Q1=round(np.nanpercentile(exp_list,25),2)
		Q3=round(np.nanpercentile(exp_list,75),2)
		mean=round(np.nanmean(exp_list),2)
		Exp[name]=map(str,[mean, Q1, Q3])
	return Exp,['meanGeneExp','Q1GeneExp','Q3GeneExp']


def writePeptideSummary(output_path, screening_result_path, gene_exp_path, IC50_cutoff, sample_HLA, sample_list):
	pep_dict={}
	screening_result_dict, screening_result_header=loadScreening(screening_result_path)
	if gene_exp_path:
		gene_exp_dict, gene_exp_header=loadGeneExp(gene_exp_path)#exp. or fpkm1. matrix
	fout_screening=open(output_path+'/epitope_summary.peptide-based.txt','w')
	if gene_exp_path:
		fout_screening.write('\t'.join(['epitope','as_event','junction_peptide_form','num_hla','num_sample']+sample_list+screening_result_header+gene_exp_header)+'\n')
	else:
		fout_screening.write('\t'.join(['epitope','as_event','junction_peptide_form','num_hla','num_sample']+sample_list+screening_result_header)+'\n')

	for l in open(output_path+'/pred_filtered.score'+str(IC50_cutoff)+'.txt'):
		ls=l.strip().split('\t')
		peptide=ls[8]
		hla_type=ls[7]
		pred_score=ls[10]
		splicing=(ls[0]+'_'+ls[1]).replace('_',':')
		form=ls[2].split(':')[5][:3]
		if peptide+'|'+splicing+form not in pep_dict:
			pep_dict[peptide+'|'+splicing+form]=[]
		pep_dict[peptide+'|'+splicing+form].append(hla_type+'|'+pred_score)
#	print len(pep_dict)
	for k in pep_dict:
		ks=k.split('|')
		peptide=ks[0]
		splicing=ks[1][:-3]
		form=ks[1][-3:]
		hla_info=pep_dict[k]
		hla={}
		for info in hla_info:
			info_split=info.split('|')
			hla_type=info_split[0]
			pred_score=info_split[1]
			hla[hla_type]=pred_score
		num_hla=len(hla)
		line= [peptide, splicing, form, str(num_hla)]
		patient_count=len(sample_list)
		for p in sample_list:
			patient_hla_info=[]
			for h in hla:
				if p in sample_HLA[h]:
					patient_hla_info.append(h+'|'+hla[h])
			if patient_hla_info==[]:
				patient_hla_info=['-']
				patient_count-=1
			patient_hla_info_line=';'.join(patient_hla_info)
			line.append(patient_hla_info_line)
		line.insert(4,str(patient_count))
		if gene_exp_path:
			gene_exp_list=gene_exp_dict[splicing.split(':')[0]]
			fout_screening.write('\t'.join(line+screening_result_dict[splicing]+gene_exp_list)+'\n')
		else:
			fout_screening.write('\t'.join(line+screening_result_dict[splicing])+'\n')

def writeJunctionSummary(output_path, screening_result_path, gene_exp_path, IC50_cutoff, sample_HLA, sample_list):
	junction_dict={}
	screening_result_dict, screening_result_header=loadScreening(screening_result_path)
	if gene_exp_path:
		gene_exp_dict, gene_exp_header=loadGeneExp(gene_exp_path)#exp. or fpkm1. matrix
	fout_screening=open(output_path+'/epitope_summary.junction-based.txt','w')
	if gene_exp_path:
		fout_screening.write('\t'.join(['as_event','junction_peptide_form','num_hla','num_sample','hla_types']+sample_list+screening_result_header+gene_exp_header)+'\n')
	else:
		fout_screening.write('\t'.join(['as_event','junction_peptide_form','num_hla','num_sample','hla_types']+sample_list+screening_result_header)+'\n')

	for l in open(output_path+'/pred_filtered.score'+str(IC50_cutoff)+'.txt'):
		ls=l.strip().split('\t')
		peptide=ls[8]
		hla_type=ls[7]
		pred_score=ls[10]
		splicing=(ls[0]+'_'+ls[1]).replace('_',':')
		form=ls[2].split(':')[5][:3]
		# if peptide+'|'+splicing+form not in pep_dict:
		# 	pep_dict[peptide+'|'+splicing+form]=[]
		# pep_dict[peptide+'|'+splicing+form].append(hla_type+'|'+pred_score)
		if splicing+form not in junction_dict:
			junction_dict[splicing+form]=[]
		junction_dict[splicing+form].append(hla_type+'|'+pred_score)	
#	print len(pep_dict)
	for k in junction_dict:
		splicing=k[:-3]
		form=k[-3:]
		hla_info=junction_dict[k]
		hla={}
		for info in hla_info:
			info_split=info.split('|')
			hla_type=info_split[0]
			pred_score=info_split[1]
			if hla_type not in hla:
				hla[hla_type]=set()
			hla[hla_type].add(pred_score)
		num_hla=len(hla)
		line= [splicing, form, str(num_hla),';'.join(sorted(hla.keys()))]
		patient_count=len(sample_list)
		for p in sample_list:
			patient_hla_info=[]
			for h in hla:
				if p in sample_HLA[h]:
					patient_hla_info.append(h+'|'+','.join(list(hla[h])))
			if patient_hla_info==[]:
				patient_hla_info=['-']
				patient_count-=1
			patient_hla_info_line=';'.join(patient_hla_info)
			line.append(patient_hla_info_line)
		line.insert(3,str(patient_count))
		if gene_exp_path:
			gene_exp_list=gene_exp_dict[splicing.split(':')[0]]
			fout_screening.write('\t'.join(line+screening_result_dict[splicing]+gene_exp_list)+'\n')
		else:
			fout_screening.write('\t'.join(line+screening_result_dict[splicing])+'\n')

def main(args):
	outdir=args.outdir
	IC50_cutoff=args.ic50_cut_off
	gene_exp_matrix=args.gene_exp_matrix

	writePositivePrediction(outdir+'/primary', IC50_cutoff)
	writePositivePrediction(outdir+'/prioritized', IC50_cutoff)
	
	sample_list,sample_HLA=loadSampleHLA(args.mhc_by_sample)
	analysis_name=[l.strip() for l in open(args.parameter_fin)][0]
	screening_result_path=outdir+'/'+analysis_name+'.primary.txt'
	writePeptideSummary(outdir+'/primary', outdir+'/'+analysis_name+'.primary.txt', gene_exp_matrix, IC50_cutoff, sample_HLA, sample_list)
	writePeptideSummary(outdir+'/prioritized', outdir+'/'+analysis_name+'.prioritized.txt', gene_exp_matrix, IC50_cutoff, sample_HLA, sample_list)

	writeJunctionSummary(outdir+'/primary', outdir+'/'+analysis_name+'.primary.txt', gene_exp_matrix, IC50_cutoff, sample_HLA, sample_list)
	writeJunctionSummary(outdir+'/prioritized', outdir+'/'+analysis_name+'.prioritized.txt', gene_exp_matrix, IC50_cutoff, sample_HLA, sample_list)
	

if __name__ == '__main__':
	main()
