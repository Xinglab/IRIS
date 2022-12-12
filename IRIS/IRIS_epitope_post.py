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

def retrievePrioritizedfromPrimary(output_path, splicing_event_type): 
	fin_list=glob.glob(output_path+'/tmp/prot.compared/skp/*.fa')
	fin_list=fin_list+glob.glob(output_path+'/tmp/prot.compared/inc/*.fa')
	pred_fin_list=[]
	retrieve_path='/'.join(output_path.split('/')[:-1])+'/'+splicing_event_type+'.tier1'
	for fin in fin_list:
		event_name_list=fin.split('/')[-1].split('.')
		prediction_result_fn=retrieve_path+'/tmp/pred/'+event_name_list[1].split('_')[0]+'/'+'.'.join(event_name_list[:-2])##WORKING ON
		prediction_result_fin_list=glob.glob(prediction_result_fn+'*')
		pred_fin_list+=prediction_result_fin_list
	return pred_fin_list

def writePositivePrediction(output_path, IC50_cutoff, splicing_event_type, retrieve_prioritized):
	print '[INFO] Collecting binding predictions:'+output_path.rstrip('/').split('/')[-1]
	pred_fin_list=[]
	if retrieve_prioritized:
		if glob.glob(output_path+'/tmp/pred/*/*')!=[] and glob.glob('/'.join(output_path.split('/')[:-1])+'/'+splicing_event_type+'.tier1'+'/tmp/pred/*/*')==[]:
			print '[INFO] Tier1 comparison result is empty. Retrieving tier2&tier3 prediction.'
			pred_fin_list=glob.glob(output_path+'/tmp/pred/*/*')
		else:
			print '[INFO] Retrieving t prediction from tier1 comparison.'
			pred_fin_list= retrievePrioritizedfromPrimary(output_path, splicing_event_type)
	else:
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

def buildJunctionKmer(seq, kmer_length, kmer_dict, name, keeper, anno):
	if anno in keeper:
		return kmer_dict

	else:
		keeper[anno]=''
		for k in kmer_length:
			for i in range(len(seq)-k+1):
				s=seq[i:i+k]
				if s not in kmer_dict:
					kmer_dict[s]=[]
				kmer_dict[s].append(name)
		return kmer_dict 

def epitopeUniquenessAnnotation(kmer_length, screening_dir):
	kmer={}
	keeper={}
	peptide_fin_dir=screening_dir+'/*.*/tmp/prot/*'
	peptide_fin_list=glob.glob(peptide_fin_dir)
	print '[INFO] Total pool of splice junctions:',len(peptide_fin_list)
	tot=len(peptide_fin_list)-1
	for m, peptide_fin in enumerate(peptide_fin_list):
		config.update_progress(m/(0.0+tot))
		pep_gene_name=peptide_fin.split('/')[-1].split('_')[0]
		anno=''
		for l in open(peptide_fin):
			if l[0]=='>':
				anno=l.strip()
				continue
			ls=l.strip().upper()
			kmer=buildJunctionKmer(ls, kmer_length, kmer, pep_gene_name, keeper, anno)
			anno=''
	print '[INFO] Total number of AS junction k mers from the analyzed sequencing data:',len(kmer)
	return kmer


def loadKmer(kmer_file):
	kmer={}
	for l in open(kmer_file):
		ls=l.strip().split('\t')
		kmer[ls[0]]=ls[1].strip(';')
	return kmer

def writePeptideSummary(output_path, screening_result_path, gene_exp_path, IC50_cutoff, sample_HLA, sample_list, match_normal, uniqueness, epitope_len_list, db_dir):
	pep_dict={}
	kmer_dict={}
	header_ext=[]
	if gene_exp_path:
		gene_exp_dict, gene_exp_header=loadGeneExp(gene_exp_path)#exp. or fpkm1. matrix
		header_ext+=gene_exp_header
	if match_normal:
		print '[INFO] Loading kmers for checking canonical proteome: '+','.join(map(str,epitope_len_list))
		for kmer_len in epitope_len_list:
			kmer_dict[kmer_len]=loadKmer(db_dir+'/resources/kmers/uniprot-all.fasta.'+str(kmer_len)+'mer_dict.txt')
		header_ext+=['canonical_match']
	if uniqueness:
		junction_kmer=epitopeUniquenessAnnotation(epitope_len_list, '/'.join(screening_result_path.split('/')[:-1]))
		header_ext+=['uniqueness']
	screening_result_dict, screening_result_header=loadScreening(screening_result_path)
	
	fout_screening=open(output_path+'/epitope_summary.peptide-based.txt','w')
	fout_screening.write('\t'.join(['epitope','as_event','junction_peptide_form','inclusion_form','num_hla','num_sample','hla_types']+sample_list+screening_result_header+header_ext)+'\n')


	for l in open(output_path+'/pred_filtered.score'+str(IC50_cutoff)+'.txt'):
		ls=l.strip().split('\t')
		peptide=ls[8]
		hla_type=ls[7]
		pred_score=ls[10]
		splicing=(ls[0]+'_'+ls[1]).replace('_',':')
		form=ls[2].split(':')[5][:3]
		inc_form=''
		if form=='inc':
			inc_form=ls[2].split(':')[5][:4]
		if peptide+'|'+splicing+form not in pep_dict:
			pep_dict[peptide+'|'+splicing+form+'|'+inc_form]=[]
		pep_dict[peptide+'|'+splicing+form+'|'+inc_form].append(hla_type+'|'+pred_score)
	for k in pep_dict:
		ks=k.split('|')
		peptide=ks[0]
		splicing=ks[1][:-3]
		form=ks[1][-3:]
		inc_form=ks[2]
		hla_info=pep_dict[k]
		hla={}
		for info in hla_info:
			info_split=info.split('|')
			hla_type=info_split[0]
			pred_score=info_split[1]
			hla[hla_type]=pred_score
		num_hla=len(hla)
		line= [peptide, splicing, form, inc_form, str(num_hla), ';'.join(sorted(hla.keys()))]
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
		optional_annotations=[]
		if gene_exp_path:
			gene_exp_list=gene_exp_dict[splicing.split(':')[0]]
			optional_annotations+=gene_exp_list
		if match_normal:
			matched_protein='NonCanonical' if peptide not in kmer_dict[len(peptide)] else kmer_dict[len(peptide)][peptide]
			optional_annotations+=[matched_protein]
		if uniqueness:
			if peptide not in junction_kmer:
				print 'error', ls[0]
			else:
				uniqueness_annotation='-'
				if len(junction_kmer[peptide])==1:
					uniqueness_annotation='unique:'+junction_kmer[peptide][0]
				else:
					uniqueness_annotation='multi:'+';'.join(junction_kmer[peptide])
				optional_annotations+=[uniqueness_annotation]
		fout_screening.write('\t'.join(line+screening_result_dict[splicing]+optional_annotations)+'\n') 


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
		if splicing+form not in junction_dict:
			junction_dict[splicing+form]=[]
		junction_dict[splicing+form].append(hla_type+'|'+pred_score)	
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
		optional_annotations=[]
		if gene_exp_path:
			gene_exp_list=gene_exp_dict[splicing.split(':')[0]]
			optional_annotations+=gene_exp_list
		fout_screening.write('\t'.join(line+screening_result_dict[splicing]+optional_annotations)+'\n')

def main(args):
	outdir=args.outdir
	splicing_event_type=args.splicing_event_type
	IC50_cutoff=args.ic50_cut_off
	gene_exp_matrix=args.gene_exp_matrix
	match_normal=True if args.no_match_to_canonical_proteome==False else False
	uniqueness=True if args.no_uniqueness_annotation==False else False
        prioritized_only=args.tier3_only
	keep_exist=args.keep_exist

	write_positive=True
	if keep_exist==True:
		if os.path.exists(outdir+'/'+splicing_event_type+'.tier2tier3/pred_filtered.score'+str(IC50_cutoff)+'.txt') or os.path.exists(outdir+'/'+splicing_event_type+'.tier1/pred_filtered.score'+str(IC50_cutoff)+'.txt'):
			write_positive=False

        sample_list,sample_HLA=loadSampleHLA(args.mhc_by_sample)
        analysis_name,db_dir=[l.strip() for l in open(args.parameter_fin)][0:2]
        db_dir='/'.join(db_dir.rstrip('/').split('/')[:-1])
        analysis_name=analysis_name+'.'+splicing_event_type
        screening_result_path=outdir+'/'+analysis_name
        if config.file_len(screening_result_path+'.tier1.txt')==1 and prioritized_only==False:
            prioritized_only=True
            print "[INFO] No tier1 comparisons (tissue-matched normal) found. Use tier2tier3 only mode. "

	epitope_len_list=map(int,args.epitope_len_list.split(','))
	if write_positive:
		if prioritized_only:
			writePositivePrediction(outdir+'/'+splicing_event_type+'.tier2tier3', IC50_cutoff, splicing_event_type, False)
		else:
			writePositivePrediction(outdir+'/'+splicing_event_type+'.tier1', IC50_cutoff, splicing_event_type, False)
			writePositivePrediction(outdir+'/'+splicing_event_type+'.tier2tier3', IC50_cutoff, splicing_event_type, True)
	else:
		print '[INFO] File(s) '+outdir+'/'+splicing_event_type+'.*/pred_filtered.score'+str(IC50_cutoff)+'.txt'+' exist. Skip the step generating this output.'
		
	if prioritized_only==False:
		writePeptideSummary(outdir+'/'+splicing_event_type+'.tier1', outdir+'/'+analysis_name+'.tier1.txt', gene_exp_matrix, IC50_cutoff, sample_HLA, sample_list, match_normal, uniqueness, epitope_len_list, db_dir)
	writePeptideSummary(outdir+'/'+splicing_event_type+'.tier2tier3', outdir+'/'+analysis_name+'.tier2tier3.txt', gene_exp_matrix, IC50_cutoff, sample_HLA, sample_list, match_normal, uniqueness, epitope_len_list, db_dir)
	
	if prioritized_only==False:
		writeJunctionSummary(outdir+'/'+splicing_event_type+'.tier1', outdir+'/'+analysis_name+'.tier1.txt', gene_exp_matrix, IC50_cutoff, sample_HLA, sample_list)
	writeJunctionSummary(outdir+'/'+splicing_event_type+'.tier2tier3', outdir+'/'+analysis_name+'.tier2tier3.txt', gene_exp_matrix, IC50_cutoff, sample_HLA, sample_list)
	

if __name__ == '__main__':
	main()
