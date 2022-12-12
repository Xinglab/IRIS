import sys, argparse, os ,datetime,logging, uuid, glob
from . import config
import numpy as np

ID = str(uuid.uuid4()).split('-')[0]


def loadFeatures(fin):
	extracelllularDict={}
	for l in open(fin):
		ls=l.strip().split('\t')
		AS_s_e=ls[0].split('-')
		if AS_s_e[0] not in extracelllularDict:
			extracelllularDict[AS_s_e[0]]=[]
		extracelllularDict[AS_s_e[0]].append((AS_s_e[1],[ls[2]]+[ls[5]]+[ls[3]]+ls[6:9]))
		if AS_s_e[1] not in extracelllularDict:
			extracelllularDict[AS_s_e[1]]=[]
		extracelllularDict[AS_s_e[1]].append((AS_s_e[0],[ls[2]]+[ls[5]]+[ls[3]]+ls[6:9]))

	return extracelllularDict

def selectJunction(AS_coord, deltaPSI_c2n, cut_off, if_select_all, splicing_event_type): 
	if splicing_event_type == 'SE':
		skp = (AS_coord[2], AS_coord[3],'skp')
		inc1 = (AS_coord[2],AS_coord[0],'inc1')
		inc2 = (AS_coord[1],AS_coord[3],'inc2')
	elif splicing_event_type == 'A3SS':
		skp = (AS_coord[5], AS_coord[2],'skp')
		inc1 = (AS_coord[5],AS_coord[0],'inc1')
		inc2 = (AS_coord[2],AS_coord[2],'inc2')
	elif splicing_event_type == 'A5SS':
		skp = (AS_coord[3], AS_coord[4],'skp')
		inc1 = (AS_coord[3],AS_coord[3],'inc1')
		inc2 = (AS_coord[1],AS_coord[4],'inc2')
	elif splicing_event_type == 'RI':
		skp = (AS_coord[3], AS_coord[4],'skp')
		inc1 = (AS_coord[3],AS_coord[3],'inc1')
		inc2 = (AS_coord[4],AS_coord[4],'inc2')

	if if_select_all:
		return [skp, inc1, inc2]
	if float(deltaPSI_c2n) < cut_off:  # tumor skipping
		return [skp]
	else:
		return [inc1, inc2]
# def selectJC(AS_coord,deltaPSI_c2n,cut_off, select_all):
# 	if select_all:
# 		return [(AS_coord[2], AS_coord[3]),(AS_coord[2],AS_coord[0]),(AS_coord[1],AS_coord[3])]
# 	if float(deltaPSI_c2n)<cut_off:# tumor skipping
# 		return [(AS_coord[2], AS_coord[3])]
# 	else:
# 		return [(AS_coord[2],AS_coord[0]),(AS_coord[1],AS_coord[3])]

def AS2FT(AS_chrom, AS_coord,AS_direction, deltaPSI_c2n, cuf_off, select_all, extracelllularDict, id_line, fout_dict, splicing_event_type):
	JC_coord_list=selectJunction(AS_coord,deltaPSI_c2n,cuf_off, select_all, splicing_event_type)
	for JC_coord in JC_coord_list:
		if JC_coord[0] in extracelllularDict:
			for exon_info in extracelllularDict[JC_coord[0]]:
				if exon_info[0].startswith(AS_chrom+':'):
					prot_name,feature,isopep_pos,feat_pos,left,right=exon_info[1]
					if feature.startswith('TOPO_DOM:Extracellular'):
						fout_dict['\t'.join([id_line,AS_chrom, JC_coord[0], prot_name,isopep_pos,feat_pos,left, right])]=''
		if AS_chrom+':'+JC_coord[1] in extracelllularDict:
			for exon_info in extracelllularDict[AS_chrom+':'+JC_coord[1]]:
				prot_name,feature,isopep_pos,feat_pos,left,right=exon_info[1]
				if feature.startswith('TOPO_DOM:Extracellular'):
					fout_dict['\t'.join([id_line,AS_chrom, JC_coord[1], prot_name,isopep_pos,feat_pos,left, right])]=''
	return JC_coord_list
def loadScreening(screening_result):
	screening_result_dict={}
	for n,l in enumerate(open(screening_result)):
		ls=l.strip().split('\t')
		if n==0:
			screening_result_dict['header']=ls[1:]
			continue
		reformat_name=ls[0].split(':')
		reformat_name=':'.join([reformat_name[0],reformat_name[1]]+reformat_name[2:]) #.split('.')[0]
		screening_result_dict[reformat_name]=ls[1:]
	return screening_result_dict

def retriveJunctionPeptide(k, screening_result_fin, splicing_event_type, pep_dir_prefix):
	IRIS_screening_result='/'.join(screening_result_fin.split('/')[:-1])
	screening_tier=screening_result_fin.split('/')[-1].split('.')[-2]

	peptide_file_name=k.replace(':','_').replace('_chr','.chr')+'.prot.fa'
	peptide_file_name_full_skp=IRIS_screening_result.rstrip('/')+'/'+splicing_event_type+'.'+screening_tier+'/tmp/'+pep_dir_prefix+'.compared/skp/skp.'+peptide_file_name
	peptide_file_name_full_inc=IRIS_screening_result.rstrip('/')+'/'+splicing_event_type+'.'+screening_tier+'/tmp/'+pep_dir_prefix+'.compared/inc/inc.'+peptide_file_name
	junction_peptide=[]
	if os.path.exists(peptide_file_name_full_skp):
		for l in open(peptide_file_name_full_skp):
			if l.startswith('>'):
				continue
			junction_peptide.append(l.strip())
	if os.path.exists(peptide_file_name_full_inc):
		for l in open(peptide_file_name_full_inc):
			if l.startswith('>'):
				continue
			junction_peptide.append(l.strip())		
	junction_peptide=';'.join(junction_peptide)
	junction_peptide='-' if junction_peptide=='' else junction_peptide
	return junction_peptide

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

def extracellularAnnotation(screening_result_fin, splicing_event_type, extracelllularDict, deltaPSI_cut_off, select_all, gene_exp_path, pep_dir_prefix):
	if select_all:
		deltaPSI_column=0
	if select_all==False:
		try:
			deltaPSI_column=int(args.deltaPSI_column)-1
		except:
			exit('[Error] choose deltaPSI column or use -u to report all Junctions.')

	screening_result_dict=loadScreening(screening_result_fin)
	#Select AS version, region and frame, translate to peptides.
	fout_dict={}
	fout=open(screening_result_fin+'.CellSurfAnno.tmp','w')
	for n,l in enumerate(open(screening_result_fin)):
		if n==0:
			continue
		ls=l.strip().split('\t')
		des=ls[0].split(':')
		JC_coord_list=AS2FT(des[2],des[4:],des[3],ls[deltaPSI_column],deltaPSI_cut_off,select_all,extracelllularDict,ls[0], fout_dict, splicing_event_type)
	for k in fout_dict.keys():
		fout.write(k+'\n')
	fout.close()

	extracellular_AS={}
	for l in open(screening_result_fin+'.CellSurfAnno.tmp'):
		ls=l.strip().split('\t')
		if ls[0] not in extracellular_AS:
			extracellular_AS[ls[0]]={}
		if ls[3]+':'+ls[5] not in extracellular_AS[ls[0]]:
			extracellular_AS[ls[0]][ls[3]+':'+ls[5]]=[]
		extracellular_AS[ls[0]][ls[3]+':'+ls[5]].append(ls[1]+':'+ls[2]+':'+ls[4])
	
	if gene_exp_path:
		gene_exp_dict, gene_exp_header=loadGeneExp(gene_exp_path)

	fout2=open(screening_result_fin+'.ExtraCellularAS.txt','w')
	if gene_exp_path:
		fout2.write('{}\t{}\t{}\t{}\t{}\t{}\n'.format('as_event','protein_domain_loc','protein_domain_loc_by_as_exon','\t'.join(screening_result_dict['header']),'\t'.join(gene_exp_header), 'junction_peptide'))
	else:
		fout2.write('{}\t{}\t{}\t{}\t{}\n'.format('as_event','protein_domain_loc','protein_domain_loc_by_as_exon','\t'.join(screening_result_dict['header']),'junction_peptide'))
	for k in sorted(extracellular_AS):
		line=''
		line=';'.join(extracellular_AS[k].keys())+'\t'
		for anno in sorted(extracellular_AS[k]):
			anno_line=anno+'<-'
			for exon in extracellular_AS[k][anno]:
				anno_line+=exon+'|'
			anno_line=anno_line.rstrip('|')
			line+=anno_line+';'
		screen_print='NA'
		if k in screening_result_dict:
			junction_peptide=retriveJunctionPeptide(k, screening_result_fin, splicing_event_type, pep_dir_prefix)
			screen_print='\t'.join(screening_result_dict[k])
		optional_annotations=''
		if gene_exp_path:
			gene_exp_list=gene_exp_dict[k.split(':')[0]]
			optional_annotations='\t'.join(gene_exp_list)
		fout2.write('\t'.join([k,line.rstrip(';'), screen_print, optional_annotations, junction_peptide])+'\n')
	fout2.close()
	os.system('rm '+screening_result_fin+'.CellSurfAnno.tmp')


def epitopePredictionPrep(outdir, hla_list_fin, analysis_name, iedb_path, epitope_len_list, task_dir, pep_dir_prefix):
	fin_list=glob.glob(outdir+'/tmp/'+pep_dir_prefix+'.compared/skp/*.fa')
	fin_list=fin_list+glob.glob(outdir+'/tmp/'+pep_dir_prefix+'.compared/inc/*.fa')
	hla_list=[]
	num=90
	for l in open(hla_list_fin):
		ls=l.strip()
		hla_list.append(ls)
	print "[INFO] Total HLA types loaded:", len(hla_list), ". Total peptide splice junctions loaded:",len(fin_list)

	analysis_name=outdir.split('/')[-1]

	script_count = 0
	def write_task_script(script_contents):
		task_script_base = 'pep2epitope_{}.{}.sh'.format(analysis_name, script_count)
		task_script = os.path.join(task_dir, task_script_base)

		with open(task_script, 'w') as f_h:
			f_h.write('#!/bin/bash\n')
			f_h.write(script_contents)

	for i in xrange(0,len(hla_list),3):
		hla_types=','.join(hla_list[i:i+3])
		n=0
		script_contents = ''
		for fin in fin_list:
			if os.stat(fin).st_size != 0:
				n+=1
				script_contents += 'echo run '+fin+'\nIRIS pep2epitope '+fin+' --hla-allele-list '+hla_types+' -o '+outdir+' --iedb-local '+iedb_path+' -e '+epitope_len_list+'\n'
				if n%num==0:
					write_task_script(script_contents)
					script_count += 1
					script_contents = ''
		if n%num!=0:
			write_task_script(script_contents)
			script_count += 1
			script_contents = ''


def main(args):
	#Define parameters
	IRIS_screening_result=args.IRIS_screening_result_path.rstrip('/') #Modified 2021
	deltaPSI_cut_off=float(args.deltaPSI_cut_off)
	splicing_event_type=args.splicing_event_type
	select_all= True if args.extracellular_anno_by_junction==False else False
	analysis_name=[l.strip() for l in open(args.parameter_fin)][0]
	analysis_name=analysis_name+'.'+splicing_event_type #group_name.type
	prioritized_only=args.tier3_only
	extracellular_only=args.extracellular_only
	gene_exp_matrix=args.gene_exp_matrix
	task_dir=args.task_dir
	all_orf=args.all_orf
	pep_dir_prefix='prot'
	if all_orf:
		pep_dir_prefix='prot_allorf'
	if not os.path.exists(task_dir):
		os.makedirs(task_dir)

	if extracellular_only==False:
		hla_list_fin=args.mhc_list
		iedb_path=args.iedb_local
		epitope_len_list=args.epitope_len_list.rstrip(',')
	

	extracelllularDict=loadFeatures(config.EXTRACELLULAR_FEATURES_UNIPROT2GTF_MAP_PATH) 
	print "[INFO] Total extracellular annotation loaded:",len(extracelllularDict)

        if config.file_len(IRIS_screening_result+'/'+analysis_name+'.tier1.txt')==1 and prioritized_only==False:
            prioritized_only=True
            print "[INFO] No tier1 comparisons (tissue-matched normal) found. Use tier2&tier3 only mode. "
	if prioritized_only:
		extracellularAnnotation(IRIS_screening_result+'/'+analysis_name+'.tier2tier3.txt', splicing_event_type, extracelllularDict, deltaPSI_cut_off, select_all, gene_exp_matrix, pep_dir_prefix)
		if extracellular_only==False:
			epitopePredictionPrep(IRIS_screening_result+'/'+splicing_event_type+'.tier2tier3',hla_list_fin, analysis_name, iedb_path, epitope_len_list, task_dir, pep_dir_prefix)
	else:
		extracellularAnnotation(IRIS_screening_result+'/'+analysis_name+'.tier1.txt', splicing_event_type, extracelllularDict, deltaPSI_cut_off, select_all, gene_exp_matrix, pep_dir_prefix)
		extracellularAnnotation(IRIS_screening_result+'/'+analysis_name+'.tier2tier3.txt', splicing_event_type, extracelllularDict, deltaPSI_cut_off, select_all, gene_exp_matrix, pep_dir_prefix)
		if extracellular_only==False:
			epitopePredictionPrep(IRIS_screening_result+'/'+splicing_event_type+'.tier1',hla_list_fin, analysis_name, iedb_path, epitope_len_list, task_dir, pep_dir_prefix)
			#epitopePredictionPrep(IRIS_screening_result+'/'+splicing_event_type+'.prioritized',hla_list_fin, analysis_name, iedb_path, epitope_len_list, task_dir) #Only primary is needed as priortized can be parsed from it
	
if __name__ == '__main__':
	main()



