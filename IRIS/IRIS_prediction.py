import sys, argparse, os ,datetime,logging, uuid, glob
from . import config

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

def selectJC(AS_coord,deltaPSI_c2n,cut_off, select_all):
	if select_all:
		return [(AS_coord[2], AS_coord[3]),(AS_coord[2],AS_coord[0]),(AS_coord[1],AS_coord[3])]
	if float(deltaPSI_c2n)<cut_off:# tumor skipping
		return [(AS_coord[2], AS_coord[3])]
	else:
		return [(AS_coord[2],AS_coord[0]),(AS_coord[1],AS_coord[3])]

def AS2FT(AS_chrom, AS_coord,AS_direction, deltaPSI_c2n, cuf_off, select_all, extracelllularDict, id_line, fout_dict):
	JC_coord_list=selectJC(AS_coord,deltaPSI_c2n,cuf_off, select_all)
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

def extracellularAnnotation(screening_result_fin, outdir, extracelllularDict, deltaPSI_cut_off, select_all):
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
		JC_pep_fasta=AS2FT(des[2],des[4:8],des[3],ls[deltaPSI_column],deltaPSI_cut_off,select_all,extracelllularDict,ls[0], fout_dict)
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
	fout2=open(screening_result_fin+'.ExtraCellularAS.txt','w')
	fout2.write('{}\t{}\t{}\t{}\n'.format('as_event','protein_domain_loc','protein_domain_loc_by_as_exon','\t'.join(screening_result_dict['header'])))
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
			screen_print='\t'.join(screening_result_dict[k])
		fout2.write(k+'\t'+line.rstrip(';')+'\t'+screen_print+'\n')
	fout2.close()
	os.system('rm '+screening_result_fin+'.CellSurfAnno.tmp')


def epitopePredictionPrep(outdir, hla_list_fin, analysis_name, iedb_path):
	fin_list=glob.glob(outdir+'/tmp/prot.compared/skp/*.fa')
	fin_list=fin_list+glob.glob(outdir+'/tmp/prot.compared/inc/*.fa')
	hla_list=[]
	num=90
	for l in open(hla_list_fin):
		ls=l.strip()
		hla_list.append(ls)
	print "[INFO] Total HLA types loaded:", len(hla_list), ". Total peptide splice junctions loaded:",len(fin_list)

	analysis_name=outdir.split('/')[-1]
	list_name='cmdlist.pep2epitope_'+analysis_name
	fout_list=open(list_name,'w')

	m=0
	for i in xrange(0,len(hla_list),3):
		hla_types=','.join(hla_list[i:i+3])
		n=0
		line=''
		for fin in fin_list:
			if os.stat(fin).st_size != 0:
				n+=1
				#line+='echo run '+fin+'\npython '+IRIS_PACKAGE_PATH+'/IRIS/IRIS.pep2epitope.py '+fin+' --hla-allele-list '+hla_types+' -o '+outdir+'\n'
				line+='echo run '+fin+'\nIRIS pep2epitope '+fin+' --hla-allele-list '+hla_types+' -o '+outdir+' --iedb-local '+iedb_path+'\n'
				line+='sleep 70\n'
				if n%num==0:
					submission_file_name=outdir+'/tmp/submit.IRIS_pep2epitope.py.'+str(n)+'.'+hla_types+'.sh'
					fout_list.write(submission_file_name+'\n')
					m+=1
					fout=open(submission_file_name,'w')
					fout.write(line)
					fout.close()
					line=''
		if n%num!=0:
			# print line
			submission_file_name=outdir+'/tmp/submit.IRIS_pep2epitope.py.'+str(n)+'.'+hla_types+'.sh'
			fout_list.write(submission_file_name+'\n')
			m+=1
			fout=open(submission_file_name,'w')
			fout=open(outdir+'/tmp/submit.IRIS_pep2epitope.py.'+str(n)+'.'+hla_types+'.sh','w')
			fout.write(line)
			fout.close()
			line=''

	fout_list.close()
	fout_qsub=open('qsub.IRIS_pep2epitope.py.'+analysis_name+'.sh','w')

	cmd='qsub -t 1-'+str(m)+':1 qsub.IRIS_pep2epitope.py.'+analysis_name+'.sh'

	fout_qsub.write('#!/bin/bash\n#$ -N IRIS_pep2epitope\n#$ -S /bin/bash\n#$ -R y\n#$ -l '+config.QSUB_PREDICTION_CONFIG+'\n#$ -V\n#$ -cwd\n#$ -j y\n#$ -m bea\n')
	fout_qsub.write('export s=`sed -n ${SGE_TASK_ID}p '+list_name+'`\necho $s\nbash $s')
	fout_qsub.close()
	print cmd

def main(args):
	#Define parameters
	IRIS_screening_result=args.IRIS_screening_result_path
	deltaPSI_cut_off=float(args.deltaPSI_cut_off)
	select_all= True if args.extracellular_anno_by_junction==False else False
	analysis_name=[l.strip() for l in open(args.parameter_fin)][0]
	hla_list_fin=args.mhc_list
	iedb_path=args.iedb_local
	extracelllularDict=loadFeatures(config.EXTRACELLULAR_FEATURES_UNIPROT2GTF_MAP_PATH) #IRIS_package_dir.rstrip('/')+'/IRIS/data/features.uniprot2gtf.ExtraCell.txt'
	print "[INFO] Total extracellular annotated loaded:",len(extracelllularDict)

	extracellularAnnotation(IRIS_screening_result+'/'+analysis_name+'.primary.txt', IRIS_screening_result, extracelllularDict, deltaPSI_cut_off, select_all)
	epitopePredictionPrep(IRIS_screening_result+'/primary',hla_list_fin, analysis_name, iedb_path)

	
	extracellularAnnotation(IRIS_screening_result+'/'+analysis_name+'.prioritized.txt', IRIS_screening_result, extracelllularDict, deltaPSI_cut_off, select_all)
	epitopePredictionPrep(IRIS_screening_result+'/prioritized',hla_list_fin, analysis_name, iedb_path)

if __name__ == '__main__':
	main()



