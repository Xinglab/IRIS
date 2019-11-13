import sys,glob,os
import numpy as np

def main(args):
	dump_all=args.dump_all
	# dump_all=False
	# if len(sys.argv)==6:
	# 	dump_all=True
	ms_search_result=args.MS_search_result_prefix
	outdir=args.outdir.rstrip('/')
	os.system('mkdir -p '+outdir)
	Qval=args.MS_Qvalue#0.1
	fin_list=glob.glob(ms_search_result+'*.tsv')

	fout=open(outdir+'/filteredPSM_'+str(Qval)+'.'+ms_search_result.rstrip('/').split('/')[-1]+'.tsv','w')
	peptide={}
	for fin in fin_list:
		header=[]
		for i,l in enumerate(open(fin)):
			if i==0:
				header=l.strip().split('\t')
				fout.write(l)
				continue
			ls=l.strip().split('\t')
			line_dict=dict(zip(header,ls))
			if float(line_dict['QValue'])<=float(Qval):
				fout.write(l)
				if peptide.has_key(line_dict['Peptide'])!=True:
					peptide[line_dict['Peptide']]=set()
				peptide[line_dict['Peptide']].add(line_dict['Protein']+'|QValue:'+line_dict['QValue'])
	fout.close()

	fout_pep=open(outdir+'/filteredPEP_'+str(Qval)+'.'+ms_search_result.rstrip('/').split('/')[-1]+'.tsv','w')
	for k in peptide.keys():
		fout_pep.write('{}\t{}\n'.format(k,';'.join(list(peptide[k]))))
	fout_pep.close()

	pred_fin=args.binding_prediction
	validated_pep={}
	validated_pep_exp={}
	hla_list={}
	for l in open(pred_fin):
		ls=l.strip().split('\t')
		if peptide.has_key(ls[8]):
			if validated_pep.has_key(ls[8])!=True:
				validated_pep[ls[8]]=set()
				validated_pep_exp[ls[8]]=set()
			validated_pep[ls[8]].add(ls[7]+'|'+str(round(float(ls[10]),2)))
			form=ls[2].split(':')[5]
			if len(ls)>=12:
				exp_value=[ls[11],ls[14],ls[17]]
				validated_pep_exp[ls[8]].add(form+'|'+':'.join(exp_value))
			hla_list[ls[7]]=''
		else:
			if dump_all:
				if validated_pep.has_key(ls[8])!=True:
					validated_pep[ls[8]]=set()
					validated_pep_exp[ls[8]]=set()
				validated_pep[ls[8]].add(ls[7]+'|'+str(round(float(ls[10]),2)))
				exp_value=[ls[11],ls[14],ls[17]]
				form=ls[2].split(':')[5]
				validated_pep_exp[ls[8]].add(form+'|'+':'.join(exp_value))
				hla_list[ls[7]]=''
	fout_pred_match_name= outdir+'/pred_match.filteredPEP_'+str(Qval)+'.'+ms_search_result.rstrip('/').split('/')[-1]+'.tsv'

	if dump_all:
		fout_pred_match_name=outdir+'/pred_annotated.filteredPEP_'+str(Qval)+'.'+ms_search_result.rstrip('/').split('/')[-1]+'.tsv'
	fout_pred_match=open(fout_pred_match_name,'w')
	for pep in sorted(validated_pep):
		MS_info=''
		MS_sum=''
		HLA_info=';'.join(list(validated_pep[pep]))
		HLA_sum= min(map(lambda x:float(x.split('|')[1]),validated_pep[pep]))
		if peptide.has_key(pep):
			MS_info=';'.join(list(peptide[pep]))
			MS_sum=min(map(lambda x:float(x.split('|')[-1].split(':')[1]),peptide[pep]))
		else:
			if dump_all:	
				MS_info='No_sig_MS_hit'
				MS_sum='NA'
			else:
				exit('error. Peptide not in MS dictionary')
		if validated_pep_exp!={}:
			Exp_info=';'.join(validated_pep_exp[pep])
			Exp_sum= np.mean(map(lambda x:float(x.split('|')[1].split(':')[2]),validated_pep_exp[pep])) 
			fout_pred_match.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(pep, HLA_info, HLA_sum, MS_info, MS_sum, Exp_info, Exp_sum))
		else:
			fout_pred_match.write('{}\t{}\t{}\t{}\t{}\n'.format(pep, HLA_info, HLA_sum, MS_info, MS_sum))
	fout_pred_match.close()

	if dump_all:
		exit()
	for hla in hla_list.keys():
		fout_pred_match_hla=open(outdir+'/pred_match.'+hla+'.filteredPEP_'+str(Qval)+'.'+ms_search_result.rstrip('/').split('/')[-1]+'.tsv','w')
		for pep in sorted(validated_pep):
			has_hla_type=False
			for hla_binding in validated_pep[pep]:
				if hla_binding.find(hla)!=-1:
					has_hla_type=True
					break
			if has_hla_type:
				fout_pred_match_hla.write('{}\t{}\t{}\t{}\n'.format(pep, ';'.join(list(validated_pep[pep])),';'.join(list(peptide[pep])),';'.join(validated_pep_exp[pep])))
		fout_pred_match_hla.close()

if __name__ == '__main__':
	main()

