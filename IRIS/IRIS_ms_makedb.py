import os, glob
from . import config
# Out of ~21000 mapper list, 19188 ENSG name exists in 60,000 gtf files, FPKM>0 genes 14,300 can be mapped out of 22,300, FPKM=0 2189 can be mapped out of 27,700.
def non_exp_gene_list(fin_list):
	non_exp_dict={}
	total_num=0
	for fin in open(fin_list):
		fin=fin.strip()
		print fin
		i=0
		total_num+=1
		for l in open(fin):
			if i==0:
				i+=1
				continue
			ls=l.strip().split('\t')
			if float(ls[9])==0:#FPKM
				ensg=ls[0].split('.')[0]
				if non_exp_dict.has_key(ensg)!=True:
					non_exp_dict[ensg]=[]
				non_exp_dict[ensg].append(ls[4]) ##gene name, gene ac
	non_exp_dict_concensus={}
	for k in non_exp_dict:
		if len(non_exp_dict[k])==total_num:
			non_exp_dict_concensus[k]=non_exp_dict[k]
	print len(non_exp_dict_concensus), len(non_exp_dict)
	return non_exp_dict_concensus

def map2Uniprot(fin):
	ID_mapper={}
	#ID_mapper2={}
	for l in open(fin):
		ls=l.strip().split('\t')
		if len(ls)==2:
			ID_mapper[ls[0]]=ls[1]
	#	ID_mapper2[ls[1]]=ls[0]
	return ID_mapper

def makeProteoTranscriptomeRef(fin_list, uniprot_fasta, outdir):
	non_exp_dict_concensus=non_exp_gene_list(fin_list)
	ID_mapper=map2Uniprot(config.UNIPROT_ENSG_ID_MAP_PATH)
	c=0
	print '[INFO] Total uniprot mapped to ENSG',len(ID_mapper)
	write=True
	outdir=outdir.rstrip('/')
	fout=open(outdir+'/'+'exp_only.'+uniprot_fasta.split('/')[-1],'w')
	for l in open(uniprot_fasta):
		if l.startswith('>'):
			write=True
			ls=l.strip().split('|')
			uniprot_ac=ls[1].split('-')[0]
			if ID_mapper.has_key(uniprot_ac)==True:
				ensg_ids=ID_mapper[uniprot_ac].split(';')
				for ensg_id in ensg_ids:
					if non_exp_dict_concensus.has_key(ensg_id):
						c+=1
						write=False
			if write:
				fout.write(l)
		else:
			if write:
				fout.write(l)
	fout.close()
	print 'removed:' ,c

def main(args):
	outdir=args.outdir.rstrip('/')
	exp_fin_list=args.exp_fin_list
	uniprot_fasta=args.uniprot_fasta
	java_path = args.java_path
	MSGF_path = args.MSGF_path
	
	print '##Creating ProteoTransicritomic Ref'
	makeProteoTranscriptomeRef(exp_fin_list, uniprot_fasta, outdir)

	print '##Creating junction peptides db'
	fout=open(outdir+'/tmp/proteome_ref_junctions.fa','w')
	fin_list=glob.glob(outdir+'/tmp/prot.compared/skp/*')
	fin_list+=glob.glob(outdir+'/tmp/prot.compared/inc/*')
	print len(fin_list)
	for fin in fin_list:
		for i,l in enumerate(open(fin)):
			if i%2==0:
				fout.write(l)
			else:
				fout.write(l.upper())

	fout.close()

	input_proteome_db=outdir+'/'+'exp_only.'+uniprot_fasta.split('/')[-1]#full path

	cmd1='cat '+outdir+'/tmp/proteome_ref_junctions.fa '+input_proteome_db+' > '+outdir+'/tmp/proteome_ref_combined.fa'
	print '##Combining proteome db with junction peptides'
	print cmd1
	os.system(cmd1)

	cmd2=java_path+' -Xmx8g -cp '+MSGF_path+' edu.ucsd.msjava.msdbsearch.BuildSA -d '+outdir+'/tmp/proteome_ref_combined.fa'
	print '##Indexing the proteotranscriptomic db'
	print cmd2
	os.system(cmd2)

if __name__ == '__main__':
	main()