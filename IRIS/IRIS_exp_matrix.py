import sys, numpy, os

def main(args):
	gene_exp_file_list=args.gene_exp_file_list
	n=0
	exp_cutoff=args.exp_cutoff
	data_name=args.data_name
	outdir=args.outdir.rstrip('/')
	matrix_dic={}
	sample_order=[]
	os.system('mkdir -p '+outdir)
	for fin in open(gene_exp_file_list):
		fin=fin.strip()
		sample_order.append(fin.split('/')[-3])
		print '[INFO] Loading', fin
		for m,l in enumerate(open(fin)):
			if m==0:
				continue
			ls=l.strip().split('\t')
			if ls[0] not in matrix_dic:
				matrix_dic[ls[0]]={}
			sample_name=fin.split('/')[-3]
			matrix_dic[ls[0]][sample_name]=ls[9]

	fout=open(outdir+'/exp.merged_matrix.'+data_name+'.txt','w')
	fout.write('\t'.join(['geneName']+sample_order)+'\n')
	fout_filter=open(outdir+'/fpkm'+str(exp_cutoff)+'.exp.merged_matrix.'+data_name+'.txt','w')
	fout_filter.write('\t'.join(['geneName']+sample_order)+'\n')
	fout_summary=open(outdir+'/sum.exp.merged_matrix.'+data_name+'.txt','w')
	fout_summary.write('\t'.join(['geneName'])+'\tmean\tQ25\tQ75\n')
	for k in sorted(matrix_dic.keys()):
		line_list=[matrix_dic[k][s] for s in sample_order]
		exp_all=map(float,line_list)
		mean=numpy.mean(exp_all)
		q25,q75=numpy.percentile(exp_all, [25 ,75])
		fout.write(k+'\t'+'\t'.join(line_list)+'\n')
		fout_summary.write(k+'\t'+str(mean)+'\t'+str(q25)+'\t'+str(q75)+'\n')
		if mean>=exp_cutoff:
			fout_filter.write(k+'\t'+'\t'.join(line_list)+'\n')
		

if __name__ == '__main__':
	main()