import sys,glob
from . import config
#only confident HLAs are selected
#when multiple top hits coexist("'", or ambiguity), only the top one/or equally top one selected. #TODO: included all ambiguious HLAs)

def main(args):
	outdir=args.outdir.rstrip('/')
	fin_list=glob.glob(outdir+'/*/hla_types/hla_types-ClassI.HLAgenotype4digits')
	HLA2patients={}
	HLA_list=set()
	for fin in fin_list:
		name=fin.split('/')[-3]
		print name
		n=0
		if name in HLA2patients:
			print 'dup'
			exit()
		HLA2patients[name]=[]
		for l in open(fin):
			if n==0:
				n+=1
				continue
			ls=l.strip().split('\t')
			if float(ls[2])<=0.05:
				HLA2patients[name].append('HLA-'+ls[1].rstrip("'"))
				HLA_list.add('HLA-'+ls[1].rstrip("'"))
			if float(ls[4])<=0.05:
				HLA2patients[name].append('HLA-'+ls[3].rstrip("'"))
				HLA_list.add('HLA-'+ls[3].rstrip("'"))
	fout1=open(outdir+'/hla_patient.tsv','w')
	for k in HLA2patients:
		fout1.write('\t'.join([k]+HLA2patients[k])+'\n')
	fout1.close()

	HLA_list=list(HLA_list)
	fout2=open(outdir+'/hla_types.list','w')
	for h in sorted(HLA_list):
		fout2.write(h+'\n')
	fout2.close()

	fin_list2=glob.glob(outdir+'/*/hla_types/hla_types-ClassI.expression')
	HLAexp2patients={}
	for fin2 in fin_list2:
		name=fin2.split('/')[-3]
		print name
		if name in HLAexp2patients:
			print 'dup'
			exit()
		HLAexp2patients[name]=[]
		for l in open(fin2):
			ls=l.strip()
			HLAexp2patients[name].append(ls)

	fout3=open(outdir+'/hla_exp.list','w')
	for k in HLAexp2patients:
		fout3.write('\t'.join([k]+HLAexp2patients[k])+'\n')

if __name__ == '__main__':
	main()
