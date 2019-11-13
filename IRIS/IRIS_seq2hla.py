import sys, os, csv, argparse, logging, datetime
from . import config
#from utilities.seq2hla import seq2HLA

def run_seq2HLA(readsFilesCaseRNA,runname,bindir):
	if os.path.exists(runname+'/hla_types-ClassI.HLAgenotype4digits')==False:
		readsFiles_split=readsFilesCaseRNA.split(',')
		os.system('mkdir -p '+runname)
		#seq2HLA.main(runname,readsFiles_split[0],readsFiles_split[1])
		cmd = 'python '+bindir+'/seq2HLA.py -1 '+readsFiles_split[0]+' -2 '+readsFiles_split[1]+' -r '+runname+'/hla_types'
		os.system (cmd)
		print cmd
		if os.path.exists(runname+'/hla_types-ClassI.HLAgenotype4digits')==False:
			sys.exit('[seq2hla] # An Error Has Occured. seq2hla Incomplete. Exit!')
	else:
		print '[seq2hla] # Skipped seq2HLA.'

	HLA_type=[]
	for n,l in enumerate(open(runname+'/hla_types-ClassI.HLAgenotype4digits')):
		if n==0:
			continue
		ls=l.strip().split('\t')
		#print ls
		if ls[2]!='NA':
			if float(ls[2])<=0.05:
				HLA_type.append('HLA-'+ls[1].strip('\''))
		if ls[4]!='NA':
			if float(ls[4])<=0.05:
				HLA_type.append('HLA-'+ls[3].strip('\''))
		continue
	# for n,l in enumerate(open(runname+'-ClassII.HLAgenotype4digits')):
	# 	if n==0:
	# 		continue
	# 	ls=l.strip().split('\t')
	# 	if ls[2]!='NA':
	# 		if float(l[2])<=0.05:
	# 			HLA_type.append('HLA-'+l[1].strip('\''))
	# 	if ls[4]!='NA':
	# 		if float(ls[4])<=0.05:
	# 			HLA_type.append('HLA-'+l[3].strip('\''))
	# 	continue

	if len(HLA_type)==0:
		sys.exit('# [INFO] No HLA type predicted. Exit.')

	HLA_type_str=','.join(list(set(HLA_type)))
	return HLA_type_str	

def main(args):
	os.system('mkdir -p '+args.sampleID_outdir)
	sampleID = args.sampleID_outdir.rstrip('/')
	runname = sampleID+'/hla_types'
	bindir = args.seq2hla_path.rstrip('/')
	print '[INFO] # Start HLA typing.'

	hla=run_seq2HLA(args.readsFilesCaseRNA, runname, bindir)

	print '[INFO] # Completed. HLA types: '+hla

if __name__ == '__main__':
	main()
