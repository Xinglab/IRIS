import sys,os
from . import config

def main(args):
	MS_fin=args.MS_fin#LCL.mzML/Subject3_rep3_021213_Fx100mM.mzML
	MS_db=args.MS_db #~/flashscratch/LCL/all_jp/tmp/proteome_ref_merged.fa
	outdir=args.outdir
	cmd1= args.java_path+' -Xmx8g -jar '+args.MSGF_path+' -s '+MS_fin+' -d '+MS_db+' -e 0 -tda 1 -maxLength 13 -minLength 7 -inst 1 -t 5ppm -o '+outdir
	print cmd1
	os.system(cmd1)
	cmd2= args.java_path+' -Xmx8g -cp '+args.MSGF_path+' edu.ucsd.msjava.ui.MzIDToTsv -i '+outdir

	print cmd2
	os.system(cmd2)

if __name__ == '__main__':
	main()