import sys,os
from . import config

def main(args):
	MS_fin=args.MS_fin#LCL.mzML/Subject3_rep3_021213_Fx100mM.mzML
	MS_db=args.MS_db #~/flashscratch/LCL/all_jp/tmp/proteome_ref_merged.fa
        jvm_max_mem = args.jvm_max_mem #Xmx8g
        if jvm_max_mem!='':
            jvm_max_mem='-'+jvm_max_mem
        MSGF_search_para_fin = args.MSGF_search_para_fin
        if MSGF_search_para_fin!='': #allow user to provide no parameters by giving an empty file
            MSGF_search_para=[x.strip() for i,x in enumerate(open(MSGF_search_para_fin)) if i==0 ]
            if MSGF_search_para==[]:
                MSGF_search_para=''
            else:
                MSGF_search_para=MSGF_search_para[0]
        else:
            MSGF_search_para='-e 0 -tda 1 -maxLength 13 -minLength 7 -inst 1 -t 5ppm'
	outdir=args.outdir
	cmd1= args.java_path+' '+jvm_max_mem+' -jar '+args.MSGF_path+' -s '+MS_fin+' -d '+MS_db+' '+MSGF_search_para+' -o '+outdir
	print cmd1
	os.system(cmd1)
	cmd2= args.java_path+' '+jvm_max_mem+' -cp '+args.MSGF_path+' edu.ucsd.msjava.ui.MzIDToTsv -i '+outdir

	print cmd2
	os.system(cmd2)

if __name__ == '__main__':
	main()
