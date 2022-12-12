import numpy as np
import os,sys,glob,argparse
from scipy import stats
import statsmodels.stats.weightstats as smw
import matplotlib
matplotlib.use('agg')
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from . import config

def fileLength(fin):
	i=sum(1 for l in open(fin))
	return i

# def indiviualPlot(psi_df, event_id, panel_list):
# 	plt.figure(figsize=(15,9))
# 	sns_plot=sns.violinplot(data=psi_df, inner="box",cut=0)
# 	#sns_plot.set_yticks(np.arange(0,1,0.2))	
# 	sns_plot.set(ylim=(0, 1))
# 	sns.despine(offset=10, trim=True)
# 	sns_plot.set_ylabel('Percent-Spliced-In',fontweight='bold',fontsize=14)
# 	sns_plot.set_xticklabels(panel_list,rotation=20,ha='right',fontsize=12)
# 	sns_plot.set_title(event_id,fontsize=15,fontweight='bold')
# 	sns_plot.figure.savefig(event_id+".pdf")

def selectJC(event_id, tumor_form, tumor_form_cutoff, splicing_event_type):
	coord=event_id.split(':')
	if splicing_event_type=='SE':
		inc1=coord[2]+':'+str(int(coord[6])+1)+':'+coord[4]
		inc2=coord[2]+':'+str(int(coord[5])+1)+':'+coord[7]
		skp=coord[2]+':'+str(int(coord[6])+1)+':'+coord[7]
		if tumor_form>tumor_form_cutoff:#inc
			return [inc1, inc2]
		else:
			return [skp]
	else:
		exit('[Error] The splicine event type is not supported currently. Exiting.')

def loadJCvalue4events(fin_plot_query, deltaPSI_col, tumor_form_cutoff, jc_result, splicing_event_type, has_header):
	JC_dict={}
	for i,l in enumerate(open(fin_plot_query)):
		if i==0 and has_header:
			continue
		ls=l.strip().split('\t')
		event_id=ls[0]
		tumor_form=float(ls[deltaPSI_col]) # negative means skipping, pos means inc
		JC_select=selectJC(event_id, tumor_form, tumor_form_cutoff, splicing_event_type)
		for k in JC_select:
			JC_dict[k]=''
	for i,l in enumerate(open(jc_result)):
		if i==0:
			continue
		ls=l.strip().split('\t')
		if ls[0] in JC_dict:
			JC_dict[ls[0]]=l
	return JC_dict

def loadParametersRow(filter_para, panel_list):
	if filter_para.strip()!='':
		para, filter_panel_list=filter_para.split(' ')
		filter_cutoff_pval, filter_cutoff_dpsi, filter_cutoff_foc, filter_cutoff_pval_PT, filter_group_cutoff =para.split(',')
		filter_cutoff_pval_PT=float(filter_cutoff_pval_PT)
		filter_cutoff_dpsi=float(filter_cutoff_dpsi)
		filter_cutoff_foc=float(filter_cutoff_foc)
		filter_group_cutoff=int(filter_group_cutoff)
		filter_panel_list=filter_panel_list.split(',')		
		panel_list+=filter_panel_list
	else:
		filter_cutoff_pval_PT, filter_cutoff_dpsi, filter_cutoff_foc, filter_group_cutoff, filter_panel_list =['','','','',[]]
	return filter_cutoff_pval_PT, filter_cutoff_dpsi, filter_cutoff_foc, filter_group_cutoff, filter_panel_list, panel_list

## screening
def main(args):
	index={}
	fin_list={}
	para_fin=args.parameter_fin
	splicing_event_type=args.splicing_event_type 
	fin_plot_query=args.event_list #IRIS screening output OR any file contain event ID and deltaPSI/indicators
	jc_result=args.jc_full_result #IRIS prev screening output
	tumor_form_col=int(args.deltaPSI_column)-1
	tumor_form_cutoff=float(args.deltaPSI_cut_off)
	step=int(args.step)
	has_header=args.header
	outdir=args.outdir.rstrip('/')

	out_prefix,db_dir,filter1_para,filter2_para,filter3_para=[l.strip() for l in open(para_fin)][:5]
	db_dir=db_dir.rstrip('/')	
	if os.path.isdir(db_dir+'_sjc'): #automatically use db_sjc if in the same dir. Otherwise, use the user input db_dir
		db_dir=db_dir+'_sjc'
	panel_list=[out_prefix]

	filter1_cutoff_pval, filter1_cutoff_dpsi, filter1_cutoff_foc, filter1_group_cutoff, filter1_panel_list, panel_list = loadParametersRow(filter1_para, panel_list)
	filter2_cutoff_pval, filter2_cutoff_dpsi, filter2_cutoff_foc, filter2_group_cutoff, filter2_panel_list, panel_list = loadParametersRow(filter2_para, panel_list)
	filter3_cutoff_pval, filter3_cutoff_dpsi, filter3_cutoff_foc, filter3_group_cutoff, filter3_panel_list, panel_list = loadParametersRow(filter3_para, panel_list)
	if filter1_panel_list==[] and filter2_panel_list==[] and filter3_panel_list==[]:
	# if filter1_panel_list==[] and filter2_panel_list==[] and filter3_panel_list==[] and test_mode[0]!='summary':
		exit("[Error] No filtering required in parameteres file. exit!")

	single_plot=False
	group_plot=True

	#Load JC values into dict based on event file
	JC_dict=loadJCvalue4events(fin_plot_query, tumor_form_col, tumor_form_cutoff, jc_result, splicing_event_type, has_header)
	#Make tmp file for JC plot based events
	query_fin_name=outdir+'/'+fin_plot_query.split('/')[-1]+'.tmpdata.txt'
	fout_tmp=open(query_fin_name,'w')
	for i,l in enumerate(open(fin_plot_query)):
		if i==0 and has_header:
			continue
		ls=l.strip().split('\t')
		event_id=ls[0]
		tumor_form=float(ls[tumor_form_col]) # negative means skipping, pos means inc
		JC_select=selectJC(event_id, tumor_form, tumor_form_cutoff, splicing_event_type)
		if tumor_form<tumor_form_cutoff:
			fout_tmp.write(event_id+'\tskp\t'+JC_dict[JC_select[0]])
		else:
			fout_tmp.write(event_id+'\tinc\t'+JC_dict[JC_select[0]].strip()+'\t'+JC_dict[JC_select[1]])
	fout_tmp.close()
	if group_plot:
		from matplotlib.backends.backend_pdf import PdfPages

	query_fin_len=fileLength(query_fin_name)
	start_point=0
	for r_start in xrange(start_point,query_fin_len,step):
		if group_plot:
			pdf = PdfPages(outdir+'/'+fin_plot_query.split('/')[-1]+'.'+str(r_start)+'.prevalencePlot.pdf')
			fig = plt.figure(figsize=(8.8,12))
			fig.suptitle('', fontsize=14,fontweight='bold')
			sns.set(style="white", color_codes=True)
		for j,k in enumerate(open(query_fin_name)):
			if j < r_start:
				continue
			if j > r_start+step-1:
				break
			ks=k.strip().split('\t')
			as_event_gene=ks[0].split(':')[1]
			if ks[1]=='skp':
				sj_id=ks[2]
				prev={}
				#percentage[out_prefix]=map(float,fetch_PsiMatrix(event_id,fin_list[out_prefix],'.','\t',index[out_prefix])[1][8:])
				prev[out_prefix]=[float(ks[4])]
				for n,group in enumerate(panel_list[1:]):
					prev[group]=[float(ks[6+n*3])]

				prev_df = pd.DataFrame.from_dict(prev, orient='index').transpose()[panel_list]
				if single_plot:
					indiviualPlot(prev_df, event_id, panel_list)
				if group_plot:
					print step, j
					ax_i = plt.subplot2grid((step,11), (j-r_start*1,0), colspan=10, rowspan=1)
					##SHOULD SET TO area. use width just for prelim
					sns_plot=sns.barplot(data=prev_df, ax=ax_i,linewidth=1.5)
					sns_plot.set_yticks(np.arange(0,1.1,0.5))     
					sns_plot.set(xticklabels=[])
					sns.despine(offset=0, trim=False)
					sns_plot.text(15.3, 0.4, as_event_gene+'\n'+sj_id, horizontalalignment='left', size='small', color='black',fontweight='bold')
			else:
				sj_id=ks[2]#TODO
				prev={}
				inc1_list=[float(ks[4])]+[float(ks[6+n*3]) for n in range(len(panel_list[1:]))]
				inc2_list=[float(ks[5+(len(panel_list)-1)*3+2])]+[float(ks[5+(len(panel_list)-1)*3+4+n*3]) for n in range(len(panel_list[1:]))]
				prev_df = pd.DataFrame({'Groups': panel_list,'Inc1': inc1_list , 'Inc2': inc2_list })
				prev_df_tidy = prev_df.melt(id_vars='Groups').rename(columns=str.title)
	
				if single_plot:
					indiviualPlot(prev_df, event_id, panel_list)
				if group_plot:
					print step, j
					ax_i = plt.subplot2grid((step,11), (j-r_start*1,0), colspan=10, rowspan=1)
					##SHOULD SET TO area. use width just for prelim
					sns_plot=sns.barplot(data=prev_df_tidy, x='Groups', y='Value',hue='Variable',ax=ax_i, linewidth=1, ci=None)
					my_pal=sns.husl_palette(len(panel_list)*2, s=.75, l=0.7)#color_palette("husl", len(panel_list))#(n_colors=len(panel_list))
					# for i, bar in enumerate(sns_plot.patches):
					# 	bar.set_color(my_pal[i%len(panel_list)])#(i-i%2)/2
					for i, bar in enumerate(sns_plot.patches):
 						if i<=len(panel_list)-1:
 							bar.set_color(my_pal[i*2])
 						else:
 							bar.set_color(my_pal[1+(i-len(panel_list))*2])
					sns_plot.legend_.remove()
					sns_plot.set_yticks(np.arange(0,1.1,0.5))     
					sns_plot.set(xticklabels=[])
					sns.despine(offset=0, trim=False)
					sns_plot.text(15.3, 0.4, as_event_gene+'\n'+sj_id, horizontalalignment='left', size='small', color='black',fontweight='bold')

		if group_plot:
			pdf.savefig(fig)
			pdf.close()

if __name__ == '__main__':
	main()

