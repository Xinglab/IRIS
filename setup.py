#!/usr/bin/env python

from setuptools import setup


def main():
	setup(
		name='IRIS',
		  
		version='1.0.0',
		  
		description='Isoform peptides from RNA splicing for Immunotherapy target Screening',

		author='Yang Pan',

		author_email='panyang@ucla.edu',

		url='',

		packages=['IRIS','IRIS.data'],

		scripts=['bin/IRIS'],

		include_package_data=True,

		package_data={'IRIS.data':[
		'brain_blacklistMay.txt',
		'features.uniprot2gtf.ExtraCell.txt',
		'UniprotENSGmap.txt',
		'uniprot2gtf.blastout.uniprotAll.txt',
		'HLA_types.least.list',
		'HLA_types.least.tsv']},
		install_requires=[]#'keras', 
				#'numpy',
				#'pyyaml',
				#'h5py',
				#'scikit-learn',
				#'scipy',
				#'tqdm>=4.14',
				#'pandas>=0.21.0',
				#'theano']
		 )
	return

if __name__ == '__main__':
	main()
