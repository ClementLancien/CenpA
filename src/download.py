#!/usr/bin/env python

# -*- coding: utf-8 -*-

"""
Created on Monday January 22 10:01:06 2018
@author: clancien
"""

""" 

Add Chip Seq Data

Next updates incoming
Choose Release Number + Organism
Use FTPLib

"""

import os
import sys
import subprocess
import logging
from logging.handlers import RotatingFileHandler
import numpy

import pandas as pd

#For Python2.7 and Python>3.0
try :
	import ConfigParser
except ImportError:
	import configparser as ConfigParser

class Download():

	def __init__(self,  link_genome_homo_sapiens, link_protein_feature, link_repeat_consensus, link_repeat_feature, number_study):

		config=ConfigParser.ConfigParser()
		config.readfp(open('../configuration.ini','r'))

		#link for download
		self.link_genome_homo_sapiens	=	link_genome_homo_sapiens
		self.link_protein_feature 		= 	link_protein_feature
		self.link_repeat_consensus 		=	link_repeat_consensus
		self.link_repeat_feature 		=	link_repeat_feature
		self.number_study 				= 	number_study

		#path to save file
		self.path_genome_homo_sapiens	=	config.get('Genome', 'homo_sapiens')
		self.path_protein_feature		=	config.get('Feature_Genome','protein_feature')
		self.path_repeat_consensus		=	config.get('Feature_Genome','repeat_consensus')
		self.path_repeat_feature		=	config.get('Feature_Genome','repeat_feature')
		self.path_study					= 	os.path.join(str(config.get('ChIP_Seq','study')), self.number_study)


		#Error File
		self.logFile		=	config.get('Error', 'logFile')
		self.logger     	=	None
		self.formatter  	=	None
		self.file_handler	=	None

		self.init_log()
		#Check if directories exist when creating the object
		self.path_exist()


	def init_log(self):
        
		# creation de l'objet logger qui va nous servir a ecrire dans les logs
		self.logger = logging.getLogger()
		# on met le niveau du logger a DEBUG, comme ca il ecrit tout
		self.logger.setLevel(logging.DEBUG)

		# creation d'un formateur qui va ajouter le temps, le niveau
		# de chaque message quand on ecrira un message dans le log
		self.formatter = logging.Formatter('%(asctime)s :: %(levelname)s :: %(message)s')
		# creation d'un handler qui va rediriger une ecriture du log vers
		# un fichier en mode 'append', avec 1 backup et une taille max de 1Mo

		self.file_handler = RotatingFileHandler(self.logFile, 'a', 1000000, 1)
		# on lui met le niveau sur DEBUG, on lui dit qu'il doit utiliser le formateur
		# cree precedement et on ajoute ce handler au logger
		self.file_handler.setLevel(logging.DEBUG)
		self.file_handler.setFormatter(self.formatter)
		self.logger.addHandler(self.file_handler)

	def path_exist(self):
		""" Check if dir exists; if not we create the path (+ dir)
			example:
				- string = dir/subdir
				- string.rsplit('/',1)[0]
		"""
		if not os.path.isdir(self.path_genome_homo_sapiens.rsplit('/',1)[0]):
			# print "1"
			# print self.path_genome_homo_sapiens.rsplit('/',1)[0]
			os.makedirs(self.path_genome_homo_sapiens.rsplit('/',1)[0])

		if not os.path.isdir(self.path_protein_feature.rsplit('/',1)[0]):
			# print "2"
			# print os.path.isdir(self.link_protein_feature.rsplit('/',1)[0])
			# print self.link_protein_feature
			# print self.path_protein_feature.rsplit('/',1)[0]
			os.makedirs(self.path_protein_feature.rsplit('/',1)[0])

		if not os.path.isdir(self.path_study,):
			os.makedirs(self.path_study)

	def get_genome(self):

		try:

			subprocess.check_output(['bash', '-c', "wget -nc " + self.link_genome_homo_sapiens + " -P " + self.path_genome_homo_sapiens.rsplit('/',1)[0] + " &>/dev/null"])

		except subprocess.CalledProcessError as error:

			self.logger.warning("Error - download.py - get_genome " + str(self.path_genome_homo_sapiens.rsplit('/',1))[1])
			self.logger.warning("Exception at the line : {}".format(sys.exc_info()[-1].tb_lineno))
			self.logger.warning(sys.exc_info())

	def get_protein_feature(self):

		try:

			subprocess.check_output(['bash', '-c', "wget -nc " + self.link_protein_feature + " -P " + self.path_protein_feature.rsplit('/',1)[0] + " &>/dev/null"])

		except subprocess.CalledProcessError as error:

			self.logger.warning("Error - download.py - get_protein_feature " +str(self.path_protein_feature.rsplit('/',1))[1])
			self.logger.warning("Exception at the line : {}".format(sys.exc_info()[-1].tb_lineno))
			self.logger.warning(sys.exc_info())

	def get_repeat_consensus(self):

		try:

			subprocess.check_output(['bash', '-c', "wget -nc " + self.link_repeat_consensus + " -P " + self.path_repeat_consensus.rsplit('/',1)[0] + " &>/dev/null"])

		except subprocess.CalledProcessError as error:

			self.logger.warning("Error - download.py - get_repeat_consensus " +str(self.path_repeat_consensus.rsplit('/',1))[1])
			self.logger.warning("Exception at the line : {}".format(sys.exc_info()[-1].tb_lineno))
			self.logger.warning(sys.exc_info())

	def get_repeat_feature(self):

		try:

			subprocess.check_output(['bash', '-c', "wget -nc " + self.link_repeat_feature + " -P " + self.path_repeat_feature.rsplit('/',1)[0] + " &>/dev/null"])

		except subprocess.CalledProcessError as error:

			self.logger.warning("Error - download.py - get_repeat_feature " +str(self.path_repeat_feature.rsplit('/',1))[1])
			self.logger.warning("Exception at the line : {}".format(sys.exc_info()[-1].tb_lineno))
			self.logger.warning(sys.exc_info())

	def get_study_info(self):
		# print "'"+"https://www.ebi.ac.uk/ena/data/warehouse/filereport?accession=" +self.number_study+ "result=read_run&fields=study_accession,sample_accession,secondary_sample_accession,experiment_accession,run_accession,tax_id,scientific_name,instrument_model,library_layout,fastq_ftp,fastq_galaxy,submitted_ftp,submitted_galaxy,sra_ftp,sra_galaxy,cram_index_ftp,cram_index_galaxy&download=txt"+"'"
		# exit(0)
		# print  "wget -nc "+"'"+"https://www.ebi.ac.uk/ena/data/warehouse/filereport?accession=" +self.number_study+ "result=read_run&fields=study_accession,sample_accession,secondary_sample_accession,experiment_accession,run_accession,tax_id,scientific_name,instrument_model,library_layout,fastq_ftp,fastq_galaxy,submitted_ftp,submitted_galaxy,sra_ftp,sra_galaxy,cram_index_ftp,cram_index_galaxy&download=txt"+"'" + " --output-document=" + os.path.join(self.path_study, self.number_study) + ".txt &>/dev/null"
		# exit(0)
		try:

			subprocess.check_output(['bash', '-c', "wget -nc "+"'"+"https://www.ebi.ac.uk/ena/data/warehouse/filereport?accession=" +self.number_study+ "&result=read_run&fields=study_accession,secondary_study_accession,sample_accession,secondary_sample_accession,experiment_accession,run_accession,submission_accession,tax_id,scientific_name,instrument_platform,instrument_model,library_name,nominal_length,library_layout,library_strategy,library_source,library_selection,read_count,base_count,center_name,first_public,last_updated,experiment_title,study_title,study_alias,experiment_alias,run_alias,fastq_bytes,fastq_md5,fastq_ftp,fastq_aspera,fastq_galaxy,submitted_bytes,submitted_md5,submitted_ftp,submitted_aspera,submitted_galaxy,submitted_format,sra_bytes,sra_md5,sra_ftp,sra_aspera,sra_galaxy,cram_index_ftp,cram_index_aspera,cram_index_galaxy,sample_alias,broker_name,sample_title,nominal_sdev,first_created&download=txt"+"'" + " --output-document=" + os.path.join(self.path_study, self.number_study) + ".txt &>/dev/null"])

		except subprocess.CalledProcessError as error:

			self.logger.warning("Error - download.py - get_study " +str(self.path_repeat_feature.rsplit('/',1))[1])
			self.logger.warning("Exception at the line : {}".format(sys.exc_info()[-1].tb_lineno))
			self.logger.warning(sys.exc_info())

	def get_study(self):
		
		filename = str(os.path.join(self.path_study, self.number_study) + ".txt")
		df = pd.read_csv(filename, header=0, sep='\t', dtype=str)
		print df

if __name__ == '__main__':
			
	import optparse
	from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter, SUPPRESS

	def get_parser():
		"""
			Get Parser Object for this script
		"""

		parser = ArgumentParser(description=__doc__,
								formatter_class=ArgumentDefaultsHelpFormatter,
								add_help=False)

		parser.add_argument("-g", "--genome",
							dest="link_genome_homo_sapiens",
							default="ftp://ftp.ensembl.org/pub/release-91/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna_sm.toplevel.fa.gz",
							type=str,
							help="Link to download your genome",
							metavar="File")

		parser.add_argument("-p", "--protein",
							dest="link_protein_feature",
							default="ftp://ftp.ensembl.org/pub/release-91/mysql/homo_sapiens_otherfeatures_91_38/protein_feature.txt.gz",
							type=str,
							help="Link to download your protein feature",
							metavar="File")

		parser.add_argument("-c", "--consensus",
							dest="link_repeat_consensus",
							default="ftp://ftp.ensembl.org/pub/release-91/mysql/homo_sapiens_otherfeatures_91_38/repeat_consensus.txt.gz",
							type=str,
							help="Link to download your repeat consensus",
							metavar="File")

		parser.add_argument("-f", "--feature",
							dest="link_repeat_feature",
							default="ftp://ftp.ensembl.org/pub/release-91/mysql/homo_sapiens_otherfeatures_91_38/repeat_feature.txt.gz",
							type=str,
							help="Link to download your repeat feature",
							metavar="File")

		parser.add_argument("-s", "--study",
							dest="link_study",
							default="PRJNA183936",
							type=str,
							help="Link to download information about the study from EMBL-EBI",
							metavar="File"
							)

		parser.add_argument('-v', '--version', 
							action='version',
							version='%(prog)s 1.0',
							help="Version 1.0")

		parser.add_argument('-h', '--help',
							action='help', 
							default=SUPPRESS,
							help='Tool to download Genome and features from Ensembl')

		return parser

	parser=get_parser()
	args=parser.parse_args()

	print "{}".format("Starting download...")
	print "{}".format("Please wait a moment...")

	download = Download(args.link_genome_homo_sapiens, args.link_protein_feature, args.link_repeat_consensus, args.link_repeat_feature, args.link_study)
	
	download.get_genome()
	download.get_protein_feature()
	download.get_repeat_consensus()
	download.get_repeat_feature()
	download.get_study_info()
	download.get_study()

	print "{}".format("Finished...")






