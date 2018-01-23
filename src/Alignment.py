# -*- coding: utf-8 -*-
"""
Created on Monday January 22 10:01:06 2018
@author: clancien
"""


import os
import sys
import subprocess
import logging
from logging.handlers import RotatingFileHandler

#For Python2.7 and Python>3.0
try :
	import ConfigParser
except ImportError:
	import configparser as ConfigParser

class Alignment():

	def __init__(self, file_genome, file_paired_1, file_paired_2):

		
		self.file_genome	=	file_genome
		self.file_paired_1	=	file_paired_1
		self.file_paired_2	=	file_paired_2

		#Error File
		self.logFile		=	config.get('Error', 'logFile')
		self.logger         =	None
		self.formatter      =	None
		self.file_handler   =	None

		self.init_log()

	def init_log(self):
        
		# création de l'objet logger qui va nous servir à écrire dans les logs
		self.logger = logging.getLogger()
		# on met le niveau du logger à DEBUG, comme ça il écrit tout
		self.logger.setLevel(logging.DEBUG)

		# création d'un formateur qui va ajouter le temps, le niveau
		# de chaque message quand on écrira un message dans le log
		self.formatter = logging.Formatter('%(asctime)s :: %(levelname)s :: %(message)s')
		# création d'un handler qui va rediriger une écriture du log vers
		# un fichier en mode 'append', avec 1 backup et une taille max de 1Mo

		self.file_handler = RotatingFileHandler(self.logFile, 'a', 1000000, 1)
		# on lui met le niveau sur DEBUG, on lui dit qu'il doit utiliser le formateur
		# créé précédement et on ajoute ce handler au logger
		self.file_handler.setLevel(logging.DEBUG)
		self.file_handler.setFormatter(self.formatter)
		self.logger.addHandler(self.file_handler)

	def is_file_exists(self):

		""" 
		Check if file exist on system
		"""
		if not os.path.isfile(self.file_genome):
			raise ValueError("File not found : {}".format(str(self.file_genome)))

		if not os.path.isfile(self.file_paired_1):
			raise ValueError("File not found : {}".format(str(self.file_paired_1)))

		if not os.path.isfile(self.file_paired_2):
			raise ValueError("File not found : {}".format(str(self.file_paired_2)))



	def __init__(self,  link_genome_homo_sapiens, link_protein_feature, link_repeat_consensus, link_repeat_feature):

		config=ConfigParser.ConfigParser()
		config.readfp(open('../configuration.ini','r'))

		#path file
		self.path_genome_homo_sapiens=config.get('Genome', 'homo_sapiens')
		self.path_protein_feature=config.get('Feature_Genome','protein_feature')
		self.path_repeat_consensus=config.get('Feature_Genome','repeat_consensus')
		self.path_repeat_feature=config.get('Feature_Genome','repeat_feature')

		#Error File
		self.logFile = config.get('Error', 'logFile')
		self.logger=None
		self.formatter=None
		self.file_handler=None

		self.init_log()
		#Check if directories exist when creating the object
		self.file_exists()

		self.get_genome()
		self.get_protein_feature()
		self.get_repeat_consensus()
		self.get_repeat_feature()

	def init_log(self):
        
		# création de l'objet logger qui va nous servir à écrire dans les logs
		self.logger = logging.getLogger()
		# on met le niveau du logger à DEBUG, comme ça il écrit tout
		self.logger.setLevel(logging.DEBUG)

		# création d'un formateur qui va ajouter le temps, le niveau
		# de chaque message quand on écrira un message dans le log
		self.formatter = logging.Formatter('%(asctime)s :: %(levelname)s :: %(message)s')
		# création d'un handler qui va rediriger une écriture du log vers
		# un fichier en mode 'append', avec 1 backup et une taille max de 1Mo

		self.file_handler = RotatingFileHandler(self.logFile, 'a', 1000000, 1)
		# on lui met le niveau sur DEBUG, on lui dit qu'il doit utiliser le formateur
		# créé précédement et on ajoute ce handler au logger
		self.file_handler.setLevel(logging.DEBUG)
		self.file_handler.setFormatter(self.formatter)
		self.logger.addHandler(self.file_handler)

	def path_exist(self):
		""" Check if dir exists; if not we create the path (+ dir)
			example:
				- string = dir/subdir
				- string.rsplit('/',1)[0]
		"""
		if not os.path.isdir(self.link_genome_homo_sapiens.rsplit('/',1)[0]):
			os.makedirs(self.link_genome_homo_sapiens.rsplit('/',1)[0])

		if not os.path.isdir(self.link_protein_feature.rsplit('/',1)[0]):
			os.makedirs(self.link_protein_feature.rsplit('/',1)[0])
		#if not os.path.isdir

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

	Download(args.link_genome_homo_sapiens, args.link_protein_feature, args.link_repeat_consensus, args.link_repeat_feature)

	print "{}".format("Finished...")






