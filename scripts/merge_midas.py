#!/usr/bin/env python

# MIDAS: Metagenomic Intra-species Diversity Analysis System
# Copyright (C) 2015 Stephen Nayfach
# Freely distributed under the GNU General Public License (GPLv3)

import os, sys, platform, argparse, numpy as np
from midas import utility

def get_program():
	""" Get program specified by user (species, genes, or snps) """
	if len(sys.argv) == 1 or sys.argv[1] in ['-h', '--help']:
		print('')
		print('Usage: merge_midas.py <command> [options]')
		print('')
		print('Note: use merge_midas.py <command> -h to view usage for a specific command')
		print('')
		print('Commands:')
		print('\tspecies\t merge abundances of bacterial species across samples')
		print('\tgenes\t merge pan-genome gene copy numbers of species across samples')
		print('\tsnps\t merge single nucleotide variants of species across samples')
		quit()
	elif sys.argv[1] not in ['species', 'genes', 'snps']:
		sys.exit("\nError: Unrecognized command: '%s'\n" % sys.argv[1])
		quit()
	else:
		return sys.argv[1]

def get_arguments(program):
	""" Get arguments for specified program """
	if program == 'species':
		args = species_arguments()
	elif program == 'genes':
		args = genes_arguments()
	elif program == 'snps':
		args = snps_arguments()
	else:
		sys.exit("\nError: Unrecognized program: '%s'\n" % program)
	return args

def species_arguments():
	parser = argparse.ArgumentParser(
		formatter_class=argparse.RawTextHelpFormatter,
		usage=argparse.SUPPRESS,
		description="""
Description: Merge species abundance files across samples

Usage: merge_midas.py species outdir [options]
""",
		epilog="""Examples:
1) provide list of paths to sample directories:
merge_midas.py species /path/to/outdir -i /path/to/samples/sample_1,/path/to/samples/sample_2 -t list

2) provide directory containing all samples:
merge_midas.py species /path/to/outdir -i /path/to/samples -t dir

3) provide file containing paths to sample directories:
merge_midas.py species /path/to/outdir -i /path/to/samples/sample_paths.txt -t file

4) run a quick test:
merge_midas.py species /path/to/outdir -i /path/to/samples -t dir --max_samples 2
""")
	parser.add_argument('program', help=argparse.SUPPRESS)
	parser.add_argument('outdir', type=str, help='Directory for output files')
	parser.add_argument('-i', type=str, dest='input', required=True,
		help="""Input to sample directories output by run_midas.py
Can be a list of directories, a directory containing all samples, or a file with paths
See '-t' for details""")
	parser.add_argument('-t', choices=['list','file','dir'], dest='intype', required=True,
		help="""list: -i is a comma-separated list (ex: /path/to/samples/sample_1,/path/to/samples/sample_2)
 dir: -i is a directory containing all samples (ex: /path/to/samples)
file: -i is a file containing paths to sample directories (ex: /path/to/sample_paths.txt)
""")
	parser.add_argument('-d', type=str, dest='db', default=os.environ['MIDAS_DB'] if 'MIDAS_DB' in os.environ else None,
		help="""Path to reference database
By default the MIDAS_DB environmental variable is used""")
	parser.add_argument('--min_cov', metavar='FLOAT', type=float, default=1.0,
		help="""Minimum marker-gene-coverage for estimating species prevalence (1.0)""")
	parser.add_argument('--max_samples', type=int, metavar='INT',
		help="""Maximum number of samples to process.
Useful for testing (use all)""")
	args = vars(parser.parse_args())
	return args

def genes_arguments():
	""" Parse command line arguments """
	parser = argparse.ArgumentParser(
		formatter_class=argparse.RawTextHelpFormatter,
		usage=argparse.SUPPRESS,
		description="""
Description: merge results from pan-genome profiling across samples

Usage: merge_midas.py genes outdir [options]
""",
		epilog="""Examples:
1) Merge results for all species. Provide list of paths to sample directories:
merge_midas.py genes /path/to/outdir -i sample_1,sample_2 -t list

2) Merge results for one species (id=Bacteroides_vulgatus_57955):
merge_midas.py genes /path/to/outdir --species_id Bacteroides_vulgatus_57955 -i sample_1,sample_2 -t list

3) Exclude low-coverage samples in output matrix:
merge_midas.py genes /path/to/outdir -i /path/to/samples -t dir --sample_depth 5.0

4) Use lenient threshold for determining gene presence-absence:
merge_midas.py genes /path/to/outdir -i /path/to/samples -t dir --min_copy 0.1

5) Run a quick test:
merge_midas.py genes /path/to/outdir -i /path/to/samples -t dir --max_species 1 --max_samples 10

""")
	parser.add_argument('program', help=argparse.SUPPRESS)
	parser.add_argument('outdir', type=str,
		help="directory for output files. a subdirectory will be created for each species_id")
	io = parser.add_argument_group('Input/Output')
	io.add_argument('-i', type=str, dest='input', required=True,
		help="""Input to sample directories output by run_midas.py
Can be a list of directories, a directory containing all samples, or a file with paths
See '-t' for details""")
	io.add_argument('-t', choices=['list','file','dir'], dest='intype', required=True,
		help="""list: -i is a comma-separated list (ex: /path/to/samples/sample_1,/path/to/samples/sample_2)
 dir: -i is a directory containing all samples (ex: /path/to/samples)
file: -i is a file containing paths to sample directories (ex: /path/to/sample_paths.txt)
""")
	io.add_argument('-d', type=str, dest='db', default=os.environ['MIDAS_DB'] if 'MIDAS_DB' in os.environ else None,
		help="""Path to reference database
By default, the MIDAS_DB environmental variable is used""")
	species = parser.add_argument_group('Species filters (select subset of species from INPUT)')
	species.add_argument('--min_samples', type=int, default=1, metavar='INT',
		help="""All species with >= MIN_SAMPLES (1)""")
	species.add_argument('--species_id', dest='species_id', type=str, metavar='CHAR',
		help="""Comma-separated list of species ids""")
	species.add_argument('--max_species', type=int, metavar='INT',
		help="""Maximum number of species to merge. Useful for testing (use all)""")
	sample = parser.add_argument_group('Sample filters (select subset of samples from INPUT)')
	sample.add_argument('--sample_depth', type=float, default=1.0, metavar='FLOAT',
		help="""Minimum read-depth across all genes with non-zero coverage (1.0)""")
	sample.add_argument('--max_samples', type=int, metavar='INT',
		help="""Maximum number of samples to process. Useful for testing (use all)""")
	gene = parser.add_argument_group('Quantification')
	gene.add_argument('--cluster_pid', type=str, dest='cluster_pid', default='95', choices=['75', '80', '85', '90', '95', '99'],
		help="""In the database, pan-genomes are defined at 6 different %% identity clustering cutoffs
CLUSTER_PID allows you to quantify gene content for any of these sets of gene clusters
By default, gene content is reported for genes clustered at 95%% identity (95)
""")
	gene.add_argument('--min_copy', type=float, default=0.35, metavar='FLOAT',
		help="""Genes >= MIN_COPY are classified as present
Genes < MIN_COPY are classified as absent (0.35)""")
	args = vars(parser.parse_args())
	return args

def snps_arguments():
	""" Parse command line arguments """
	parser = argparse.ArgumentParser(
		formatter_class=argparse.RawTextHelpFormatter,
		usage=argparse.SUPPRESS,
		description="""
Description: merge single-nucleotide variant results across samples

Usage: merge_midas.py snps outdir [options]
""",
		epilog="""Examples:
1) Merge results for all species. Provide list of paths to sample directories:
merge_midas.py snps /path/to/outdir -i sample_1,sample_2 -t list

2) Merge results for one species (id=Bacteroides_vulgatus_57955):
merge_midas.py snps /path/to/outdir --species_id Bacteroides_vulgatus_57955 -i sample_1,sample_2 -t list

3) Only use samples with >15x average depth and only use sites covered by >=10 reads in at least >=95% of samples:
merge_midas.py snps /path/to/outdir -i /path/to/samples -t dir --sample_depth 15 --site_depth 10 --site_prev 0.95

4) Run a quick test:
merge_midas.py snps /path/to/outdir -i /path/to/samples -t dir --max_species 1 --max_samples 10 --max_sites 1000

""")
	parser.add_argument('program', help=argparse.SUPPRESS)
	parser.add_argument('outdir', type=str,
		help="Directory for output files. a subdirectory will be created for each species_id")
	parser.add_argument('--threads', type=int, default=1, metavar='INT',
		help="Number of CPUs to use for merging files (1)\nIncreases speed when merging many species")
	parser.add_argument('--sites_per_iter', type=int, default=20000, metavar='INT',
		help=argparse.SUPPRESS)
	io = parser.add_argument_group('Input/Output')
	io.add_argument('-i', type=str, dest='input', required=True,
		help="""Input to sample directories output by run_midas.py
Can be a list of directories, a directory containing all samples, or a file with paths
See '-t' for details""")
	io.add_argument('-t', choices=['list','file','dir'], dest='intype', required=True,
		help="""list: -i is a comma-separated list (ex: /path/to/samples/sample_1,/path/to/samples/sample_2)
 dir: -i is a directory containing all samples (ex: /path/to/samples)
file: -i is a file containing paths to sample directories (ex: /path/to/sample_paths.txt)
""")
	io.add_argument('-d', type=str, dest='db', default=os.environ['MIDAS_DB'] if 'MIDAS_DB' in os.environ else None,
		help="""Path to reference database
By default, the MIDAS_DB environmental variable is used""")
	presets = parser.add_argument_group("Presets (option groups for easily...)")
	
	species = parser.add_argument_group("Species filters (select subset of species from INPUT)")
	species.add_argument('--min_samples', type=int, default=1, metavar='INT',
		help="""All species with >= MIN_SAMPLES (1)""")
	species.add_argument('--species_id', dest='species_id', type=str, metavar='CHAR',
		help="""Comma-separated list of species ids""")
	species.add_argument('--max_species', type=int, metavar='INT',
		help="""Maximum number of species to merge (use all)""")
	sample = parser.add_argument_group("Sample filters (select subset of samples from INPUT)")
	sample.add_argument('--sample_depth', dest='sample_depth', type=float, default=5.0, metavar='FLOAT',
		help="""Minimum average read depth per sample (5.0)""")
	sample.add_argument('--fract_cov', dest='fract_cov', type=float, default=0.4, metavar='FLOAT',
		help="""Fraction of reference sites covered by at least 1 read (0.4)""")
	sample.add_argument('--max_samples', type=int, metavar='INT',
		help="""Maximum number of samples to process. useful for quick tests (use all)""")
	snps = parser.add_argument_group("Site filters (select subset of genomic sites from INPUT)")
	snps.add_argument('--site_depth', type=int, default=3, metavar='INT',
		help="""Minimum number of mapped reads per site (3)
A high value like 20 will result in accurate allele frequencies, but may discard many sites.
A low value like 1 will retain many sites but may not result in accurate allele frequencies""")
	snps.add_argument('--site_ratio', type=float, default=3.0, metavar='FLOAT',
		help="""Maximum ratio of site depth to genome depth (3.0)
This filter helps to eliminate genomic sites with abnormally high read depth""")
	snps.add_argument('--site_prev', type=float, default=0.95, metavar='FLOAT',
		help="""Genomic site is present in at least <site_prev> proportion of samples (0.95)
Site presence/absence is determined by --site_depth and --site_ratio.
A value of 1.0 will select sites that are present in all samples.
A value of 0.0 will select all sites, including those with low coverage""")
	snps.add_argument('--site_maf', type=float, default=0.0, metavar='FLOAT',
		help="""Minimum pooled minor allele frequency of site (0.0)
This filter helps to eliminate invariant genomic sites.
Values above zero (e.g. 0.01, 0.02, 0.05) will only keep common variants""")
	snps.add_argument('--site_multi_freq', type=float, default=0.01, metavar='FLOAT',
		help="""Maximum pooled frequency of 3rd and 4th alleles (0.01)
This filter helps to eliminate sites with more than 2 common variants""")
	snps.add_argument('--max_sites', type=int, default=float('Inf'), metavar='INT',
		help="""Maximum number of sites to include in output (use all)
Useful for quick tests """)
	args = vars(parser.parse_args())
	return args

def check_arguments(program, args):
	""" Run program specified by user (species, genes, or snps) """
	if program in ['species', 'snps', 'genes']:
		if not os.path.isdir(args['outdir']): os.mkdir(args['outdir'])
		check_input(args)
		utility.check_database(args)
	else:
		sys.exit("\nError: Unrecognized program: '%s'\n" % program)
	if platform.system() not in ['Linux', 'Darwin']:
		sys.exit("\nError: Operating system '%s' not supported\n" % system())

def check_input(args):
	args['indirs'] = []
	error = "\nError: specified input %s does not exist: %s\n"
	if args['intype'] == 'dir':
		if not os.path.isdir(args['input']):
			sys.exit(error % (args['intype'], os.path.abspath(args['input'])))
		else:
			for dir in os.listdir(args['input']):
				args['indirs'].append(os.path.join(args['input'], dir))
	elif args['intype'] == 'file':
		if not os.path.isfile(args['input']):
			sys.exit(error % (args['intype'], os.path.abspath(args['input'])))
		else:
			for line in open(args['input']):
				dir = line.rstrip().rstrip('/')
				if not os.path.isdir(dir): sys.exit(error % ('dir', dir))
				else: args['indirs'].append(dir)
	elif args['intype'] == 'list':
		for dir in args['input'].split(','):
			if not os.path.isdir(dir): sys.exit(error % ('dir', dir))
			else: args['indirs'].append(dir)

def print_arguments(program, args):
	""" Run program specified by user (species, genes, or snps) """
	if program == 'species':
		print_species_arguments(args)
	elif program == 'genes':
		print_genes_arguments(args)
	elif program == 'snps':
		print_snps_arguments(args)
	else:
		sys.exit("\nError: Unrecognized program: '%s'\n" % program)

def print_species_arguments(args):
	print ("===========Parameters===========")
	print ("Command: %s" % ' '.join(sys.argv))
	print ("Script: merge_midas.py species")
	print ("Database: %s" % args['db'])
	print ("Input: %s" % args['input'])
	print ("Input type: %s" % args['intype'])
	print ("Output directory: %s" % args['outdir'])
	print ("Minimum coverage for estimating prevalence: %s" % args['min_cov'])
	if args['max_samples']: print ("Keep <= %s samples" % args['max_samples'])
	print ("===============================")
	print ("")

def print_genes_arguments(args):
	print ("===========Parameters===========")
	print ("Command: %s" % ' '.join(sys.argv))
	print ("Script: merge_midas.py genes")
	print ("Database: %s" % args['db'])
	print ("Input: %s" % args['input'])
	print ("Input type: %s" % args['intype'])
	print ("Output directory: %s" % args['outdir'])
	print ("Species selection criteria:")
	if args['species_id']: print ("  keep species ids: %s" % args['species_id'].split(','))
	else: print ("  keep species with >= %s samples" % args['min_samples'])
	if args['max_species']: print ("  keep <= %s species" % args['max_species'])
	print ("Sample selection criteria:")
	print ("  keep samples with >=%s mean coverage across genes with non-zero coverage" % args['sample_depth'])
	if args['max_samples']: print ("  keep <= %s samples" % args['max_samples'])
	print ("Gene quantification criterea:")
	print ("  quantify genes clustered at %s%% identity" % args['cluster_pid'])
	print ("  present (1): genes with copy number >= %s" % args['min_copy'])
	print ("  absent (0): genes with copy number < %s" % args['min_copy'])
	print ("===============================")
	print ("")

def print_snps_arguments(args):
	print ("===========Parameters===========")
	print ("Command: %s" % ' '.join(sys.argv))
	print ("Script: merge_midas.py snps")
	print ("Database: %s" % args['db'])
	print ("Input: %s" % args['input'])
	print ("Input type: %s" % args['intype'])
	print ("Output directory: %s" % args['outdir'])
	print ("Species selection criteria:")
	if args['species_id']: print ("  keep species ids: %s" % args['species_id'].split(','))
	else: print ("  keep species with >= %s samples" % args['min_samples'])
	if args['max_species']: print ("  keep <= %s species" % args['max_species'])
	print ("Sample selection criteria:")
	print ("  keep samples with >= %s mean coverage across sites with non-zero coverage" % args['sample_depth'])
	print ("  keep samples where >= %s percent of sites have non-zero coverage" % (100*args['fract_cov']))
	if args['max_samples']: print ("  keep <= %s samples" % args['max_samples'])
	print ("Site selection criteria:")
	print ("  keep sites covered by >= %s reads across >= %s percent of samples" % (args['site_depth'], 100*args['site_prev']))
	if args['max_sites'] != float('Inf'): print ("  keep <= %s sites" % (args['max_sites']))
	print ("Number of CPUs to use: %s" % args['threads'])
	print ("===============================")
	print ("")

def run_program(program, args):
	""" Run program specified by user (species, genes, or snps) """
	if program == 'species':
		from midas.merge import merge_species
		merge_species.run_pipeline(args)
	elif program == 'genes':
		from midas.merge import merge_genes
		merge_genes.run_pipeline(args)
	elif program == 'snps':
		from midas.merge import merge_snps
		merge_snps.run_pipeline(args)
	else:
		sys.exit("\nError: Unrecognized program: '%s'\n" % program)

if __name__ == '__main__':
	program = get_program()
	args = get_arguments(program)
	check_arguments(program, args)
	utility.print_copyright()
	print_arguments(program, args)
	run_program(program, args)



