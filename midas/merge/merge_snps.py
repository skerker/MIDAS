#!/usr/bin/env python

# MIDAS: Metagenomic Intra-species Diversity Analysis System
# Copyright (C) 2015 Stephen Nayfach
# Freely distributed under the GNU General Public License (GPLv3)

import sys, os, numpy as np, gzip, csv, random, Bio.SeqIO
from time import time
from operator import itemgetter
from midas import utility
from midas.merge import merge
import multiprocessing as mp

def format_dict(d):
	""" Format dictionary. ex: 'A:SYN|C:NS|T:NS|G:NS' """
	return '|'.join(['%s:%s' % (x, y) for x, y in d.items()])

def format_list(x, y):
	""" Format list. ex: 'A:SYN|C:NS|T:NS|G:NS' """
	return '|'.join(['%s:%s' % (i, j) for i, j in zip(x,y)])

def write_readme(species, outdir, dbdir):
	outfile = open('%s/%s/readme.txt' % (outdir, species.id), 'w')
	outfile.write("""
Description of output files and file formats from 'merge_midas.py snps'

Output files
############
snps_freq.txt
  frequency of minor allele per genomic site and per sample
snps_depth.txt
  number of reads mapped to genomic site per sample
snps_info.txt  
  metadata for genomic site
snps_summary.txt
  alignment summary statistics per sample
snps_log.txt
  log file containing parameters used

Output formats
############
snps_freq.txt and snps_depth.txt
  tab-delimited matrix files
  field names are sample ids
  row names are genome site ids
snps_summary.txt
  sample_id: sample identifier
  genome_length: number of base pairs in representative genome
  covered_bases: number of reference sites with at least 1 mapped read
  fraction_covered: proportion of reference sites with at least 1 mapped read
  mean_coverage: average read-depth across reference sites with at least 1 mapped read
snps_info.txt
  site_id: genomic site_id; format: ref_id|ref_pos|ref_allele
  mean_freq: average frequency of reference allele across samples
  mean_depth: average read-depth across samples
  site_prev: proportion of samples where site_id was covered with sufficient depth
  allele_props: pooled frequency of 4 nucleotides
  site_type: NC (non-coding), 1D, 2D, 3D, 4D (degeneracy)
  gene_id: gene that intersects site
  amino_acids: protein affect of all 4 possible nucleotides
  snps: SYN/NS for all 4 possible nucleotides

Additional information for species can be found in the reference database:
 %s/rep_genomes/%s
""" % (dbdir, species.id) )
	outfile.close()

def read_genes(db, species_id, contigs):
	""" Read in gene coordinates from features file """
	genes_path = '%s/rep_genomes/%s/genome.features.gz' % (db, species_id)
	genes = []
	for gene in utility.parse_file(genes_path):
		if gene['gene_type'] == 'RNA':
			continue
		else:
			gene['start'] = int(gene['start'])
			gene['end'] = int(gene['end'])
			gene['seq'] = get_gene_seq(gene, contigs[gene['scaffold_id']])
			genes.append(gene)
	return genes

def read_genome(db, species_id):
	""" Read in representative genome from reference database """
	inpath = '%s/rep_genomes/%s/genome.fna.gz' % (db, species_id)
	infile = utility.iopen(inpath)
	genome = {}
	for r in Bio.SeqIO.parse(infile, 'fasta'):
		genome[r.id] = r.seq.upper()
	infile.close()
	return genome

def get_gene_seq(gene, contig):
	""" Fetch nucleotide sequence of gene from genome """
	seq = contig[gene['start']-1:gene['end']] # 2x check this works for + and - genes
	if gene['strand'] == '-':
		return(rev_comp(seq))
	else:
		return(seq)

def complement(base):
	""" Complement nucleotide """
	d = {'A':'T', 'T':'A', 'G':'C', 'C':'G'}
	if base in d: return d[base]
	else: return base

def rev_comp(seq):
	""" Reverse complement sequence """
	return(''.join([complement(base) for base in list(seq[::-1])]))

def translate(codon):
	""" Translate individual codon """
	codontable = {
	'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
	'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
	'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
	'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
	'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
	'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
	'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
	'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
	'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
	'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
	'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
	'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
	'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
	'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
	'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
	'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W',
	}
	return codontable[str(codon)]

def index_replace(codon, allele, pos, strand):
	""" Replace character at index i in string x with y"""
	bases = list(codon)
	bases[pos] = allele if strand == '+' else complement(allele)
	return(''.join(bases))

class GenomicSite:
	def __init__(self, nsamples):
		self.id = None
		self.alleles = {'A':0, 'T':1, 'C':2, 'G':3} # <dic>, mapping of allele to list index
		
		# site annotations
		self.ref_id = None
		self.ref_pos = None
		self.site_type = None
		self.gene_id = None
		self.amino_acids = None
		self.ref_codon = None
		self.codon_pos = None
		
		# site info
		self.ref_allele = None
		self.major_allele = None
		self.minor_count = None
		self.minor_freq = None
		self.prevalence = None
		
		# pooled statistics
		self.pooled_counts = [0,0,0,0] # <list>, count of each allele across all samples
		self.pooled_depth = None # <int>, count of all alleles from all samples
		self.pooled_freq = None # <list>, allele frequency from pooled samples
		
		# per-sample statistics
		self.sample_counts = [] # <list>, [A,T,C,G] counts per sample
		self.sample_freqs = [] # <list>, minor allele frequencies
		self.sample_depths = [] # <list>, count of major+minor allele frequencies
				
	def compute_pooled_counts(self):
		""" Compute pooled nucleotide statistics (counts, depth, freq) at GenomicSite """
		for counts in self.sample_counts:
			for i in range(4):
				self.pooled_counts[i] += counts[i]
		self.pooled_depth = sum(self.pooled_counts)
		if self.pooled_depth > 0:
			self.pooled_freq = [self.pooled_counts[i]/float(self.pooled_depth) for i in range(4)]
		else:
			self.pooled_freq = [0.0]*4

	def call_major_minor_alleles(self):
		""" Call major and minor alleles at GenomicSite """
		if self.pooled_depth > 0:
			allele_counts = [(i,j) for i,j in zip(list('ATCG'), self.pooled_counts)]
			sorted_counts = sorted(allele_counts, key=itemgetter(1), reverse=True)
			self.major_allele, self.major_count = sorted_counts[0]
			self.minor_allele, self.minor_count = sorted_counts[1]
			self.minor_freq = float(self.minor_count)/(self.minor_count+self.major_count)
		else:
			self.major_allele, self.minor_allele = random.sample(list('ATCG'), 2)
			self.minor_freq = 0.0
	
	def compute_multi_freq(self):
		""" Compute combined frequency of 3rd and 4rd alleles """
		self.multi_freq = 0
		for allele in list('ATCG'):
			if not allele in [self.major_allele, self.minor_allele]:
				self.multi_freq += self.pooled_freq[self.alleles[allele]]

	def compute_minor_freqs(self):
		""" Compute per-sample depth (major + minor allele) and minor allele freq at GenomicSite """
		for counts in self.sample_counts:
			count_major = counts[self.alleles[self.major_allele]]
			count_minor = counts[self.alleles[self.minor_allele]]
			sample_depth = count_major + count_minor
			sample_freq = float(count_minor)/sample_depth if sample_depth > 0 else 0.0
			self.sample_freqs.append(sample_freq)
			self.sample_depths.append(sample_depth)

	def compute_prevalence(self, species, min_depth, max_ratio):
		""" Compute the fraction of samples where site passes all filters """
		pass_qc = []
		for mean_depth, site_depth in zip(species.sample_depth, self.sample_depths):
			if site_depth < min_depth:
				pass_qc.append(0)
			elif site_depth/mean_depth > max_ratio:
				pass_qc.append(0)
			else:
				pass_qc.append(1)
		self.prevalence = sum(pass_qc)/float(len(pass_qc))

	def filter(self, min_maf=0.0, min_prev=0.0, max_multi_freq=1.0):
		""" Filter genomic site based on MAF and prevalence """
		if self.minor_freq < min_maf:
			self.remove = True
			self.flag = 'min_maf'
		elif self.prevalence < min_prev:
			self.remove = True
			self.flag = 'min_prev'
		elif self.multi_freq > max_multi_freq:
			self.remove = True
			self.flag = 'max_multi_freq'
		else:
			self.remove = False
			self.flag = None

	def annotate(self, genes, gene_index, contigs):
		""" Annotate variant and reference site """
		# genes: list of genes, each gene contains info
		# contig: contig sequence
		# gene_index: current position in list of genes; global variable
		self.amino_acids = {}
		while True:
			# 1. fetch next gene
			#    if there are no more genes, snp must be non-coding so break
			if gene_index[0] < len(genes):
				gene = genes[gene_index[0]]
			else:
				self.site_type = 'NC'; self.gene_id = ''
				return
			# 2. if snp is upstream of next gene, snp must be non-coding so break
			if (self.ref_id < gene['scaffold_id'] or
			   (self.ref_id == gene['scaffold_id'] and self.ref_pos < gene['start'])):
				self.site_type = 'NC'; self.gene_id = ''
				return
			# 3. if snp is downstream of next gene, pop gene, check (1) and (2) again
			if (self.ref_id > gene['scaffold_id'] or
			   (self.ref_id == gene['scaffold_id'] and self.ref_pos > gene['end'])):
				gene_index[0] += 1
				continue
			# 4. otherwise, snp must be in gene
			#    annotate site (1D-4D)
			else:
				self.gene_id = gene['gene_id']
				self.ref_codon, self.codon_pos = self.fetch_ref_codon(gene)
				if not all([_ in ['A','T','C','G'] for _ in self.ref_codon]): # check for invalid bases in codon
					self.site_type = 'NA'; self.gene_id = ''
				else:
					for allele in ['A','T','C','G']: # + strand
						codon = index_replace(self.ref_codon, allele, self.codon_pos, gene['strand']) # +/- strand
						self.amino_acids[allele] = translate(codon)
					unique_aa = set(self.amino_acids.values())
					degeneracy = 4 - len(unique_aa) + 1
					self.site_type = '%sD' % degeneracy
					# AA's identical: degeneracy = 4 - 1 + 1 = 4
					# AA's all different, degeneracy = 4 - 4 + 1 = 1
				return

	def fetch_ref_codon(self, gene):
		""" Fetch codon within gene for given site """
		# position of site in gene
		gene_pos = self.ref_pos-gene['start'] if gene['strand']=='+' else gene['end']-self.ref_pos
		# position of site in codon
		codon_pos = gene_pos % 3
		# gene sequence (oriented start to stop)
		ref_codon = gene['seq'][gene_pos-codon_pos:gene_pos-codon_pos+3]
		return ref_codon, codon_pos
						
	def store(self, species):
		""" Store data for GenomicSite in Species"""
		# snps_info
		atcg_counts = ','.join([str(_) for _ in self.pooled_counts])
		atcg_aas = ','.join([self.amino_acids[_] for _ in list('ATCG')]) if self.amino_acids != {} else ''
		info = '\t'.join([self.id, str(self.prevalence), self.ref_allele,
					self.major_allele, self.minor_allele, str(self.minor_freq),
					atcg_counts, self.site_type, atcg_aas, self.gene_id])+'\n'
		species.files['info'].write(info)
		# snps_freq
		freq = self.id + '\t' + '\t'.join([str(freq) for freq in self.sample_freqs])+'\n'
		species.files['freq'].write(freq)
		# snps_depth
		depth = self.id + '\t' + '\t'.join([str(depth) for depth in self.sample_depths])+'\n'
		species.files['depth'].write(depth)

def init_sample_offsets(species):
	""" Store sample.path, sample.offset, and sample.header
		species is list of Species objects
	    Species.samples is a list of Sample objects """
	for sample in species.samples:
		sample.path = '%s/snps/output/%s.snps.gz' % (sample.dir, species.id)
		file = gzip.open(sample.path)
		line = next(file)
		sample.offset = len(line)
		sample.header = line.rstrip('\n').split('\n')
		file.close()

def read_data(sample, nsites):
	""" Fetch block of nsites from sample """
	file = gzip.open(sample.path)
	file.seek(sample.offset)
	data = [next(file) for i in range(nsites)]
	offset = sample.offset + len(''.join(data))
	file.close()
	return offset, data

def load_data(species, sites, output):
	""" Loop over data from multiple samples and load into GenomicSite objects """
	# fetch site info from 1st sample
	for site, line in zip(sites, output[0]):
		values = line.split('\t')
		site.ref_id = values[0]
		site.ref_pos = int(values[1])
		site.ref_allele = values[2]
		site.id = values[0]+'|'+values[1]
	# fetch site counts from all sample
	for i in range(len(sites)):
		sites[i].sample_counts = [ [int(_) for _ in o[i].split('\t')[-1].split(',')] for o in output ]


#def load_data(species, sites, output):
#	""" Loop over data from multiple samples and load into GenomicSite objects 
#	
#		output: list of lists, each element
#	"""
#	# fetch site info from 1st sample
#	for site, line in zip(sites, output[0]):
#		values = line.split('\t')
#		site.ref_id = values[0]
#		site.ref_pos = int(values[1])
#		site.ref_allele = values[2]
#		site.id = values[0]+'|'+values[1]
#	# fetch site counts from all sample
#	for i in range(len(sites)):
#		sites[i].sample_counts = [ [int(_) for _ in o[i].split('\t')[-1].split(',')] for o in output ]


def parallel_read_data(species, nsites, threads):
	""" """
	pool = mp.Pool(processes=threads)
	results = [pool.apply_async(read_data, args=(s, nsites)) for s in species.samples]
	offsets = [result.get()[0] for result in results]
	data = [result.get()[1] for result in results]
	for offset, sample in zip(offsets, species.samples):
		sample.offset = offset
	return data

def open_outfiles(species, outdir):
	""" Open output files for species """
	species.outdir = os.path.join(outdir, species.id)
	species.files = {}
	for type in ['info', 'freq', 'depth']:
		species.files[type] = open('%s/snps_%s.txt' % (species.outdir, type), 'w')
	for type in ['freq', 'depth']:
		species.files[type].write('\t'.join(['site_id']+[s.id for s in species.samples])+'\n')
	info_fields = ['site_id', 'site_prev', 'ref_allele', 'major_allele', 'minor_allele', 'minor_freq', 'atcg_counts', 'site_type', 'atcg_aas', 'gene_id']
	species.files['info'].write('\t'.join(info_fields)+'\n')

def split_sites(genome_length, sites_per_iter, max_sites=float('Inf')):
	""" Create list of integers indicating the number of genomic sites to process at one time """
	max_sites = min(genome_length, max_sites)
	total_sites = 0
	site_batches = []
	while True:
		site_batches.append(sites_per_iter)
		total_sites += sites_per_iter
		if total_sites >= max_sites:
			site_batches[-1] -= (total_sites - max_sites)
			break
	return site_batches

def merge_snps(args, species):
	"""
	Merge nucleotide variants for specified species

	args = dictionary of command-line arguments:
	  'outdir' = path to output directory
	  'db' = path to midas database
	  'max_sites' = maximum number of sites to process
	  'site_depth' = minimum depth for calling a site present in a sample
	  'site_ratio' = minimum ratio of site_depth:sample depth
	  'site_maf' = minimum pooled minor allele frequency for site
	  'site_prev' = minimum fraction of samples where site passes filters
	  'sites_per_iter' = number of sites to store in memory at one time
	species = Species object:
	  species.id = species identifier
	  species.samples = list of Sample objects:
	    sample.dir = path to sample directory
		sample.id = sample identifier
		sample.info = dictionary of summary stats
	  sample.sample_depth = list of average read-depths of samples for species
	"""
		
	write_readme(species, args['outdir'], args['db'])
	species.write_sample_info(type='snps', outdir=args['outdir'])
	open_outfiles(species, args['outdir'])
	init_sample_offsets(species)
	
	contigs = read_genome(args['db'], species.id)
	genes, gene_index = read_genes(args['db'], species.id, contigs), [0]

	# sites are split into blocks to limit memory usage
	site_blocks = split_sites(species.genome_length, species.sites_per_iter, args['max_sites'])
	for nsites in site_blocks:

		# 1. read raw data from files in parallel (1 file per thread)
		# data = [ [l1_s1, l2_s1, ...], [l1_s2, l2_s2, ...] ]
		data = parallel_read_data(species, nsites, args['threads'])

		# 2. parse raw data in parallel (1 block of sites per thread)
		# each site block split up a 2 time
		# this is done to increase throughput
		
		
		### pick up here
		
		start_stop = []
		for msites in split_sites(nsites, nsites/args['threads']):
			if len(start_stop) == 0:
				start_stop.append([0,msites])
			else:
				start_stop.append([start_stop[-1][0]+msites, start_stop[-1][1]+msites])

		nsamples = len(species.samples)

		start, stop = start_stop[0]

		data_block = [_[start:stop] for _ in data]

		sites = [GenomicSite(nsamples) for i in range(stop-start)]

		for site, line in zip(sites, data_block[0]):
			values = line.split('\t')
			site.ref_id = values[0]
			site.ref_pos = int(values[1])
			site.ref_allele = values[2]
			site.id = values[0]+'|'+values[1]

		for m, site in enumerate(sites):
			site.sample_counts = [ [int(_) for _ in o[m].split('\t')[-1].split(',')] for o in data_block ]
			site.compute_pooled_counts()
			site.call_major_minor_alleles()
			site.compute_multi_freq()
			site.compute_minor_freqs()
			site.compute_prevalence(species, args['site_depth'], args['site_ratio'])
			site.filter(args['site_maf'], args['site_prev'], args['site_multi_freq'])
			if not site.remove:
				site.annotate(genes, gene_index, contigs)

		quit()

def id_sites_per_iter(num_samples, max_gb, max_sites_per_iter):
	""" Set the variable 'sites_per_iter' to keep memory below 'max_gb' 
		num_samples: number of samples in species
		max_gb: do not exceed this much ram
		max_sites_per_iter: maximum value for sites_per_iter 
	"""
	if not max_gb:
		return max_sites_per_iter
	else:
		overhead = 2e8 # 200 Mb
		max_bytes = max_gb * 1e9 - overhead  # max mem in bytes
		bytes_per_sample_per_site = 450  # num of bytes per sample per site
		bytes_per_site = bytes_per_sample_per_site * num_samples # num of bytes per site
		sites_per_iter = int(max_bytes/bytes_per_site)
		return min(max_sites_per_iter, sites_per_iter)


def run_pipeline(args):

	print("Identifying species")
	species = merge.select_species(args, data_type='snps')
	for sp in species:
		num_samples=len(sp.samples)
		sp.sites_per_iter = id_sites_per_iter(num_samples, args['max_gb'], args['sites_per_iter'])

	print("Merging snps")
	merge_snps(args, species[0])







			
