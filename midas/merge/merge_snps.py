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

def read_genes(species_id, db):
	""" Read in gene coordinates from features file """
	contigs = read_genome(db, species_id)
	genes = []
	path = '%s/rep_genomes/%s/genome.features.gz' % (db, species_id)
	for gene in utility.parse_file(path):
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
	def __init__(self, line):
	
		# initialize
		values = line.split('\t')
		self.id = values[0]+'|'+values[1]
		self.ref_id = values[0]
		self.ref_pos = int(values[1])
		self.ref_allele = values[2]
		
		# site annotations
		self.site_type = None
		self.gene_id = None
		self.amino_acids = None
		self.ref_codon = None
		self.codon_pos = None
		
		# site info
		self.major_allele = None
		self.minor_count = None
		self.minor_freq = None
		self.prevalence = None
		
		# pooled statistics
		self.alleles = {'A':0, 'T':1, 'C':2, 'G':3} # <dic>, mapping of allele to list index
		self.pooled_counts = [0,0,0,0] # <list>, count of each allele across all samples
		
		# per-sample statistics
		self.sample_counts = [] # <list>, [A,T,C,G] counts per sample
		self.sample_freqs = []  # <list>, minor allele frequencies
		self.sample_depths = [] # <list>, count of major+minor allele frequencies
				
	def compute_pooled_counts(self):
		""" Compute pooled nucleotide statistics (counts, depth, freq) at GenomicSite """
		for counts in self.sample_counts:
			for i in range(4):
				self.pooled_counts[i] += counts[i]

	def call_major_minor_alleles(self):
		""" Call major and minor alleles at GenomicSite """
		if sum(self.pooled_counts) > 0:
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
		counts = 0
		for count, allele in zip(self.pooled_counts, list('ATCG')):
			if not allele in [self.major_allele, self.minor_allele]:
				counts += count
		if counts > 0:
			self.multi_freq = counts/float(sum(self.pooled_counts))
		else:
			self.multi_freq = 0.0

	def compute_minor_freqs(self):
		""" Compute per-sample depth (major + minor allele) and minor allele freq at GenomicSite """
		for counts in self.sample_counts:
			count_major = counts[self.alleles[self.major_allele]]
			count_minor = counts[self.alleles[self.minor_allele]]
			sample_depth = count_major + count_minor
			sample_freq = float(count_minor)/sample_depth if sample_depth > 0 else 0.0
			self.sample_freqs.append(sample_freq)
			self.sample_depths.append(sample_depth)

	def compute_prevalence(self, mean_depths, min_depth, max_ratio):
		""" Compute the fraction of samples where site passes all filters """
		pass_qc = []
		for mean_depth, site_depth in zip(mean_depths, self.sample_depths):
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

	def annotate(self, genes, gene_index):
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
						
	def write(self, species):
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

def open_infiles(species):
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

def parallel_read_data(species, nsites, threads):
	""" Read data in parallel """
	# read nsites from each sample in parallel
	pool = mp.Pool(processes=threads)
	results = [pool.apply_async(read_data, args=(s, nsites)) for s in species.samples]
	pool.close()
	data = [result.get()[1] for result in results]
	
	# update file offset for each sample
	offsets = [result.get()[0] for result in results]
	for offset, sample in zip(offsets, species.samples):
		sample.offset = offset
		
	# split data into batches of msites for parallel processing
	data_blocks = []
	starts, stops = fetch_start_stops(nsites, threads)
	for start, stop in zip(starts, stops):
		data_blocks.append([sample[start:stop] for sample in data])

	return data_blocks

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

def fetch_start_stops(nsites, threads):
	starts, stops = [], []
	for msites in split_sites(nsites, nsites/threads):
		if len(starts) == 0:
			starts.append(0)
			stops.append(msites)
		else:
			starts.append(starts[-1]+msites)
			stops.append(stops[-1]+msites)
	return starts, stops

def parse_counts(site_index, data_block):
	counts = [ [int(_) for _ in sample[site_index].split('\t')[-1].split(',')] for sample in data_block ]
	return counts
		
def process_sites(data_block, mean_depths, args):
	sites = [GenomicSite(_) for _ in data_block[0]]
	for site_index, site in enumerate(sites):
		site.sample_counts = parse_counts(site_index, data_block)
		site.compute_pooled_counts()
		site.call_major_minor_alleles()
		site.compute_multi_freq()
		site.compute_minor_freqs()
		site.compute_prevalence(mean_depths, args['site_depth'], args['site_ratio'])
		site.filter(args['site_maf'], args['site_prev'], args['site_multi_freq'])
	return sites

def parallel_process_data(species, data, args):
	pool = mp.Pool(processes=args['threads'])
	results = [pool.apply_async(process_sites, args=(msites, species.sample_depth, args)) for msites in data]
	pool.close()
	sites = []
	for result in results:
		sites += result.get()
	return sites

def merge_snps(args, species):
	""" Main function to merge nucleotide variants for specified species """
		
	# write readme and sample info files
	write_readme(species, args['outdir'], args['db'])
	species.write_sample_info(type='snps', outdir=args['outdir'])
	
	# setup input and output files
	open_infiles(species)
	open_outfiles(species, args['outdir'])
	
	# read in list of genes for SNP annotation
	# genes = [gene1, gene2]; gene1 = {'start':int, 'end':int, 'seq':str}
	genes = read_genes(species.id, args['db'])
	gene_index = [0]

	# split sites into blocks to limit peak memory usage
	sites_per_iter = id_sites_per_iter(species.nsamples, args['max_gb'], args['sites_per_iter'])
	for nsites in split_sites(species.genome_length, sites_per_iter, args['max_sites']):

		# read nsites from each sample in parallel (1 thread per file)
		data = parallel_read_data(species, nsites, args['threads'])

		# parse raw data in parallel (1 block of sites per thread)
		sites = parallel_process_data(species, data, args)
	
		# annotate sites and write to output files
		for site in sites:
			site.annotate(genes, gene_index)
			site.write(species)

def run_pipeline(args):

	print("Identifying species")
	species = merge.select_species(args, data_type='snps')
		
	print("Merging snps")
	for sp in species:
		merge_snps(args, sp)

