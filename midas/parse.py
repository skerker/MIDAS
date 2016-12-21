#!/usr/bin/env python

# MIDAS: Metagenomic Intra-species Diversity Analysis System
# Copyright (C) 2015 Stephen Nayfach
# Freely distributed under the GNU General Public License (GPLv3)

import sys, os, gzip, numpy as np, random, csv
from midas.utility import print_copyright

class Sample:
	""" Base class for sample """
	def __init__(self, info):
		self.id = info['sample_id']
		self.info = info
		self.mean_depth = float(self.info['mean_coverage'])
		self.fract_cov = float(self.info['fraction_covered'])

	def filter(self, mean_depth, fract_cov):
		if self.fract_cov < fract_cov:
			return True
		elif self.mean_depth < mean_depth:
			return True
		else:
			return False

class Species:
	""" Base class for species """
	
	def __init__(self, dir):
		
		self.dir = dir
		self.id = os.path.basename(self.dir)
		self.init_paths()
		self.init_files()
		self.init_samples()

	def init_paths(self):
		self.paths = {}
		for type in ['freq', 'depth', 'info', 'summary']:
			self.paths[type] = '%s/snps_%s.txt' % (self.dir, type)

	def init_files(self):
		self.files = {}
		for type in ['freq', 'depth', 'info', 'summary']:
			file = open(self.paths[type])
			if type in ['info', 'summary']:
				self.files[type] = csv.DictReader(file, delimiter='\t')
			else:
				self.files[type] = csv.reader(file, delimiter='\t')

	def init_samples(self):
		self.sample_ids = None
		for file in ['freq', 'depth']:
			self.sample_ids = next(self.files[file])[1:]


class GenomicSite:
	""" Base class for genomic sites """
	def __init__(self, species, samples):
		try:
			# fetch site info
			info = next(species.files['info'])
			self.id = info['site_id']
			self.ref_allele = info['ref_allele']
			self.minor_allele = info['minor_allele']
			self.major_allele = info['major_allele']
			self.gene_id = info['gene_id']
			self.site_type = info['site_type']
						
			# copy samples
			self.samples = samples
			
			# fetch site data from freq and depth matrixes
			self.fetch_row(species)

		except StopIteration:
			self.id = None
		
	def fetch_row(self, species):
		""" Fetch next row from freq and depth matrices
		    Store in sample objects: sample.freq, sample.depth """
		freqs = next(species.files['freq'])[1:]
		depths = next(species.files['depth'])[1:]
		for sample in self.samples.values():
			self.samples[sample.id].freq = float(freqs[sample.index])
			self.samples[sample.id].depth = int(depths[sample.index])
				
	def flag_samples(self, site_depth, site_ratio):
		""" Filter samples at site based on coverage
			Set flag: sample.keep = [True/False] """
		for sample in self.samples.values():
			if (sample.depth >= site_depth and
					sample.depth/sample.mean_depth <= site_ratio):
				sample.keep = True
			else:
				sample.keep = False

	def filter(self, site_prev=None, site_maf=None, site_type=None):
		""" determine if site passes quality control """
		if self.ref_allele not in ['A','T','C','G']:
			self.keep = False
		if site_prev and self.prevalence < max(1e-6, site_prev):
			self.keep = False
		elif site_maf and self.minor_freq < site_maf:
			self.keep = False
		elif site_type and self.site_type not in site_type:
			self.keep = False
		else:
			self.keep = True

	def compute_prevalence(self):
		self.count_samples = sum([1 for s in self.samples.values() if s.keep])
		return float(self.count_samples)/len(self.samples)

	def summary_stats(self):
		""" compute summary stats for site across samples """
		self.prevalence = self.compute_prevalence()
		self.pooled_maf = self.compute_pooled_maf()

	def compute_pooled_maf(self):
		""" compute average frequency of reference allele with optional weighting of samples """
		if self.count_samples == 0:
			return 0.0
		else:
			return np.mean([s.freq for s in self.samples.values() if s.keep])

	def resample_reads(self, rand_reads, replace_reads):
		""" resample random number of reads per sample """
		index = 0
		for sample in self.samples.values():
			sample.depth = rand_reads
			if sample.freq > 0 and sample.freq < 1:
				count_minor = int(round(sample.freq * sample.depth))
				count_major = sample.depth - count_minor
				alleles = np.random.choice([1]*count_minor+[0]*count_major, rand_reads, replace=replace_reads)
				sample.freq = np.mean(alleles)
			index += 1

def fetch_samples(species, mean_depth=0, fract_cov=0, max_samples=float('inf'),
                  keep_samples=None, exclude_samples=None, rand_samples=None):
	""" Select samples from input
		mean_depth: filter samples based on average genome/gene depth
		fract_cov: filter samples based on the number of genomic sites covered by >=1 read
		max_samples: stop when this numbe of samples reached
		keep_samples: only keep samples in this list
		exclude_samples: exclude any sample in this list
		rand_samples: select random subset of samples
	"""
	samples = {}
	# select samples that pass filters
	for index, info in enumerate(species.files['summary']):
		sample = Sample(info)
		sample.index = index
		if sample.filter(mean_depth, fract_cov):
			continue
		if keep_samples and sample.id not in keep_samples:
			continue
		if exclude_samples and sample.id in exclude_samples:
			continue
		if len(samples) >= max_samples:
			continue
		samples[sample.id] = sample
	# check there's at least 1 sample
	if len(samples) == 0:
		error = "\nError: no samples satisfied your selection criteria.\n"
		error += "Try running again with more lenient parameters\n"
		sys.exit(error)
	# select random subset of samples
	if rand_samples:
		if rand_samples > len(samples):
			error = "\nError: --rand_samples cannot exceed the number of samples\n"
			sys.exit(error)
		ids = np.random.choice(samples.keys(), rand_samples, replace=False)
		for id in samples.keys():
			if id not in ids:
				del samples[id]
	return samples


def fetch_sites(species, samples):
	""" yield genomic sites from species across samples """
	index = 0
	while True:
		site = GenomicSite(species, samples)
		if not site.id:
			break
		else:
			index += 1
			yield site

