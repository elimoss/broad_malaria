# util_vcf.py - This gives a number of useful quick methods for dealing with
#	VCF files.
#
# requires python >= 2.5
#
# dpark@broadinstitute.org
# $Id: util_vcf.py 7616 2013-04-12 18:35:23Z dpark $

import os, shutil, logging, sqlite3
import pysam		## must be installed by user
import util_files

log = logging.getLogger(__name__)

def make_intervals(i, n, fasta, chr_prefix='', verbose=False):
	''' Divide a sorted genome into n equally sized parts and return the i'th
		part.  We will return a list of intervals: chr, start, stop.  It may
		contain multiple chromosomes in order to fill out a total length equal
		to the other parts.  Each part will be adjacent and non-overlapping with
		the next part. i must be a number from 1 to n.
	'''
	assert 1 <= i <= n
	
	# read genome dict file
	tot = 0
	chrlens = []
	for c,c_len in get_chrlens(fasta):
		if c.startswith(chr_prefix):
			chrlens.append((c,c_len,tot))
			tot += c_len
	
	# define our chunk by gpos:
	part_size = tot/n
	g_start = 1 + part_size * (i-1)
	g_stop = part_size * i
	if i==n:
		g_stop = tot
	
	# find the genomic intervals that correspond to our gpos window
	out = []
	for c, c_len, c_g_start in chrlens:
		c_g_stop = c_g_start + c_len
		c_g_start += 1
		if c_g_stop >= g_start and c_g_start <= g_stop:
			start = max(g_start, c_g_start) - c_g_start + 1
			stop  = min(g_stop,  c_g_stop)  - c_g_start + 1
			out.append((c, start, stop))
	
	if verbose:
		log.info("Dividing the %d bp genome into %d chunks of %d bp each.  The %dth chunk contains the following %d intervals: %s" % (
			tot, n, part_size, i, len(out), ', '.join(["%s:%d-%d"%x for x in out])))
	return out

def vcf_bgzip_index(inVcf, outVcf, tabixPath=None, vcftoolsPath=None):
	assert (inVcf.endswith('.vcf') or inVcf.endswith('.vcf.gz')) and outVcf.endswith('.vcf.gz')
	''' Compress (bgzip) and index (tabix and/or vcftools) a VCF file.'''
	
	if inVcf.endswith('.vcf'):
		log.info("compressing output to %s" % outVcf)
		cmdline = "%s/bgzip -c %s > %s" % (tabixPath, inVcf, outVcf)
		assert not os.system(cmdline)
	elif inVcf==outVcf:
		log.info("leaving %s, already compressed" % inVcf)
	else:
		log.info("copying compressed vcf to %s" % outVcf)
		shutil.copy(inVcf, outVcf)
	
	log.info("indexing with tabix")
	cmdline = "%s/tabix %s -f -p vcf" % (tabixPath, outVcf)
	assert not os.system(cmdline)
	
	if vcftoolsPath:
		log.info("indexing with vcftools")
		tmpFile = util_files.mkstempfname(prefix='vcftools-log-', suffix='.vcf')
		cmdline = "%s/vcftools --gzvcf %s --out %s --force-index-write" % (vcftoolsPath, outVcf, tmpFile)
		assert not os.system(cmdline)
		os.unlink(tmpFile)
		os.unlink(tmpFile+'.log')
	return outVcf

class GenomePosition:
	''' Provide a mapping from chr:pos to genomic position.
		Read chromosome lengths and order from either a Picard/GATK-index for
		a FASTA file (a .dict file) or from a VCF header.
	'''
	def __init__(self, seqDb):
		self.gpos_map = {}
		self.chrs = []
		totlen = 0
		for c,clen in get_chrlens(seqDb):
			self.chrs.append((c,clen))
			self.gpos_map[c] = totlen
			totlen += clen
	def get_gpos(self, c, p):
		return p + self.gpos_map[c]
	def get_chr_pos(self, gpos):
		totlen = 0
		for c,clen in self.chrs:
			if gpos < totlen+clen:
				break
			totlen += clen
		return (c,gpos-totlen)

def get_chrlens(inFile):
	''' Read chromosome lengths and order from either a Picard/GATK-index for
		a FASTA file (a .dict file) or from "contig" rows in the VCF header.
	'''
	chrlens = []
	if hasattr(inFile, 'chrlens'):
		chrlens = inFile.chrlens()
	elif inFile.endswith('.fasta'):
		inFile = inFile[:-len('.fasta')] + '.dict'
	elif inFile.endswith('.fa'):
		inFile = inFile[:-len('.fa')] + '.dict'
	if inFile.endswith('.dict'):
		with open(inFile, 'rt') as inf:
			for line in inf:
				row = line.rstrip('\n').split('\t')
				if row[0]=='@SQ':
					assert row[1].startswith('SN:') and row[2].startswith('LN:')
					c = row[1][3:]
					c_len = int(row[2][3:])
					chrlens.append((c,c_len))
	elif inFile.endswith('.vcf') or inFile.endswith('.vcf.gz'):
		with util_files.open_or_gzopen(inFile, 'rt') as inf:
			for line in inf:
				line = line.rstrip('\n')
				if line.startswith('##contig=<ID=') and line.endswith('>'):
					line = line[13:-1]
					c = line.split(',')[0]
					clen = int(line.split('=')[1])
					chrlens.append((c,clen))
				elif line.startswith('#CHROM'):
					break
	else:
		raise AssertionError("unrecognized file type %s" % inFile)
	assert chrlens, "no sequence data found in %s % inFile"
	return chrlens

def get_chroms(inVcf, tabixPath=None):
	''' Get a list of unique chromosomes for this genome, in sort order.
		(Use tabix to do it quickly)
	'''
	tmpFile = util_files.mkstempfname(prefix='chrnames-', suffix='.txt')
	cmdline = "%s/tabix -l %s > %s" % (tabixPath, inVcf, tmpFile)
	assert not os.system(cmdline)
	chroms = []
	with open(tmpFile, 'rt') as inf:
		for line in inf:
			chroms.append(line.rstrip('\r\n'))
	os.unlink(tmpFile)
	return chroms

def vcf_sample_names(inVcf):
	''' Return the list of sample names in a given VCF file (quickly). '''
	samples = None
	with util_files.open_or_gzopen(inVcf, 'rt') as inf:
		for line in inf:
			if line.startswith('#CHROM'):
				row = line.rstrip('\r\n').split('\t')
				samples = row[9:]
				break  # we must properly close the file.. not sure if "return" does that
	assert samples, "error: no header line!"
	return samples

def vcf_subset(inVcf, c, start_stop=None, outVcf=None, keepHeader=False, tabixPath=None):
	''' Pull just a piece of a VCF file into a new VCF file (create a temp
		file if outVcf is not specified).
	'''
	if outVcf==None:
		outVcf = util_files.mkstempfname(prefix='vcf_subset-%s-'%c, suffix='.vcf')
	assert inVcf.endswith('.vcf.gz') and outVcf.endswith('.vcf')
	cmdline = "%s/tabix" % tabixPath
	if keepHeader:
		cmdline += ' -h'
	cmdline += ' %s %s' % (inVcf, c)
	if start_stop:
		cmdline += ':%d-%d' % start_stop
	cmdline += ' > %s' % outVcf
	assert not os.system(cmdline)
	return outVcf

def vcf_subset_rows(inVcf, c=None, start_stop=None, tabixPath=None):
	''' Pull just a piece of a VCF file and return as an iterator of 4-tuples '''
	if c:
		sub_vcf = vcf_subset(inVcf,c,start_stop=start_stop,keepHeader=True,tabixPath=tabixPath)
		for x in vcf_rows(sub_vcf):
			yield x
		os.unlink(sub_vcf)
	else:
		for x in vcf_rows(inVcf):
			yield x

def vcf_rows(inVcf):
	''' Read a VCF file and return contents as an iterator.
		Return each row as a 4-tuple:
			1: [first 9 columns]
			2: [all genotype columns (10th col and onward)]
			3: {dict of sample name : genotype column}
			4: [list of sample names]
	'''
	with util_files.open_or_gzopen(inVcf, 'rt') as inf:
		samples = []
		for line in inf:
			if not line.startswith('#'):
				row = line.rstrip('\n').split('\t')
				yield (row[:9],
					row[9:],
					dict([(samples[i],row[9+i]) for i in range(len(samples))]),
					samples)
			elif line.startswith('#CHROM'):
				samples = line.rstrip('\n').split('\t')[9:]

def vcf_haploid_iterator(inVcf, sample_list=None,
	drop_indels=False, drop_monomorphic=False, drop_multiallelic=False,
	interval=None, tabixPath=None):
	'''	Read a VCF file and return contents as an iterator with parsed
		contents (assuming haploid genotypes).  Each row is returned as
		a 4-tuple:
			1: chr
			2: pos (int)
			3: list of allele strings (in order)
			4: list of genotypes as 2-tuples: (sample, int allele)
	'''
	if tabixPath or interval==None:
		c, start_stop = (None,None)
		if interval:
			iparts = interval.split(':')
			c = iparts[0]
			if len(iparts)>1:
				start_stop = map(int, iparts[1].split('-'))
		for info,genos,genoMap,rowsamp in vcf_subset_rows(inVcf, c, start_stop=start_stop, tabixPath=tabixPath):
			c,p = info[:2]
			p = int(p)
			alleles = [info[3]] + info[4].split(',')
			if drop_monomorphic and len(alleles)==1:
				continue
			if drop_indels and not all([len(a)==1 for a in alleles]):
				continue
			if sample_list!=None:
				rowsamp = sample_list
			genos = [(s,int(genoMap[s][0])) for s in rowsamp if genoMap[s][0]!='.']
			n_alleles = len(set([a for s,a in genos]))
			if drop_monomorphic and n_alleles<2 or drop_multiallelic and n_alleles>2:
				continue
			yield (c, p, alleles, genos)
	else:
		with VcfReader(inVcf) as vcf:
			for c,p,alleles,genos in vcf.get_range(region=interval, as_strings=False):
				if drop_monomorphic and len(alleles)==1:
					continue
				if drop_indels and not all([len(a)==1 for a in alleles]):
					continue
				if sample_list!=None:
					genos = dict(genos)
					genos = [(s,genos[s]) for s in sample_list if s in genos]
				n_alleles = len(set([a for s,a in genos]))
				if drop_monomorphic and n_alleles<2 or drop_multiallelic and n_alleles>2:
					continue
				yield (c, p, alleles, genos)

class SnpDb:
	''' Present genotype data as a sqlite database.  Load from a VCF file. '''
	def __init__(self, inVcf,
		enforce_unique_pos=True, sample_list=None,
		drop_indels=False, drop_monomorphic=False, drop_multiallelic=False,
		interval=None):
		assert inVcf.endswith('.vcf.gz')
		clens = get_chrlens(inVcf)
		self.clens = dict(clens)
		self.contigs = [c for c,l in clens]
		self.sample_names = vcf_sample_names(inVcf)
		root_fname = inVcf[:-7].split('/')[-1]
		self.dbFile = util_files.mkstempfname(prefix='%s-'%root_fname,suffix='.db')
		self.conn = sqlite3.connect(self.dbFile, isolation_level='DEFERRED')
		self.cur = self.conn.cursor()
		self.cur.execute("""create table snp (
			chr string not null,
			pos integer not null,
			alleles string not null)""")
		self.cur.execute("""create table geno (
			chr string not null,
			pos integer not null,
			sample string not null,
			allele integer not null)""")
		self.cur.execute("create %s index snp_idx on snp(chr,pos)" % (
			enforce_unique_pos and "unique" or ""))
		self.cur.execute("create %s index geno_idx on geno(chr,pos,sample)" % (
			enforce_unique_pos and "unique" or ""))
		self.conn.commit()
		self.cur.executemany("insert into snp (chr,pos,alleles) values (?,?,?)",
			[(c,p,','.join(alleles))
				for c,p,alleles,genos in vcf_haploid_iterator(inVcf,
					sample_list=sample_list,
					drop_indels=drop_indels, drop_monomorphic=drop_monomorphic,
					drop_multiallelic=drop_multiallelic,
					interval=interval)])
		self.conn.commit()
		self.cur.executemany("insert into geno (chr,pos,sample,allele) values (?,?,?,?)",
			self._geno_iterator(vcf_haploid_iterator(inVcf, sample_list=sample_list,
				drop_indels=drop_indels, drop_monomorphic=drop_monomorphic,
				drop_multiallelic=drop_multiallelic,
				interval=interval)))
		self.conn.commit()
	def _geno_iterator(self, vcf_iterator):
		for c,p,alleles,genos in vcf_iterator:
			for s,a in genos:
				yield (c,p,s,a)
	def close(self):
		self.cur.close()
		self.conn.close()
		os.unlink(self.dbFile)
	def __enter__(self):
		return self
	def __exit__(self, exc_type, exc_val, exc_tb):
		self.close()
		return 0
	def get_conn(self):
		return self.conn
	def chroms(self):
		return self.contigs
	def samples(self):
		return self.sample_names
	def chrlens(self):
		return self.clens
	def get_snp_genos(self, c, p, as_strings=True):
		self.cur.execute("select sample,allele from geno where chr=? and pos=?", [c,p])
		out = [(s,a) for s,a in self.cur]
		if as_strings and out:
			self.cur.execute("select alleles from snp where chr=? and pos=?", [c,p])
			alleles = [x[0] for x in self.cur]
			if not len(alleles)==1:
				out = []
			else:
				alleles = alleles[0].split(',')
				out = [(s,alleles[a]) for s,a in out]
		return dict(out)


class SnpDbOneSamp:
	''' Present genotype data as a sqlite database.  Load from a VCF file. '''
	def __init__(self, inVcf,
		enforce_unique_pos=True, sample=None, interval=None):
		assert inVcf.endswith('.vcf.gz')
		clens = get_chrlens(inVcf)
		self.clens = dict(clens)
		self.contigs = [c for c,l in clens]
		self.sample_names = vcf_sample_names(inVcf)
		root_fname = inVcf[:-7].split('/')[-1]
		self.dbFile = util_files.mkstempfname(prefix='%s-'%root_fname,suffix='.db')
		self.conn = sqlite3.connect(self.dbFile, isolation_level='DEFERRED')
		self.cur = self.conn.cursor()
		self.cur.execute("""create table cons (
			chr string not null,
			pos integer not null,
			allele string not null)""")
		self.cur.execute("create %s index cons_idx on cons(chr,pos)" % (
			enforce_unique_pos and "unique" or ""))
		self.conn.commit()
		self.cur.executemany("insert into cons (chr,pos,allele) values (?,?,?)",
			[(c,p,alleles[genos[0][1]])
				for c,p,alleles,genos in vcf_haploid_iterator(inVcf,
					sample_list=sample and [sample] or None,
					interval=interval)])
		self.conn.commit()
	def close(self):
		self.cur.close()
		self.conn.close()
		os.unlink(self.dbFile)
	def __enter__(self):
		return self
	def __exit__(self, exc_type, exc_val, exc_tb):
		self.close()
		return 0
	def get_conn(self):
		return self.conn
	def chroms(self):
		return self.contigs
	def samples(self):
		return self.sample_names
	def chrlens(self):
		return self.clens
	def get_geno(self, c, p):
		self.cur.execute("select allele from cons where chr=? and pos=?", [c,p])
		out = [a[0] for a in self.cur]
		return len(out)==1 and out[0] or None


class TabixReader(pysam.Tabixfile):
	''' A wrapper around pysam.Tabixfile that provides a context and
		allows us to query using 1-based coordinates.
	'''
	def __init__(self, inFile, parser=pysam.asTuple()):
		pysam.Tabixfile.__init__(self, inFile)
		self.parser = parser
	def __enter__(self):
		return self
	def __exit__(self, exc_type, exc_val, exc_tb):
		self.close()
		return 0
	def close(self):
		if hasattr(pysam.Tabixfile, 'close'):
			pysam.Tabixfile.close(self)
		else:
			log.warn("warning: pysam-0.5 lacks a pysam.Tabixfile.close() method.  The input file may not be closed.")
			# NOTE: this does not exist in pysam-0.5 (the Broad's installed version)
			# I'm not sure if that means it actually doesn't close the file ever.
	def chroms(self):
		return self.contigs
	def get(self, chrom=None, start=None, stop=None, region=None):
		if start!=None:
			start -= 1
		return pysam.Tabixfile.fetch(self, reference=chrom, start=start, end=stop,
			region=region, parser=self.parser)

class VcfReader(TabixReader):
	''' Same as TabixReader with a few more perks for VCF files:
		- emit results parsed as pysam VCF rows
		- provide self.chrlens(), a dict mapping chrom name to chrom length
		- provide self.samples(), a list of sample names in order of appearance
		- provide get_range(c,start,stop) and get_snp_genos(c,pos)
	'''
	def __init__(self, inFile, ploidy=1, parser=pysam.asVCF()):
		TabixReader.__init__(self, inFile, parser=parser)
		assert ploidy in (1,2)
		self.ploidy = ploidy
		self.clens = []
		self.sample_names = None
		for line in self.header:
			if line.startswith('##contig=<ID=') and line.endswith('>'):
				line = line[13:-1]
				c = line.split(',')[0]
				clen = int(line.split('=')[1])
				self.clens.append((c,clen))
			elif line.startswith('#CHROM'):
				row = line.split('\t')
				self.sample_names = row[9:]
		self.clens = dict(self.clens)
		assert self.sample_names
	def samples(self):
		return self.sample_names
	def chrlens(self):
		return self.clens
	def get_range(self, c=None, start=None, stop=None, region=None, as_strings=True):
		'''	Read a VCF file (optionally just a piece of it) and return contents
			as an iterator with parsed contents.  Each row is returned as
			a 4-tuple:
				1: chr
				2: pos (int)
				3: list of allele strings (in order)
				4: list of genotypes as 2-tuples:
					haploid: (sample, allele)
					diploid: (sample, [allele, allele])
			If as_strings, the alleles will be actual alleles.  Otherwise,
			alleles will be integers.
		'''
		for snp in self.get(c,start,stop,region):
			alleles = [snp.ref] + snp.alt.split(',')
			if self.ploidy==1:
				genos = [(self.sample_names[i], int(snp[i][0]))
					for i in range(len(self.sample_names))
					if snp[i][0] != '.']
				if as_strings:
					genos = [(s,alleles[a]) for s,a in genos]
			else:
				genos = [(self.sample_names[i], [int(snp[i][j*2]) for j in range(self.ploidy)])
					for i in range(len(self.sample_names))
					if snp[i][0] != '.']
				if as_strings:
					genos = [(s,[alleles[a] for a in g]) for s,g in genos]
			yield (snp.contig, snp.pos, alleles, genos)
	def get_snp_genos(self, c, p, as_strings=True):
		''' Read a single position from a VCF file and return the genotypes
			as a sample -> allele map.  If there is not exactly one mathing
			row in the VCF file at this position (if there are none or multiple)
			then we return an empty map: {}.
		'''
		snps = [x for x in self.get_range(c,p,p,as_strings=as_strings)]
		return len(snps)==1 and dict(snps[0][3]) or {}
	def getFullSequence(self, chr, start, stop=None, sample=None,
		na='-', refSeq=None, ignoreIndels=False):
		''' chr - chromosome name
			start - start position
			stop - default = start
			sample - default = REF allele
			if refSeq is a string with length = stop-start+1, then use this as the
				base sequence for all missing data, otherwise, use the na string
		'''
		assert chr in self.chroms()
		assert sample==None or sample in self.samples()
		if stop==None:
			stop = start
		assert 1<=start<=stop
	
		if refSeq:
			assert len(refSeq)==(stop-start+1)
			seq = list(refSeq)
		else:
			seq = list(na * (stop-start+1))
	
		for c,p,alleles,genos in self.get_range(chr, start, stop, as_strings=True):
			i = p-start
			if sample==None:
				allele = alleles[0]
				alleles = [allele]
			else:
				fields = row[8].split(':')
				samp_geno = dict(genos).get(sample)
				if not samp_geno:
					continue
				if not isinstance(samp_geno, basestring):
					log.warn("TO DO: add code to turn hets into IUPAC ambiguity codes (%s:%s %s = %s)." % (c,p,sample, '/'.join(samp_geno)))
					continue
				allele = samp_geno
			if allele == None:
				# nothing here
				pass
			elif ignoreIndels and (len(alleles[0])!=1 or len(allele)!=1):
				# skip this line because we were told to
				pass
			elif len(alleles[0])==1:
				# the most common cases: SNP, novar, or indel replacing one ref base
				seq[i] = allele
			else:
				# more complicated case: something replacing multiple ref bases
				# TODO: beware--this is all untested!
				for j in range(max(len(alleles[0]), len(allele))):
					if i+j < (stop-start+1):
						if j<len(allele):
							if j==len(alleles[0])-1:
								# new allele is >= ref length, fill out the rest of the bases
								seq[i+j] = allele[j:]
							else:
								seq[i+j] = allele[j]
						else:
							# new allele is shorter than ref, so delete extra bases
							seq[i+j] = ''
		seq = ''.join(seq)
		assert len(seq)==(stop-start+1) or not ignoreIndels
		return seq

class TabixWriter:
	''' A file-writing context that writes tab delimited text to a temp file
		and, upon exit, bgzip compresses to its final destination and tabix
		indexes it.
		tabixPath is only required if
			- header_prefix='' and there are header lines
			- or if you don't have pysam installed
	'''
	def __init__(self, outFile, header_prefix='#',
			col_chr=1, col_start=2, col_stop=2, preset=None, zerobased=False,
			tabixPath=None):
		assert outFile.endswith('.gz')
		assert not os.access(outFile, os.F_OK) or os.access(outFile, os.W_OK)
		self.outFile = outFile
		root_fname = outFile.split('/')[-1].split('.')[0]
		self.tmpFile = util_files.mkstempfname(prefix='tabix-%s-'%root_fname,suffix='.txt')
		self.tmp_f = open(self.tmpFile, 'wt')
		self.params = {'seq_col':col_chr, 'start_col':col_start, 'end_col':col_stop,
			'preset':preset, 'zerobased':zerobased, 'meta_char':header_prefix,
			'line_skip':0}
		self.header = None
		self.n_header = 0
		self.n_rows = 0
		self.tabixPath=tabixPath
		self.last_chr_start = [None,None]
		self.seen_chrs = set()
	def __enter__(self):
		return self
	def __exit__(self, exc_type, exc_val, exc_tb):
		self.close(error=(exc_type!=None))
		return 0
	def close(self, error=False):
		close(self.tmp_f)
		if not error:
			if self.params['meta_char']=='':
				self.params['line_skip'] = self.n_header
			bgzip_index(self.tmpFile, self.outFile, params, tabixPath=self.tabixPath)
			os.unlink(self.tmpFile)
	def write_header(self, header):
		''' header is a list of strings to be joined by tabs '''
		assert not self.n_rows, "cannot write more header again after data has started"
		self.header = header
		self.tmp_f.write(self.params['meta_char']+'\t'.join(header)+'\n')
		self.n_header += 1
	def write_list(self, row):
		''' row is a list of data to be converted to strings and joined by tabs '''
		self.validate(row)
		self.tmp_f.write('\t'.join([x!=None and str(x) or '' for x in row])+'\n')
		self.n_rows += 1
	def write_map(self, row):
		''' row is a dict of data which we key into with the names of the last
			header row. '''
		self.write_list([row.get(h) for h in self.header])
	def validate(self, row):
		c, start, stop = [row[self.params[x]-1] for x in ('seq_col','start_col','end_col')]
		start,stop = (int(start), int(stop))
		assert stop>=start, "error: stop (%d) < start (%s) on chr %s" % (stop,start,c)
		if self.last_chr_start[0]==c:
			assert start>=self.last_chr_start[1], "error: data must be sorted, %s:%d seen after %s:%d" %(c,start,self.last_chr_start[0],self.last_chr_start[1])
		else:
			assert c not in self.seen_chrs, "error: data must be sorted, %s seen in discontinuous order" % c
			self.seen_chrs.add(c)
		self.last_chr_start = [c,start]
		


def bgzip_index(inFile, outFile, params={}, tabixPath=None):
	assert not inFile.endswith('.gz') and outFile.endswith('.gz')
	
	log.debug("compressing with bgzip %s -> %s" % (inFile, outFile))
	if tabixPath:
		cmdline = "%s/bgzip -c %s > %s" % (tabixPath, inFile, outFile)
		assert not os.system(cmdline)
	else:
		pysam.tabix_compress(self.tmpFile, self.outFile, force=True)
	
	log.debug("indexing with tabix: %s" % outFile)
	if tabixPath:
		cmdline = "%s/tabix %s -f" % (tabixPath, outFile)
		if params.get('seq_col')!=None:
			cmdline += ' -s %d' % params['seq_col']
		if params.get('start_col')!=None:
			cmdline += ' -b %d' % params['start_col']
		if params.get('end_col')!=None:
			cmdline += ' -e %d' % params['end_col']
		if params.get('preset')!=None:
			cmdline += ' -p %s' % params['preset']
		if params.get('meta_char')!=None:
			cmdline += ' -c %s' % params['meta_char']
		if params.get('line_skip')!=None:
			cmdline += ' -S %d' % params['line_skip']
		if params.get('zerobased')!=None:
			cmdline += ' -0'
		assert not os.system(cmdline)
	else:
		assert not params.get('line_skip'), "error: pysam does not seem to support this option, even though their documentation talks about it"
		pysam.tabix_index(self.outFile, force=True,
			seq_col=params.get('seq_col'),
			start_col=params.get('start_col'), end_col=params.get('end_col'),
			preset=params.get('preset'), meta_char=params.get('meta_char','#'),
			zerobased=params.get('zerobased',False))
	return outFile
