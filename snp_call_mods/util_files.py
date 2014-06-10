# util_files.py - This gives a number of useful quick methods for dealing with
#	tab-text files and gzipped files.
#
# requires python >= 2.5
#
# dpark@broadinstitute.org
# $Id: util_files.py 7616 2013-04-12 18:35:23Z dpark $

import os, gzip, tempfile, logging

log = logging.getLogger(__name__)

def histogram(items):
	''' I count the number of times I see stuff and return a dict of counts. '''
	out = {}
	for i in items:
		out.setdefault(i, 0)
		out[i] += 1
	return out

def unique(items):
	''' Return unique items in the same order as seen in the input. '''
	seen = set()
	out = []
	for i in items:
		if i not in seen:
			seen.add(i)
			out.append(i)
	return out

def mkstempfname(suffix='', prefix='tmp', dir=None, text=False):
	''' There's no other one-liner way to securely ask for a temp file by filename only.
		This calls mkstemp, which does what we want, except that it returns an open
		file handle, which causes huge problems on NFS if we don't close it.  So close
		it first then return the name part only.
	'''
	fd, fn = tempfile.mkstemp(prefix=prefix, suffix=suffix, dir=dir, text=text)
	os.close(fd)
	return fn

def open_or_gzopen(fname, mode):
	return fname.endswith('.gz') and GZFile(fname, mode) or open(fname, mode)

class GZFile(gzip.GzipFile):
	''' This is needed for python 2.5 and 2.6.  Python 2.7 natively supports
		this via gzip.GzipFile. '''
	def __enter__(self):
		return self
	def __exit__(self, exc_type, exc_val, exc_tb):
		self.close()
		return 0

def cat_files_and_header(inFiles, outFile, add_header=None, strip_in_rows=0):
	''' This is like cat... except that we optionally drop header rows from
		the input files and/or optionally add a new header row to the output
		file.
	'''
	with open_or_gzopen(outFile, 'wt') as outf:
		if add_header:
			outf.write(add_header+'\n')
		for f in inFiles:
			with open_or_gzopen(f, 'rt') as inf:
				for i in range(strip_in_rows):
					inf.readline()
				for line in inf:
					outf.write(line)
	return outFile

def readFlatFileHeader(filename, headerPrefix='#', delim='\t'):
	with open_or_gzopen(filename, 'rt') as inf:
		header = inf.readline().rstrip('\n').split(delim)
	if header and header[0].startswith(headerPrefix):
		header[0] = header[0][len(headerPrefix):]
	return header

class FlatFileParser:
	'''	Generic flat file parser that parses tabular text input	'''
	def __init__(self, lineIter=None, name=None, outType='dict', readHeader=True,
		headerPrefix='#', delim='\t', requireUniqueHeader=False):
		self.lineIter = lineIter
		self.header = None
		self.name = name
		self.headerPrefix = headerPrefix
		self.readHeader = readHeader
		self.delim = delim
		self.requireUniqueHeader = requireUniqueHeader
		self.line_num=0
		assert outType in ('dict','arrayStrict', 'arrayLoose','both')
		self.outType=outType
		assert readHeader or outType in ('arrayStrict', 'arrayLoose')
	def __enter__(self):
		return self
	def __exit__(self, exc_type, exc_val, exc_tb):
		return 0
	def __iter__(self):
		assert self.lineIter
		for row in self.lineIter:
			out = self.parse(row)
			if out != None:
				yield out
	def parse(self, row):
		self.line_num += 1
		try:
			line = row.rstrip('\n').split(self.delim)
			if self.readHeader:
				if self.headerPrefix and row.startswith(self.headerPrefix):
					line[0] = line[0][len(self.headerPrefix):]
					assert not (self.requireUniqueHeader and self.header)
					self.parseHeader(line)
					return None
				elif not self.header:
					self.parseHeader(line)
					return None
				else:
					return self.parseRow(line)
			else:
				return self.parseRow(line)
		except Exception as e:
			if self.name:
				log.exception("Exception parsing file at line %d.  Line contents: '%s'  File: %s" % (self.line_num, row, self.name))
			else:
				log.exception("Exception parsing file at line %d.  Line contents: '%s'" % (self.line_num, row))
			raise
	def parseHeader(self, row):
		assert row
		self.header = row
		if self.outType != 'arrayLoose':
			assert len(row) == len(dict([(x,0) for x in row]))
	def parseRow(self, row):
		assert self.outType=='arrayLoose' or self.header and len(self.header)==len(row)
		if self.outType=='arrayLoose' or self.outType=='arrayStrict':
			return row
		out = dict([(self.header[i],row[i]) for i in range(len(self.header))])
		if self.outType=='both':
			for i in range(len(self.header)):
				out[i] = row[i]
		return out

class FlatFileMaker:
	def __init__(self, header, mapIter, comment='#', delim='\t', nullchar='', noHeaderline=False):
		self.comment = comment
		self.delim = delim
		self.header = header
		self.nullchar = nullchar
		self.maps = mapIter
		self.noHeaderline = noHeaderline
	def __iter__(self):
		if not self.noHeaderline:
			yield self.comment + self.delim.join(self.header) + "\n"
		for x in self.maps:
			out = [str(x.get(h,self.nullchar)) for h in self.header]
			yield self.delim.join(out) + "\n"

class FastaMaker:
	def __init__(self, seqIter=[], linewidth=60):
		self.seqs = seqIter
		self.linewidth = linewidth
		assert linewidth>0
	def add(self, id, seq):
		self.seqs.append((id,seq))
	def __iter__(self):
		for id,seq in self.seqs:
			yield ">"+id+"\n"
			while len(seq)>self.linewidth:
				line = seq[:self.linewidth]
				seq = seq[self.linewidth:]
				yield line+"\n"
			if seq:
				yield seq+"\n"

class IteratorMerger:
	def __init__(self, iterIter=[]):
		self.queue = iterIter
	def addIter(self, iter):
		self.queue.append(iter)
	def addItem(self, item):
		self.queue.append([item])
	def __iter__(self):
		for i in self.queue:
			for j in i:
				yield j

