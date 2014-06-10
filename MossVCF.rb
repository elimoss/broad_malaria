#!/usr/bin/env ruby


#A collection of utility functions that are useful for other scripts.

GATK = '/humgen/gsa-hpprojects/GATK/bin/GenomeAnalysisTK-2.4-9-g532efad/GenomeAnalysisTK.jar'



	def make_intervals(i, n, fasta)
=begin
		Divide a sorted genome into n equally sized parts and return the i'th part.
		We will return a list of intervals: chr, start, stop.  It may contain multiple
		chromosomes in order to fill out a total length equal to the other parts.  Each
		part will be adjacent and non-overlapping with the next part.
		i must be a number from 1 to n.
=end
		i = i.to_i
		n = n.to_i
		
		raise "Index is beyond allowable range" unless 1 <= i and i<= n
		raise "Necessary partitioning argument(s) omitted" unless i != nil and n != nil
		
		# read genome dict file
		tot = 0
		chrlens = []
		raise "Invalid fasta file." unless fasta.end_with?('.fasta') or fasta.end_with?('.fa')
		File.open(fasta.split(".")[0] + '.dict', 'r').read.each_line("\n") do |rawline|
			row = rawline.split("\t")
			if row[0]=='@SQ'
				raise "Malformed fasta dictionary row" unless row[1].start_with?('SN:') and row[2].start_with?('LN:')
				c = row[1][3...row[1].length]
				c_len = row[2][3...row[2].length].to_i
				chrlens << [c,c_len,tot]
				tot += c_len
			end
		end
		
		# define our chunk by gpos:
		part_size = tot/n
		g_start = 1 + part_size * (i-1)
		g_stop = part_size * i
		g_stop = tot if i==n

		# find the genomic intervals that correspond to our gpos window
		out = []
		chrlens.each do |c, c_len, c_g_start|
			c_g_stop = c_g_start + c_len
			c_g_start += 1
			if c_g_stop >= g_start and c_g_start <= g_stop
				start = [g_start, c_g_start].max - c_g_start + 1
				stop  = [g_stop,  c_g_stop].min  - c_g_start + 1
				out << [c, start, stop]
			end
		end
		
		return out
	end
