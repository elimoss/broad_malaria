#!/usr/bin/env ruby

def coalesce(haplos_with_reads, min_ratio, verbose = false, output_fasta)
	haplos = haplos_with_reads.transpose.second.map{|t| t.split('')}
	haplo_list = haplos.uniq
	
	dists = Hash.new
	haplo_list.each_with_index{|h1, idx1| haplo_list.each_with_index{|h2, idx2| dists[[idx1, idx2]] = hamming_distance(h1, h2)}}

	ratios = Hash.new
	haplos_with_reads.each_with_index{|h1, idx1| 
		haplos_with_reads.each_with_index{|h2, idx2| 
			c1 = h1.first.to_i
			c2 = h2.first.to_i
			ratios[[idx1, idx2]] = [c1, c2]
		}
	}
=begin
#	print ratios
#	puts
	print dists
	puts
=end
#	puts "min haplo: #{"
	while(true)
		#coordinates to coalesce

		#the next to coalesce will be chosen as the two haplotypes having the least imbalanced abundance ratio of those that pass the filtration step.  The filtration step is intended to remove from consideration haplotype pairs satisfying any of the following: identity pairs, pairs more than hamming distance one apart, and pairs where either member is hamming distance 1 away from multiple other haplotypes (to encourage coalescence of sequential errors and avoid amputation of such series from the originating haplotype and subsequent creation of a spurious haplotype)
		
		ratios_filtered_keys = ratios.keys.select{|a,b| dists[[a,b]] == 1 and a != b and dists.keys.select{|k| k.include?(a)}.map{|k| dists[k]}.count(1) == 2}
#		print ratios_filtered_keys
#		puts
		ratios_computed = ratios_filtered_keys.map{|a,b| ratios[[a,b]].min/ratios[[a,b]].max.to_f}
		
		break if ratios_computed.empty? or ratios_computed.min > min_ratio #nothing left that's imbalanced enough
		
		
		x,y = ratios_filtered_keys[ratios_computed.index(ratios_computed.max)] #coordinates in ratios table of two clusters to be merged
		x,y = [y,x] if ratios[[x,y]][0] > ratios[[x,y]][1]
		contaminant, real_haplotype = x,y
		
		#remove from distance table
		dists.keys.select{|key| key.include? contaminant}.each{|key| dists.delete key}
		#remove from ratios table
		deleted_cluster_size = ratios[[contaminant,contaminant]][0]
		merged_cluster_size = ratios[[contaminant,real_haplotype]].reduce(:+)
#		puts "Merging #{[contaminant,real_haplotype]} "
		ratios.keys.select{|key| key.include? contaminant}.each{|key| ratios.delete key}
			
		#recalculate remaining ratios table
		
		ratios.keys.select{|key| key.include? real_haplotype}.each{|key|  #each row and column intersecting with the haplotype being kept
			if key == [real_haplotype,real_haplotype]							#if this is the location corresponding with the haplotype's identity ratio...
				ratios[key] = [merged_cluster_size, merged_cluster_size] 			#assign it the new size
			elsif key[0] == real_haplotype											#otherwise if the real haplotype is the first position in the ratio...
				ratios[key] = [merged_cluster_size, ratios[key][1]]		#increment the first position of the ratio
			elsif key[1] == real_haplotype											#otherwise if the real haplotype is the second position in the ratio...
				ratios[key] = [ratios[key][0], merged_cluster_size]		#increment the denominator of the ratio
			end
		}
		
=begin	
	ratios.keys.select{|a,b| a == b}.each do |remaining_haplo|
		h = haplo_list[remaining_haplo.first]
		puts "#{ratios[remaining_haplo][0]}\t#{haplo_list[remaining_haplo.first].join}"
	end
=end
=begin
		puts "ratios"
		print ratios
		puts
		puts "dists"
		print dists
		puts
		puts '=================='
#		$stdin.gets
=end
	end

	major_constituent = ratios.values.transpose.first.max
	coi = 0
	counter = 0
	
	out = ratios.keys.select{|a,b| a == b}.map{ |remaining_haplo|
		counter += 1
		fraction = (ratios[remaining_haplo][0] / major_constituent.to_f)
		h = haplo_list[remaining_haplo.first]
		coi += 1
		[ratios[remaining_haplo][0], haplo_list[remaining_haplo.first].join]
	}
	
	if verbose
		puts "Clustering removed #{haplos_with_reads.length - out.length} haplotypes" 
	end
	
	return out
end

def generate_fasta(haplo_info, interval, sample)
	return haplo_info.map{|h|
		">#{sample.split('/').last}__#{interval.first}:#{interval.second.first}-#{interval.second.last}__#{h.first}" + "\n" + h.second
	}
end

def split_cigar(cig)
	split = cig.scan /\p{Alpha}+|\p{Digit}+/u
	split = split.map{|thing|thing.scan(/\p{Digit}/).empty? ? thing : thing.to_i}
	grouped = [split.values_at(* split.each_index.select {|i| i.even?}),
	split.values_at(* split.each_index.select {|i| i.odd?})].transpose
	return grouped
end


def parse_locs_table(filepath, idx)
	locs = dice_file(filepath)
	if locs[0].length == 5 #'answer sheet' format with expected per-region haplotype count in fifth column
		locs = locs.map{|l| [l[0], (l[1].to_i..l[2].to_i).to_a, l[4].to_i]}
	else #normal BED
		locs = locs.map{|l| [l[0], (l[1].to_i..l[2].to_i).to_a]}
	end
	
	if idx.include? ','
		start = ((idx.split(',')[0].to_i - 1) * (locs.length / idx.split(',')[1].to_i))
		finish = idx.split(',')[0].to_i * (locs.length / idx.split(',')[1].to_i)
		return locs[start..finish]
	end
		
	if idx == '-'
		return locs
	else
		if idx.include? '-'
			return locs[(idx.split('-')[0].to_i - 1)..(idx.split('-')[1].to_i - 1)]
		else
			return [locs[idx.to_i - 1]]
		end
	end
end

def get_reads(sample, location, verbose, min_mq, min_read_occurrence)
	
	chromosome = location[0]
	coords = location[1]
	reads = dice_text(`samtools view -F 4 #{sample} #{chromosome}:#{coords.min}-#{coords.max}`).select{|r| r[4].to_i > min_mq}.map{|r| [r[5], r[9], r[3].to_i]}
	puts "#{reads.length} reads retrieved" if verbose
	reads_downsampled = reads.uniq.map{|r| [reads.count(r), r]}
	initial_count = reads_downsampled.length
	reads_downsampled.reject!{|r| r.first <= min_read_occurrence}
	#just use reads spanning the interval
#	puts print reads_downsampled.first; abort
	frequency_count = reads_downsampled.length
	reads_downsampled = reads_downsampled.select{|r| r.second.last < coords.min and r.second.last + r.second.second.length > coords.max}
	return nil if reads_downsampled == nil
	reads_downsampled.sort!{|a,b| b.first <=> a.first}
	puts "#{initial_count} distinct haplotypes retrieved." if verbose
	puts "#{initial_count - frequency_count} singleton reads removed." if verbose and min_read_occurrence > 0
	puts "#{frequency_count - reads_downsampled.length} non-spanning haplotypes removed." if verbose
	return reads_downsampled
end

def n_filter(reads, verbose)
	initial_length = reads.length
	reads.reject!{|r| r.second.second.include? 'N' unless r.second.second.split('').last == 'N' and r.second.second.split('').count('N') == 1}
	puts "N-Filtered #{initial_length - reads.length} haplotypes" if verbose
end

def assemble_haplotypes(loci, min_occurrence, subset_size)
	return [] if loci.length == 0
	haplos = []
	#assemble haplotypes
	loci[loci.keys[0]].keys.each do |read_name| #each read...
		haplotype = loci.keys.map{ |key| loci[key][read_name]}
		haplos << haplotype if not haplotype.include? nil
	end
	
	passing_haplos = haplos.uniq.select{|h| (haplos.count(h) / haplos.length.to_f) > min_occurrence and haplos.count(h) > 1}
	haplos_high_occurrence = haplos.select{|h| passing_haplos.include? h}

	return haplos_high_occurrence
end

def get_aligned_read(seq, split_cigar, interval, mapped_start)
	aligned_read = ''
#	qual_string = ''
	
	split_cigar.each do |cig|
		case cig[1]
			when 'M', 'N', 'I'
				aligned_read += seq[0...cig[0]]
			#	qual_string += qual[0...cig[0]]
				seq = seq[cig[0]...seq.length]
#			when 'I'
#				return nil
#				seq = seq[cig[0]...seq.length]
#				qual_string += qual[0...cig[0]]
#				aligned_read += seq[0...cig[0]]
			when 'D'
				aligned_read += (['-'] * cig[0]).join
			#	qual_string += ([' '] * cig[0]).join
			when 'S', 'H'
				seq = seq[cig[0]...seq.length]
			#	qual_string += qual[0...cig[0]]

		end
	end
#	qual_i = []
#	qual_string.each_byte { |c| qual_i << c.to_i - 33}
	desired_aligned_read_segment = aligned_read[(interval.last.first - mapped_start)..(interval.last.last - mapped_start)]
	return [desired_aligned_read_segment]
end

#loci: {locus1 => {read1 => nuc, read2=> nuc}, locus2 => etc etc}
def get_base(reads, locs, min_quality, trim = TRUE, min_allele_occurrence)
	all_genomic_positions = {}
	interval = [locs.first, [locs.last.min, locs.last.max]]
	aligned_haplotypes = []
	reads.each do |read|
		s = read
		cigar = split_cigar(s.second.first)
		seq = s.second.second
		mapped_start = s.second.third
		aligned_haplotypes += [[read.first, get_aligned_read(seq, cigar, interval, mapped_start)].flatten]
	end
	aligned_haplotypes = aligned_haplotypes.select{|c, h| h.length == locs.last.max - locs.last.min + 1}
	uniq_seqs = aligned_haplotypes.map{|h| h.second}.uniq
	condensed = []
	uniq_seqs.each{|h| condensed += [[aligned_haplotypes.select{|ah| ah.second == h}.map{|ah| ah.first}.sum, h]]}
	return condensed
end

def frequency_threshold(haplos, min_frequency, remove_singletons, verbose)
	initial_count = haplos.length
	total_haplos = haplos.map{|h| h.first}.sum
	haplos.reject!{|h| h.first < min_frequency * total_haplos}	
	
	singleton_count = haplos.length
	haplos.reject!{|h| h.first == 1} if remove_singletons
	
	puts "Frequency filtered #{initial_count - singleton_count} haplotypes (of #{initial_count})" if verbose
	puts "Removed #{singleton_count - haplos.length} singleton haplotypes." if verbose and remove_singletons

end

def hamming_distance(s1, s2)
	dist = 0
	s1.each_index{|pos| dist += 1 if s1[pos] != s2[pos]}
	return dist
end


















=begin
def downsample_reads(reads)
	out_reads = []
	counts = {}
	reads.each do |r|
		seq = r.split("\t")[9]
		if counts[seq] == nil
			counts[seq] = 1
		else
			counts[seq] += 1
		end
		
		if counts[seq] < 200
		end
	end
	return
	
end
=end

=begin
def headers(simple_output)
	if not simple_output
		puts "Chromosome\tStart\tEnd\tLoci\tReads\tCOI"
	else
		puts "Chromosome\tStart\tEnd\tLoci\tReads\t1st/1st\t2nd/1st\t3rd/1st\t..."
	end
end
=end


=begin
def filter_monomorphic_biallelic_loci(loci, minor_allele_freqs, max_freq, biallelic, verbose)
	passing_loci = loci.select{ |locus| (minor_allele_freqs[locus].select{|r| r.nan?}.length == 0 ) and minor_allele_freqs[locus].max <= max_freq and not loci[locus].values.include? '-'}
	polymorphism = passing_loci.length
	
	passing_loci = passing_loci.select{|locus| minor_allele_freqs[locus].select{|r| !r.nan?}.select{|t| t>0}.length == 2} if biallelic
	biallelism = passing_loci.length
	
#	puts "#{passing_loci.length} loci passing monomorphism filter" if verbose
	puts "Polymorphic loci: #{polymorphism}, biallelic loci: #{biallelism}" if verbose
	if polymorphism == 0
		return loci
	else
		return passing_loci
	end
end
=end

=begin
def remove_singleton_haplotypes(haplos)
	return haplos.select{|h| haplos.count(h) > 1}
end

def report_haplo_frequencies(location, haplos, freq_threshold, freq_fraction_threshold, depth_threshold, loci, verbose)
	final_output = []
	nreads = haplos.length
	out = [nreads]
	freqs = {}
	sorted_haplos = haplos.sort{|h1, h2| haplos.count(h2) <=> haplos.count(h1) }
	sorted_haplos.uniq.each do |h|
		if (haplos.count(h)/haplos.length.to_f).round(3) > freq_threshold 
			freqs[h.join] = haplos.count(h)
		end
	end

	filt_freqs = {}
	freqs.keys.each{|k| filt_freqs[k] = freqs[k] if freqs[k] > (freqs.values.max.to_f * freq_fraction_threshold)}
	freqs = {}
	filt_freqs.keys.map{|t| freqs[t] = (filt_freqs[t].to_f / filt_freqs.values.sum).round(3)}
#	filt_freqs.keys.map{|t| freqs[t] = (filt_freqs[t].to_f / filt_freqs.values.max).round(3)}
	
	out += [freqs.values.sort.reverse]
	
	if nreads >= depth_threshold
		final_output << out.flatten.fill(0, out.flatten.length..10).map{|t| t.to_f}
		if not verbose
			puts [location[0], location[1].first.to_s, location[1].last.to_s, loci.length.to_s, nreads.to_s, out.second.length.to_s].join("\t")
		else			
			puts [location[0], location[1].first.to_s, location[1].last.to_s, loci.length.to_s, out.join("\t")].join("\t")
		end

	else
		puts location[0] + "\t" + location[1].first.to_s + "\t" + location[1].last.to_s + "\t" + loci.length.to_s + "\t\tDP < #{depth_threshold}"
	end
	
	return final_output
end

def report_allele_frequencies(loci, location, silent = false)
	mafs = {}
	loci.keys.each do |key| 
		bases = loci[key]
		puts "Nucleotide frequencies at #{location[0]}:#{key}:" unless silent
		puts "Total count: #{bases.length}" unless silent
		
		fractions = ['A', 'C', 'G', 'T', '-'].map {|nuc|
			fraction = (bases.values.count(nuc)/bases.values.length.to_f).round(3)
			puts "#{nuc}: #{fraction}" unless silent
			fraction
		}
		puts '--' unless silent
		
		depth = bases.length
		#[key, fractions, depth]
		mafs[key] = fractions
	end
	return mafs
end
=end




#cruft
=begin
def filter_locus_depth_and_gene(loci, location, gff, min_freq, min_depth, min_sites)
	
	nonuniform_deep = loci.select{|l| nonzero = l[1].select{|foo| foo > 0}; nonzero.min > min_freq and l[2] > min_depth and nonzero.length > 1 }.transpose[0]
	genes = `tabix #{gff} #{location[0]}:#{location[1][0]}-#{location[1][1]}`
	if nonuniform_deep != nil and nonuniform_deep.length > min_sites and not (genes.include? 'stevor' or genes.include? 'rif' or genes.include? 'var')
		return nonuniform_deep
	else
		return nil
	end
end


def is_unique_enough(location, max_bit_score)
	formatted_location = location[0]+':'+location[1].join('-')
	tmp = '/tmp/tmp.fa'
	`samtools faidx ~/9/genome.fasta #{formatted_location} > #{tmp}`
	blast = `blastn -db ~/9/blast/Pfal_9.0.blast -query #{tmp} -dust no -outfmt 7 -max_target_seqs 10 | sed 's/\\#\\ Fields\\:\\ //g' | sed 's/,\\ /\\t/g' |  grep -v \\#`.split("\n").map{|t| t.split("\t")}
	#puts blast.map{|t| t.join("\t")}[0..5]
	blast = blast.select{|b| not (b[1] == location[0] and ((b[8].to_i >= location[1][0] and b[8].to_i <= location[1][1]) or b[9].to_i <= location[1][1] and b[9].to_i >= location[1][0]))}
	blast = blast.sort{|a, b| b.last.to_i <=> a.last.to_i}
	#puts "\n\n\n==="
	#puts blast.map{|t| t.join("\t")}[0..5]
	

	if blast[0].last.to_i > max_bit_score
		return false
	else
		puts blast[0].join("\t")
		return true
	end
end


def filter(reads, silent = false)
	puts "#{reads.length} total reads found."  unless silent
	filtered = reads.select{ |read|
		flag = read.split("\t")[1].to_i#.to_s(2).split('')
		flag_binary = read.split("\t")[1].to_i.to_s(2).reverse.split('').map{|blar| blar.to_i}
		
		
		flag < 256 and flag_binary[2] == 0
	#	flag[0] == 0 and
	#	flag[1] == 1 and
	#	flag[2] == 0 and
	#	flag[1] == 0 and
#		flag[0] == 0
=begin

0 0x1 template having multiple segments in sequencing
1 0x2 each segment properly aligned according to the aligner
2 0x4 segment unmapped
3 0x8 next segment in the template unmapped
4 0x10 SEQ being reverse complemented
5 0x20 SEQ of the next segment in the template being reversed
6 0x40 the rst segment in the template
7 0x80 the last segment in the template
8 0x100 secondary alignment
9 0x200 not passing quality controls
10 0x400 PCR or optical duplicate


		
	}
	puts "#{filtered.length} reads remaining after filtration." unless silent
	return filtered
end

=end
