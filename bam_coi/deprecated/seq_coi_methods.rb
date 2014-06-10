#!/usr/bin/env ruby
def get_base(file, location, min_quality)
	all_bases = []
	location[1].each do |coord|
		pu = `samtools mpileup #{file} -r #{location[0]}:#{coord}-#{coord} -d 1000000`.split("\t")
		bases = pu[4].chomp.split('')
		all_bases += bases
	end
	haplos = all_bases.transpose
	print haplos
	abort
end


def headers(simple_output)
	if not simple_output
		puts "Chromosome\tStart\tEnd\tLoci\tReads\tCOI"
	else
		puts "Chromosome\tStart\tEnd\tLoci\tReads\t1st/1st\t2nd/1st\t3rd/1st\t..."
	end
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

def assemble_haplotypes(loci)
	return [] if loci.length == 0
	haplos = []
	#assemble haplotypes
	loci[loci.keys[0]].keys.each do |read_name| #each read...
		haplotype = loci.keys.map{ |key| loci[key][read_name]}
		haplos << haplotype if not haplotype.include? nil
	end
	
	return haplos
end
 
			
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
	filt_freqs.keys.map{|t| freqs[t] = (filt_freqs[t].to_f / filt_freqs.values.max).round(3)}
	
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
