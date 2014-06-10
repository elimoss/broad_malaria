#!/usr/bin/env ruby
require('/home/unix/emoss/scripts/kit/util.rb')
#input: key relating R- codes to strain names, directory containing fasta files for strains, list of strains in TEP mock mixture R- code
#output: disagreeing polymorphic sites, base pairs occurring at those sites
abort("Usage: clashing_polymorphism_anticipator.rb code-strain_key(or -) clustal_omega_seq_alignment amplicon.ped windows.ped strain1,strain2,strain3,..(or -)") if ARGV.length != 5

key = Hash[*dice_file(ARGV[0]).drop(1).flatten] if not ARGV[0] == '-' #create nice and tidy hash with code keys and strain values 
aln = ARGV[1]
amplicon = dice_file(ARGV[2])[0]
windows = dice_file ARGV[3]
strains = ARGV[4].split(',')

#parse the clustal omega output into a fasta_label => sequence hash
t = File.open(aln).read
t.gsub!('  ', ' ') while t.include? '  '
t = dice_text(t.gsub(' ', "\t")).select{|r| r.length > 0}
seqs = {}
t.drop(1).each do |chunk| 
	if chunk.first != ''
		if seqs[chunk.first] == nil
			seqs[chunk.first] = chunk.last
		else
			seqs[chunk.first] += chunk.last
		end
	end
end

if strains[0] == '-'
	strains = seqs.keys
	key = Hash[*strains.map{|s| [s,s]}.flatten]
end

#find sites that aren't monomorphic in this collection as well as haplotypes

puts "Chromosome\tStart\tEnd\t\tSites\tHaplos"
#puts "Position\t#{strains.join("\t")}"

#amplicon
windows.each do |window|
	haplotypes = []
	polymorphic_sites = 0
	((window[1].to_i - amplicon[1].to_i)..(window[2].to_i - amplicon[1].to_i)).each do |pos|
		site_polymorphism = strains.map{|s| seqs[seqs.keys.find{|e| Regexp.new(key[s]) =~ e}][pos]}
		polymorphic_sites += 1 if site_polymorphism.uniq.length != 1
		haplotypes << site_polymorphism
	end
#	haplotypes.transpose.each {|h| print h.join(); puts}
	haplotype_count = haplotypes.transpose.map{|t| t[1..t.length].join}.uniq.length
	puts window.join("\t") +"\t"+ polymorphic_sites.to_s + "\t" +haplotype_count.to_s + "\t"# + haplotypes.last.join("\t") 
end

