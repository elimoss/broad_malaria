#!/usr/bin/env ruby

#filter out rif/var/stevor hits in a .qassoc file

header = true
ARGF.read.each_line("\n") do |l|
	if header
		header = false
		next
	end
	s = l.split(' ')
	locus = s[1].split("-")[0]+':'+s[1].split("-")[1]+'-'+s[1].split("-")[1]
	gene = `tabix ~/9/genes.gff.gz #{locus}`
	puts l if not (gene.downcase.include? 'var' or gene.downcase.include? 'stevor' or gene.downcase.include? 'rif')
end
