#!/usr/bin/env ruby

#handy utilities for the rest of the code I write and use.


def dice_file(file_loc)
	return File.open(file_loc).read.split("\n").map{|t| t.split("\t").length > 1 ? t.split("\t") : t}
end

def dice_text(input)
	return input.split("\n").map{|t| t.split("\t")}
end

def but(message)
	abort(message)
end

class Array
	def second
		return self[1]
	end
	
	def third
		return self[2]
	end
	
	def fourth
		return self[3]
	end
		
	def mean
		return self.map{|t| t.to_f}.reduce(:+) / self.length.to_f
	end
	
	def sum
		return self.map{|t| t.to_f}.reduce(:+)
	end
end

def dice_gff(file_loc)
	#return a hash containing hash[gene_id] => [exon1,exon2...]
	raw = File.open(file_loc).read.split("\n")
	raw = raw.select{|t| (not t.include? '#') and t.split("\t").length > 4}
	raw = dice_text raw.join("\n")
	gff = {}
	raw.each do |l|
		id = l[8].split(";").first.split("=").last.match(/P.*_[0123456789]{7}/).to_s
		if l[2] == 'exon'
			if not gff.keys.include? id
				gff[id] = [l]
			else
				gff[id] += [l]
			end
		end
	end
	gff
end

def dice_fasta(file_loc)
	file = File.open(file_loc)
	line = file.gets
	out = {}
	prev_id = ''
	while(line != nil)
		line.chomp!
		if line.include? '>'
			abort("Duplicate sequence identifiers") if out.keys.include? line.delete('>')
			prev_id = line.delete('>')
			out[prev_id] = ''
		else
			out[prev_id] += line
		end
		line = file.gets
	end
	out
end

class String
	def include_with_mismatches?(query, num_mismatches)
		return true if self.include? query
		patterns = [query]
				
		num_mismatches.times do
			patterns = patterns.map{|p| 
				p.split('').each_index.map{|idx| 
					p[0...idx]+'.'+p[(idx+1)..p.length]
				}
			}.flatten.uniq
		end
		
		patterns.each do |p| 
			return true if self.match(p) != nil
		end
		return false
	end
	
	def hamming(query)
		return nil if query.length != self.length
		count = 0
		query.length.times{|idx| count += 1 if query[idx] != self[idx]}
		count
	end
	
	def include_with_mismatches2?(query, max_mismatches, substrings)
		substrings.each do |sub|
			return true if sub.hamming(query) <= max_mismatches
		end
		return false
	end
	
	def all_substrings(k)
		subs = []
		(self.length - k+1).times do |idx|
			str = self.split('')[idx..(idx+k - 1)].join
			subs << str
		end
		subs.uniq!
		return subs
	end
	
	def revcomp()
		self.split('').map{|n| 
			case n
				when 'A','a'
					'T'
				when 'C','c'
					'G'
				when 'G','g'
					'C'
				when 'T','t'
					'A'
			end
		}.join.reverse
	end
end



amino_acids = dice_text("A	R	N	D	C	Q	E	G	H	I	L	K	M	F	P	S	T	W	Y	V").flatten

grantham_table = dice_text \
"0	112	111	126	195	91	107	60	86	94	96	106	84	113	27	99	58	148	112	64
112	0	86	96	180	43	54	125	29	97	102	26	91	97	103	110	71	101	77	96
111	86	0	23	139	46	42	80	68	149	153	94	142	158	91	46	65	174	143	133
126	96	23	0	154	61	45	94	81	168	172	101	160	177	108	65	85	181	160	152
195	180	139	154	0	154	170	159	174	198	198	202	196	205	169	112	149	215	194	192
91	43	46	61	154	0	29	87	24	109	113	53	101	116	76	68	42	130	99	96
107	54	42	45	170	29	0	98	40	134	138	56	126	140	93	80	65	152	122	121
60	125	80	94	159	87	98	0	98	135	138	127	127	153	42	56	59	184	147	109
86	29	68	81	174	24	40	98	0	94	99	32	87	100	77	89	47	115	83	84
94	97	149	168	198	109	134	135	94	0	5	102	10	21	95	142	89	61	33	29
96	102	153	172	198	113	138	138	99	5	0	107	15	22	98	145	92	61	36	32
106	26	94	101	202	53	56	127	32	102	107	0	95	102	103	121	78	110	85	97
84	91	142	160	196	101	126	127	87	10	15	95	0	28	87	135	81	67	36	21
113	97	158	177	205	116	140	153	100	21	22	102	28	0	114	155	103	40	22	50
27	103	91	108	169	76	93	42	77	95	98	103	87	114	0	74	38	147	110	68
99	110	46	65	112	68	80	56	89	142	145	121	135	155	74	0	58	177	144	124
58	71	65	85	149	42	65	59	47	89	92	78	81	103	38	58	0	128	92	69
148	101	174	181	215	130	152	184	115	61	61	110	67	40	147	177	128	0	37	88
112	77	143	160	194	99	122	147	83	33	36	85	36	22	110	144	92	37	0	55
64	96	133	152	192	96	121	109	84	29	32	97	21	50	68	124	69	88	55	0"

@grantham = {}

amino_acids.each_with_index do |aa1, i1|
	amino_acids.each_with_index do |aa2, i2|	
		@grantham[[aa1, aa2]] = grantham_table[i1][i2].to_i
		@grantham[[aa2, aa1]] = grantham_table[i2][i1].to_i
	end
end

def grantham_distance(seq1, seq2)
	abort("Sequence length mismatch") if seq1.length != seq2.length
	seq1.length.times.map{|i| @grantham[[seq1[i], seq2[i]]]}.sum
end
