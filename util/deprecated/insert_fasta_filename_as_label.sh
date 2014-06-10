#!/usr/bin/env ruby

require('/home/unix/emoss/scripts/kit/util.rb')

fasta = dice_file(ARGV[0])
fasta[0] = ">#{ARGV[0]}"

File.open(ARGV[0], 'w').write(fasta[0] + "\n" + fasta[1..fasta.length].join + "\n")
