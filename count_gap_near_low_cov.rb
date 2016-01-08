#!/bin/env ruby
require 'bio'

cov_r_file = ARGV.shift
fa_file = ARGV.shift
rhash=Hash.new
open(cov_r_file).each_line do |l|
  l =~ /(\w+):(\d+)-(\d+)/
  scaff=$1
  rstart=$2
  rend=$3
  rhash[scaff] = Array.new if rhash[scaff] == nil
  rhash[scaff] << [rstart,rend]
end

Bio::FlatFile.open(nil,fa_file).each do |fe|
  if rhash[fe.entry_id] != nil
    rhash[fe.entry_id].each do |ar|
      rstart = ar[0].to_i-10000
      rend = ar[1].to_i+10000
      region_seq=fe.seq.subseq(rstart,rend)
      if region_seq.composition["N"] < 5000
        puts "#{fe.entry_id}:#{ar[0]}-#{ar[1]} #{region_seq.composition}"
      end
    end
  end
end
