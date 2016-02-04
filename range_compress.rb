#!/bin/env ruby
#
lastref=nil
lastpos=0
min_depth=100000
cur_range_start = nil

ARGF.each_line do |l|
  a=l.chomp.split
  ref=a[0]
  pos=a[1].to_i
  depth=a[2].to_i
  min_depth = depth if depth < min_depth
  if lastpos == 0
    cur_range_start = pos
    lastref = ref
    lastpos = pos
    next
  end
  if ref!=lastref or pos != lastpos + 1
    puts "#{lastref}:#{cur_range_start}-#{lastpos}\t#{min_depth}"
    lastref=ref
    cur_range_start=pos
    min_depth=100000
  end
  lastpos=pos
end

