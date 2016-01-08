#!/bin/env ruby
#
lastref=nil
lastpos=0
cur_range_start = nil

ARGF.each_line do |l|
  a=l.chomp.split
  ref=a[0]
  pos=a[1].to_i
  if lastpos == 0
    cur_range_start = pos
    lastref = ref
    lastpos = pos
    next
  end
  if ref!=lastref or pos != lastpos + 1
    puts "#{lastref}:#{cur_range_start}-#{lastpos}"
    lastref=ref
    cur_range_start=pos
  end
  lastpos=pos
end

