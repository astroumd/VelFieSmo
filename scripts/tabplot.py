#! /usr/bin/env python
#
#  Simple table plotter, e.g.
#         tabplot.py *.tab
#  would plot each table in a new color.
#
#  based on:
#  http://stackoverflow.com/questions/5260593/is-there-a-command-line-tool-for-data-visualization-and-analysis
#
#  http://matplotlib.org/faq/usage_faq.html#matplotlib-pyplot-and-pylab-how-are-they-related
#
#  @todo    have a batch mode

import sys
import matplotlib.pyplot as plt

for filename in sys.argv[1:]:
   with open(filename,'rt') as sf:
      table = []
      for line in sf:
         if line[0] == '#': continue
         table.append( [float(val) for val in line.split()] )
      table = [ row for row in table if len(row) ] ## remove empty rows
      if len(table[0]) == 1 :
         plt.plot( [y[0] for y in table ] )
      for x in range(1,len(table[0])):
         plt.plot( [ y[0] for y in table ], [ y[x] for y in table ] )
plt.xlabel("Frequency (GHz)")
plt.ylabel("Flux")
#plt.title()
plt.show()
