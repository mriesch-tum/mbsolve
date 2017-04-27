#!/bin/bash

grep -e "Performance: \|Thread count: " | tr -d '\n' | sed -e 's/Performance: /,/g' | sed -e 's/ MGPU\/s/\n/g' | sed -e 's/Thread count: //g'

# | gnuplot performance.gnuplot
