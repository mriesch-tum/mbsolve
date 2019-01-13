#!/bin/bash

grep -e "Performance: \|Number of threads: " | tr -d '\n' | sed -e 's/Performance: /,/g' | sed -e 's/ MGPU\/s/\n/g' | sed -e 's/Number of threads: //g'
