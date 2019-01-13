#!/bin/bash

grep -e "Performance:" | tr -d '\n' | sed -e 's/Performance: //g' | sed -e 's/ MGPU\/s/, /g'
echo ""
