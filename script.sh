#!/bin/bash

declare -a NS=($(seq 2 2 10))

num=100
file="../networks"
mean=0.5
std=1000

for value in "${NS[@]}"
do 
  ./solve_gridgraph "$value" "$num" "$file" "$mean" "$std"
done
