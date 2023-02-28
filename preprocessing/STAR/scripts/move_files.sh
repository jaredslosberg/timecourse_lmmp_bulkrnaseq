#!/bin/bash

cd /mnt/morbo/Data/Users/jslosberg/aged_bulkrna/preprocessing/STAR

while IFS= read -r line; do
    echo "Text read from file: $line"
    mkdir -p $line
    for filename in ./*"$line"*; do
         echo $filename
         
       [ ! -f "$filename" ] && continue
       mv -i "$filename" "$line"/
   done
done < ../../samples_merged.txt