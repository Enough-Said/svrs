#!/usr/bin/env bash

FILES=$(cat clean/splitData.txt)

for file in $FILES; do
    cat $file.zip.?? > "$file.zip"
    rm $file.zip.??
    gunzip -S .zip "$file.zip"
done


