#!/bin/bash
sed '/^$/d' $1 > tt.$2
grep -v "[\#:@]" tt.$2 > ttt.$2
mv ttt.$2 $1
rm -f tt.$2
awk '{print $1, $1, $2, "X C", $4, $5, $6, $8, $3, $7}' $1 > uploads/tmp.$2;
