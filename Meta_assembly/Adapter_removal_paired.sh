awk -F'\t' '
$1=="79RR" || $1=="YG15"||$1=="YG01"||$1=="80R" {
    print "AdapterRemoval --file1 " $3 \
          " --file2 " $4 \
          " --basename " $2 \
          " --threads 4" \
          " --minadapteroverlap 1" \
          " --adapter1 AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC" \
          " --adapter2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT" \
          " --minlength 30 --gzip --trimns --trimqualities"
}' ../config/samples.tsv > adapterremoval_cmds.txt
