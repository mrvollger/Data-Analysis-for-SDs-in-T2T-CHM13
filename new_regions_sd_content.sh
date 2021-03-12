
for f in ../Assembly_analysis/SEDEF/*.SDs.bed ; do
    echo $f
    awk '$19 >= 10000 && $24 >= 0.95{print $0}' $f \
        | grep -P  "^chr(\d+|X|Y)\s" \
        | bedtools merge -i - |  awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print "10k: ",SUM/10^6}'
    awk '$19 >= 20000 && $24 >= 0.95{print $0}' $f \
        | grep -P "^chr(\d+|X|Y)\s" \
        | bedtools merge -i - |  awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print "20k: ",SUM/10^6 }'
    echo
done

