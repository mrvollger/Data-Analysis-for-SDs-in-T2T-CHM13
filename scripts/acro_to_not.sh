#!/usr/bin/env bash
set -euo pipefail




grep -P "chr(13|14|15|21|22)\s" \
  ../Assembly_analysis/SEDEF/chm13.draft_v1.0_plus38Y.SDs.bed | \
  grep -P "chr(1|2|3|4|5|6|7|8|9|10|11|12|16|17|18|19|20|Y|X)\s" \
  > tmp.sd.acro.bed
bedtools merge -i tmp.sd.acro.bed > tmp.sd.acro.merged.bed


~/projects/assembly_workflows/scripts/enriched.py --scale 7 -S 3000001 --drop -3000000 \
  tmp.sd.acro.merged.bed \
  ../Assembly_analysis/SEDEF/chm13.draft_v1.0_plus38Y_masked.fasta.fai | \
  sed 's/ /\t/g' > ./large.enriched.bed




#grep -P "^chr(1|2|3|4|5|6|7|8|9|10|11|12|16|17|18|19|20|Y|X)\s" tmp.sd.acro.bed

