minimap2 -N 1000 -ax asm20 --eqx -Y --secondary=yes -p 0.5\
  LPA.all.fasta  tmp.lpa.motif.fasta | samIdentity.py --header /dev/stdin
