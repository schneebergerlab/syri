# Run nucmer
nucmer --maxmatch -t 10 -c 100 -b 500 -l 50 tair10_chr4.fa.gz ler_chr4.fa.gz       # Whole genome alignment. Any other alignment can also be used.
delta-filter -m -i 90 -l 100 out.delta > out.filtered.delta     # Remove small and lower quality alignments
show-coords -THrd out.filtered.delta > out.filtered.coords      # Convert alignment information to a .TSV format as required by SyRI

# Run minimap2
minimap2 -ax asm5 -t 10 --eqx tair10_chr4.fa.gz ler_chr4.fa.gz > out.sam
samtools view -O BAM out.sam > out.bam


minimap2 -cx asm5 -t 10 --eqx tair10_chr4.fa.gz ler_chr4.fa.gz > out.paf