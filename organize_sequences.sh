# Create the FASTA file with unaligned recombinant sequences
cp inputs/sample_xbc.1/ph_xbc.1.fasta inputs/sample_xbc.1/ph_recom_unaligned.fasta

# Create the FASTA file with unaligned nonrecombinant sequences
cat inputs/sample_xbc.1/ph_b.1.617.2.fasta inputs/sample_xbc.1/ph_ba.2.fasta > inputs/sample_xbc.1/ph_nonrecom_unaligned.fasta

# Create the FASTA file with unaligned mixed (recombinant and nonrecombinant) sequences
cat inputs/sample_xbc.1/ph_nonrecom_unaligned.fasta inputs/sample_xbc.1/ph_recom_unaligned.fasta > inputs/sample_xbc.1/ph_mixed_unaligned.fasta