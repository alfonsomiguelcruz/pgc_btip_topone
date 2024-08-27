# chmod +x nextclade-x86_64-unknown-linux-gnu
nextclade-x86_64-unknown-linux-gnu run inputs/sample_xbc.1/ph_recom_unaligned.fasta \ 
    --input-ref inputs/wuhan_reference_genome.fasta \
    --output-fasta inputs/sample_xbc.1/ph_recom_aligned.fasta

nextclade-x86_64-unknown-linux-gnu run inputs/sample_xbc.1/ph_nonrecom_unaligned.fasta \ 
    --input-ref inputs/wuhan_reference_genome.fasta \
    --output-fasta inputs/sample_xbc.1/ph_nonrecom_aligned.fasta

nextclade-x86_64-unknown-linux-gnu run inputs/sample_xbc.1/ph_mixed_unaligned.fasta \ 
    --input-ref inputs/wuhan_reference_genome.fasta \
    --output-fasta inputs/sample_xbc.1/ph_mixed_aligned.fasta
