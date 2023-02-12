#!/bin/bash

# Generate figures
echo "Generate framome stats..."
Rscript framome_stats.R > /dev/null 2>&1

echo "Generate in silico immunogenic NOPs..."
Rscript in_silico_immunogenic_NOPs.R > /dev/null 2>&1

echo "Generate lrRNA sequencing stats..."
Rscript lrRNA_seq_stats.R > /dev/null 2>&1

echo "Generate lung binding assays..."
Rscript lung_binding_assays.R > /dev/null 2>&1

echo "Generate NOP expression..."
Rscript nop_expression.R > /dev/null 2>&1

echo "Generate sequencing stats..."
Rscript seq_stats.R > /dev/null 2>&1

echo "Generate Ribo-seq signaling..."
Rscript riboseq_signaling.R > /dev/null 2>&1

echo "Done"
