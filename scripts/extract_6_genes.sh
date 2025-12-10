#!/bin/bash
#SBATCH --job-name=extract_6genes
#SBATCH --output=logs/extract_6genes_%j.out
#SBATCH --error=logs/extract_6genes_%j.err
#SBATCH --partition=short
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4
#SBATCH --time=01:00:00

module load miniconda3
source activate ecoli_final

cd /scratch/evangelatos.s/final_project

# Create fresh output directory
mkdir -p core_genes_v2
cd core_genes_v2

echo "****************************************************"
echo "Extracting 6 core genes from Prokka annotations"
echo "Started at: $(date)"
echo "****************************************************"

# All 6 core genes
CORE_GENES=("rpoB" "gyrB" "recA" "dnaK" "tuf" "fusA")

for gene in "${CORE_GENES[@]}"; do
    echo ""
    echo "=== Extracting $gene ==="
    mkdir -p ${gene}
    
    found=0
    
    for annot_dir in ../annotations/*/; do
        sample=$(basename $annot_dir)
        tsv_file="${annot_dir}${sample}.tsv"
        faa_file="${annot_dir}${sample}.faa"
        
        # Find locus tag for this gene
        locus_tag=$(awk -F'\t' -v gene="$gene" '$4 ~ gene {print $1; exit}' "$tsv_file")
        
        if [ ! -z "$locus_tag" ]; then
            # Extract sequence using locus tag
            grep -A 1 "^>$locus_tag " "$faa_file" | \
                sed "s/>.*/>$sample/" > ${gene}/${sample}.faa
            
            found=$((found + 1))
        fi
    done
    
    echo "Found $gene in $found/50 genomes"
done

echo ""
echo "****************************************************"
echo "Extraction complete at: $(date)"
echo "****************************************************"

# Summary
echo ""
echo "=== FINAL SUMMARY ==="
for gene in "${CORE_GENES[@]}"; do
    n_files=$(ls ${gene}/*.faa 2>/dev/null | wc -l)
    echo "$gene: $n_files genomes"
done

echo ""
echo "Total sequences: $(find . -name '*.faa' | wc -l)"
