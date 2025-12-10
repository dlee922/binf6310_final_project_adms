#!/bin/bash
#SBATCH --job-name=align_tree
#SBATCH --output=logs/align_tree_%j.out
#SBATCH --error=logs/align_tree_%j.err
#SBATCH --partition=short
#SBATCH --mem=32G
#SBATCH --cpus-per-task=16
#SBATCH --time=04:00:00

module load miniconda3
source activate ecoli_final

cd /scratch/evangelatos.s/final_project/core_genes_v2

echo "****************************************************"
echo "Step 1: Aligning core genes with MAFFT"
echo "Started at: $(date)"
echo "****************************************************"

CORE_GENES=("rpoB" "gyrB" "recA" "dnaK" "tuf" "fusA")

# Align each gene
for gene in "${CORE_GENES[@]}"; do
    echo ""
    echo "=== Aligning $gene ==="
    
    # Check if gene directory exists and has files
    if [ ! -d "$gene" ] || [ -z "$(ls -A ${gene})" ]; then
        echo "  WARNING: No sequences found for $gene, skipping"
        continue
    fi
    
    # Concatenate all sequences for this gene
    cat ${gene}/*.faa > ${gene}_all.faa
    
    # Count sequences
    n_seqs=$(grep -c "^>" ${gene}_all.faa)
    echo "  Sequences to align: $n_seqs"
    
    if [ $n_seqs -lt 10 ]; then
        echo "  WARNING: Too few sequences for $gene, skipping"
        continue
    fi
    
    # Align with MAFFT
    echo "  Running MAFFT..."
    mafft --auto --thread 16 ${gene}_all.faa > ${gene}_aligned.faa
    
    echo "  ✓ $gene alignment complete"
done

echo ""
echo "****************************************************"
echo "Step 2: Concatenating alignments"
echo "Started at: $(date)"
echo "****************************************************"

# Python script to concatenate alignments properly
python3 << 'PYTHON_SCRIPT'
from collections import defaultdict

genes = ["rpoB", "gyrB", "recA", "dnaK", "tuf", "fusA"]
alignments = defaultdict(dict)

# Read all aligned sequences
for gene in genes:
    filename = f"{gene}_aligned.faa"
    try:
        with open(filename, 'r') as f:
            current_sample = None
            current_seq = []
            
            for line in f:
                line = line.strip()
                if line.startswith('>'):
                    # Save previous sequence
                    if current_sample:
                        alignments[current_sample][gene] = ''.join(current_seq)
                    # Start new sequence
                    current_sample = line[1:]  # Remove '>'
                    current_seq = []
                else:
                    current_seq.append(line)
            
            # Save last sequence
            if current_sample:
                alignments[current_sample][gene] = ''.join(current_seq)
        
        print(f"✓ Loaded {gene}: {len(alignments)} samples")
    except FileNotFoundError:
        print(f"✗ Warning: {filename} not found, skipping {gene}")

# Write concatenated alignment in FASTA format
print(f"\nConcatenating alignments for {len(alignments)} samples...")

with open("concatenated_alignment.fasta", "w") as out:
    for sample in sorted(alignments.keys()):
        # Concatenate all genes in order
        concat_seq = ""
        for gene in genes:
            if gene in alignments[sample]:
                concat_seq += alignments[sample][gene]
        
        # Write in FASTA format
        out.write(f">{sample}\n")
        # Write sequence in lines of 60 characters
        for i in range(0, len(concat_seq), 60):
            out.write(concat_seq[i:i+60] + "\n")

print(f"✓ Concatenated alignment written: {len(alignments)} samples")
print(f"  Alignment length: {len(concat_seq)} amino acids")
PYTHON_SCRIPT

echo ""
echo "****************************************************"
echo "Alignment and concatenation complete at: $(date)"
echo "****************************************************"

# Check output
echo ""
echo "=== Output Files ==="
ls -lh *_aligned.faa
ls -lh concatenated_alignment.fasta

echo ""
echo "=== Alignment Statistics ==="
echo "Number of samples:"
grep -c "^>" concatenated_alignment.fasta

echo "Alignment length (first sequence):"
grep -v "^>" concatenated_alignment.fasta | head -1 | wc -c
