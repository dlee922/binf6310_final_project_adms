#!/bin/bash
#SBATCH --job-name=download_ecoli
#SBATCH --output=logs/download_%j.out
#SBATCH --error=logs/download_%j.err
#SBATCH --partition=short
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=8G
#SBATCH --time=06:00:00

# Load conda
module load miniconda
source activate ecoli_project

# Navigate to project
cd /scratch/evangelatos.s/final_project

# Create assemblies directory
mkdir -p assemblies
cd assemblies

echo "Starting downloads at $(date)"
echo "Total accessions: $(wc -l < ../FINAL_accession_list_v3.txt)"
echo ""

# Download all assemblies
# The CP format means RefSeq - need to download each contig separately
count=0

while read accession_range; do
    count=$((count + 1))
    echo "=== Processing $count/50: $accession_range ==="
    
    # Parse range (e.g., CP163825-CP163828)
    start=$(echo $accession_range | cut -d'-' -f1)
    end=$(echo $accession_range | cut -d'-' -f2)
    
    # Extract numeric parts
    start_num=$(echo $start | sed 's/CP//')
    end_num=$(echo $end | sed 's/CP//')
    
    # Create a directory for this isolate
    isolate_dir="${accession_range}"
    mkdir -p "$isolate_dir"
    
    # Download each contig in the range
    for num in $(seq $start_num $end_num); do
        accession="CP${num}"
        echo "  Downloading $accession..."
        
        # Use NCBI efetch API
        curl -s "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=${accession}&rettype=fasta&retmode=text" \
            -o "${isolate_dir}/${accession}.fna"
        
        # Check file size
        size=$(stat -f%z "${isolate_dir}/${accession}.fna" 2>/dev/null || stat -c%s "${isolate_dir}/${accession}.fna")
        if [ $size -lt 1000 ]; then
            echo "    WARNING: ${accession} download may have failed (${size} bytes)"
        else
            size_mb=$(echo "scale=2; $size/1000000" | bc)
            echo "    Success: ${accession} (${size_mb} Mb)"
        fi
        
        sleep 0.5  # Be nice to NCBI
    done
    
    # Concatenate all contigs for this isolate into one file
    cat ${isolate_dir}/*.fna > ${accession_range}.fna
    echo "  Combined into ${accession_range}.fna"
    
done < ../FINAL_accession_list_v3.txt

echo ""
echo "=== Download Summary ==="
echo "Completed at: $(date)"
echo "Total .fna files created: $(ls *.fna 2>/dev/null | wc -l)"
echo "Total size: $(du -sh . | cut -f1)"
echo ""

# List first 10 files with sizes
echo "First 10 genome files:"
ls -lh *.fna | head -10
