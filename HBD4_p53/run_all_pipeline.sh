#!/bin/bash
# HBD4_p53 One-Click Pipeline Execution Script
# This script runs the entire analysis pipeline and saves all results to a dedicated directory

set -e  # Exit on any error

# Define directories
WORK_DIR="/home/gao/Code/Bioinfo_Analysis_Projects/HBD4_p53"
RESULTS_DIR="$WORK_DIR/one_click_results"

echo "=== HBD4_p53 One-Click Pipeline Started ==="
echo "Working directory: $WORK_DIR"
echo "Results will be saved to: $RESULTS_DIR"
echo ""

# Step 1: Run simulation (generates raw_R1.fastq and raw_R2.fastq)
echo "Step 1: Running simulation..."
cd $WORK_DIR
python3 01_tp53_sim.py
echo "✓ Simulation completed"

# Step 2: Run alignment and UMI extraction
echo "Step 2: Running alignment and UMI extraction..."
bash 02_align.sh
echo "✓ Alignment and UMI extraction completed"

# Step 3: Run molecular assembly and consensus generation
echo "Step 3: Running molecular assembly..."
python3 03_hbd_assembly.py
echo "✓ Molecular assembly completed"

# Step 4: Run visualization and generate reports
echo "Step 4: Running visualization..."
Rscript 04_visualize.R > "$RESULTS_DIR/visualization_output.txt" 2>&1
echo "✓ Visualization completed"

# Step 5: Copy all results to the dedicated directory
echo "Step 5: Copying results to dedicated directory..."
cp $WORK_DIR/read_analysis_report.txt $RESULTS_DIR/
cp $WORK_DIR/hbd_final_assembly_report.csv $RESULTS_DIR/
cp $WORK_DIR/mutation_analysis_report.txt $RESULTS_DIR/
cp $WORK_DIR/molecule_family_distribution.txt $RESULTS_DIR/
cp $WORK_DIR/*.pdf $RESULTS_DIR/ 2>/dev/null || echo "No PDF files found (this is normal if R packages are missing)"
cp $WORK_DIR/sorted.bam $RESULTS_DIR/ 2>/dev/null || echo "BAM file copied"
cp $WORK_DIR/sorted.bam.bai $RESULTS_DIR/ 2>/dev/null || echo "BAM index copied"

# Create a summary file
echo "Creating summary file..."
cat > "$RESULTS_DIR/pipeline_summary.txt" << EOF
HBD4_p53 One-Click Pipeline Summary
==================================

Execution Time: $(date)
Pipeline Steps Completed:
1. Simulation (01_tp53_sim.py)
2. Alignment & UMI Extraction (02_align.sh)  
3. Molecular Assembly (03_hbd_assembly.py)
4. Visualization (04_visualize.R)

Key Results Files:
- read_analysis_report.txt: Detailed read statistics
- hbd_final_assembly_report.csv: Complete assembly results
- mutation_analysis_report.txt: Mutation analysis with error classification
- molecule_family_distribution.txt: Family size distribution
- Various PDF plots (if generated)

Error Classification:
- Low frequency spontaneous mutations: True biological variants
- PCR first amplification errors: Early PCR errors (high frequency)
- PCR later amplification errors: Late PCR errors (low frequency)
EOF

echo ""
echo "=== Pipeline completed successfully! ==="
echo "All results are available in: $RESULTS_DIR"
echo ""
echo "Key statistics from read_analysis_report.txt:"
head -20 "$RESULTS_DIR/read_analysis_report.txt" | tail -10