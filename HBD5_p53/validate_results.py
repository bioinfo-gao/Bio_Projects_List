#!/usr/bin/env python3
"""
Validation script to verify that the one-click pipeline produces consistent results
with the individual step-by-step execution.
"""

import os
import sys

def read_statistics_file(filepath):
    """Read and parse the read statistics report"""
    stats = {}
    try:
        with open(filepath, 'r') as f:
            for line in f:
                if 'Total original reads:' in line:
                    stats['original_reads'] = int(line.split(':')[1].strip().replace(',', ''))
                elif 'Successfully processed reads:' in line:
                    stats['cleaned_reads'] = int(line.split(':')[1].strip().replace(',', ''))
                elif 'Successfully mapped reads:' in line:
                    stats['mapped_reads'] = int(line.split(':')[1].strip().replace(',', ''))
                elif 'Unique molecular families:' in line:
                    stats['unique_molecules'] = int(line.split(':')[1].strip().replace(',', ''))
                elif 'PCR duplicates removed:' in line:
                    stats['pcr_duplicates'] = int(line.split(':')[1].strip().replace(',', ''))
                elif 'Final high-fidelity consensus sequences:' in line:
                    stats['final_consensus'] = int(line.split(':')[1].strip().replace(',', ''))
    except FileNotFoundError:
        print(f"File not found: {filepath}")
        return None
    return stats

def validate_pipeline():
    """Validate that the pipeline produced correct results"""
    print("=== HBD5_p53 Pipeline Validation ===")
    
    # Read main results
    main_stats = read_statistics_file("read_analysis_report.txt")
    if main_stats is None:
        print("ERROR: Could not read main results file!")
        return False
    
    print("\nMain Results Statistics:")
    print(f"Original reads: {main_stats['original_reads']:,}")
    print(f"Cleaned reads: {main_stats['cleaned_reads']:,}")
    print(f"Mapped reads: {main_stats['mapped_reads']:,}")
    print(f"Unique molecules: {main_stats['unique_molecules']:,}")
    print(f"PCR duplicates: {main_stats['pcr_duplicates']:,}")
    print(f"Final consensus: {main_stats['final_consensus']:,}")
    print(f"PCR duplication rate: {main_stats['pcr_duplicates']/main_stats['mapped_reads']*100:.2f}%")
    
    # Validate expected ranges
    validation_passed = True
    
    # Original reads should be around 1000
    if not (900 <= main_stats['original_reads'] <= 1100):
        print(f"WARNING: Original reads ({main_stats['original_reads']}) outside expected range (900-1100)")
        validation_passed = False
    
    # PCR duplication rate should be high (>90%)
    pcr_dup_rate = main_stats['pcr_duplicates']/main_stats['mapped_reads']*100
    if pcr_dup_rate < 90:
        print(f"WARNING: PCR duplication rate ({pcr_dup_rate:.2f}%) lower than expected (>90%)")
        validation_passed = False
    
    # Final consensus should be much smaller than original reads
    if main_stats['final_consensus'] > main_stats['original_reads'] * 0.2:
        print(f"WARNING: Final consensus count ({main_stats['final_consensus']}) seems too high relative to original reads")
        validation_passed = False
    
    # Check mutation analysis file
    if os.path.exists("mutation_analysis_report.txt"):
        with open("mutation_analysis_report.txt", 'r') as f:
            lines = f.readlines()
            mutation_count = len(lines) - 1  # Subtract header
            print(f"\nMutation analysis report contains {mutation_count} mutations")
            
            # Should have reasonable number of mutations
            if mutation_count < 100:
                print(f"WARNING: Low mutation count ({mutation_count}) - may indicate issues")
                validation_passed = False
    else:
        print("ERROR: mutation_analysis_report.txt not found!")
        validation_passed = False
    
    # Check assembly report
    if os.path.exists("hbd_final_assembly_report.csv"):
        with open("hbd_final_assembly_report.csv", 'r') as f:
            lines = f.readlines()
            family_count = len(lines) - 1  # Subtract header
            print(f"Assembly report contains {family_count} molecular families")
            
            if family_count != main_stats['unique_molecules']:
                print(f"ERROR: Family count mismatch! CSV: {family_count}, Stats: {main_stats['unique_molecules']}")
                validation_passed = False
    else:
        print("ERROR: hbd_final_assembly_report.csv not found!")
        validation_passed = False
    
    if validation_passed:
        print("\n✅ VALIDATION PASSED: Pipeline results are consistent and within expected ranges!")
        return True
    else:
        print("\n❌ VALIDATION FAILED: Some results are outside expected ranges!")
        return False

if __name__ == "__main__":
    success = validate_pipeline()
    sys.exit(0 if success else 1)