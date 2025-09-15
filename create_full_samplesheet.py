#!/usr/bin/env python3

import os
import csv

def create_full_samplesheet():
    assemblies_dir = "/fastscratch/tylerdoe/ks_samples_project/assemblies"
    output_file = "full_kansas_samplesheet.csv"
    
    samples = []
    
    # Map directory names to organism names for AMRFinder+
    organism_map = {
        'campylobacter': 'Campylobacter',
        'salmonella': 'Salmonella', 
        'ecoli': 'Escherichia'
    }
    
    for organism_dir in ['campylobacter', 'salmonella', 'ecoli']:
        full_path = os.path.join(assemblies_dir, organism_dir)
        if os.path.exists(full_path):
            for sample_dir in os.listdir(full_path):
                scaffolds_file = os.path.join(full_path, sample_dir, "scaffolds.fasta")
                if os.path.exists(scaffolds_file):
                    # Extract sample ID from directory name (remove _spades suffix)
                    sample_id = sample_dir.replace("_spades", "")
                    samples.append({
                        'sample': sample_id,
                        'fasta': scaffolds_file,
                        'organism': organism_map[organism_dir]
                    })
    
    # Write the samplesheet
    with open(output_file, 'w', newline='') as csvfile:
        fieldnames = ['sample', 'fasta', 'organism']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        for sample in samples:
            writer.writerow(sample)
    
    print(f"Created {output_file} with {len(samples)} samples:")
    for organism in organism_map.values():
        count = len([s for s in samples if s['organism'] == organism])
        print(f"  {organism}: {count} samples")

if __name__ == "__main__":
    create_full_samplesheet()
