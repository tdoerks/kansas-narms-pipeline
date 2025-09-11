# Kansas NARMS Bacterial Genomics Pipeline

## Overview
This pipeline extracts and analyzes Kansas-specific bacterial isolates from the NARMS (National Antimicrobial Resistance Monitoring System) surveillance data. The workflow covers Campylobacter, Salmonella, and E. coli samples.

## Prerequisites
- Access to NCBI SRA Run Selector
- Python 3 with pandas library
- Understanding of NARMS sample naming format: `YYStateMMTTNN-S` where State = "KS" for Kansas
- HPC cluster access (or equivalent computational resources)

---

## Step 1: Metadata Download and Parsing

### Overview
Download metadata files from NCBI SRA for NARMS bacterial surveillance data and extract Kansas-specific samples based on standardized sample naming conventions.

### 1.1 Download Metadata from NCBI

Navigate to the NCBI SRA Run Selector and download metadata for each pathogen:

**For Campylobacter:**
1. Go to https://www.ncbi.nlm.nih.gov/sra
2. Search for BioProject: `PRJNA292664`
3. Click "Send results to Run selector"
4. Download metadata as CSV: `SraRunTableCampylobacterPRJNA292664.csv`

**For Salmonella:**
1. Search for BioProject: `PRJNA292661`
2. Download metadata as CSV: `SraRunTableSalmonellaPRJNA292661.csv`

**For E. coli:**
1. Search for BioProject: `PRJNA292663`
2. Download metadata as CSV: `SraRunTableEscherichiacoliPRJNA292663.csv`

### 1.2 Create Project Directory Structure

```bash
# Create main project directory
mkdir /path/to/ks_samples_project
cd /path/to/ks_samples_project

# Create subdirectories for organization
mkdir campylobacter_ks salmonella_ks ecoli_ks
mkdir metadata raw_data assemblies
```

### 1.3 Parse Kansas Samples from Metadata

Create a Python script to extract Kansas samples based on sample naming convention:

```python
# Create file: metadata/parse_kansas_samples.py
#!/usr/bin/env python3

import pandas as pd
import sys

def find_kansas_samples(csv_file, output_prefix):
    # Read the CSV file
    print(f"Reading {csv_file}...")
    df = pd.read_csv(csv_file)
    print(f"Total samples: {len(df)}")
    
    # Check available sample name columns
    sample_name_cols = [col for col in df.columns if 'sample' in col.lower() and 'name' in col.lower()]
    print(f"Sample name columns found: {sample_name_cols}")
    
    # Look for KS pattern in sample names (format: YYStateMMTTNN-S)
    # Where State would be "KS" for Kansas
    kansas_mask = pd.Series([False] * len(df))
    
    for col in sample_name_cols:
        if col in df.columns:
            # Look for pattern like 19KS07CB01-C (YY+KS+MM+rest)
            kansas_mask |= df[col].astype(str).str.contains(r'\d{2}KS\d{2}', case=False, na=False)
    
    # Also check other potential sample columns
    for col in ['Sample Name', 'Sample_name', 'isolate', 'strain']:
        if col in df.columns:
            kansas_mask |= df[col].astype(str).str.contains(r'\d{2}KS\d{2}', case=False, na=False)
    
    # Filter Kansas samples
    kansas_df = df[kansas_mask]
    print(f"Kansas samples found: {len(kansas_df)}")
    
    if len(kansas_df) > 0:
        # Save Kansas samples
        kansas_df.to_csv(f"{output_prefix}_kansas_samples.csv", index=False)
        
        # Extract just SRR accessions for easy downloading
        srr_list = kansas_df['Run'].tolist()
        with open(f"{output_prefix}_kansas_srr_list.txt", 'w') as f:
            for srr in srr_list:
                f.write(srr + '\n')
        
        print(f"Saved {len(kansas_df)} Kansas samples to:")
        print(f"  - {output_prefix}_kansas_samples.csv")
        print(f"  - {output_prefix}_kansas_srr_list.txt")
        
        # Show sample names to verify the pattern
        print("\nSample Kansas sample names:")
        for col in sample_name_cols:
            if col in kansas_df.columns:
                sample_names = kansas_df[col].dropna().unique()[:10]
                print(f"{col}: {sample_names}")
                
        # Parse the sample names to show the structure
        print("\nParsing sample name structure:")
        for col in sample_name_cols:
            if col in kansas_df.columns:
                for name in kansas_df[col].dropna().unique()[:5]:
                    if pd.notna(name) and 'KS' in str(name):
                        print(f"  {name} -> Year: 20{str(name)[:2]}, State: {str(name)[2:4]}, Month: {str(name)[4:6]}")
                        
    else:
        print("No Kansas samples found!")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python3 parse_kansas_samples.py <csv_file> <output_prefix>")
        sys.exit(1)
    
    csv_file = sys.argv[1]
    output_prefix = sys.argv[2]
    find_kansas_samples(csv_file, output_prefix)
```

### 1.4 Run Kansas Sample Extraction

```bash
# Navigate to metadata directory
cd metadata

# Parse each pathogen's metadata for Kansas samples
python3 parse_kansas_samples.py SraRunTableCampylobacterPRJNA292664.csv campylobacter
python3 parse_kansas_samples.py SraRunTableSalmonellaPRJNA292661.csv salmonella
python3 parse_kansas_samples.py SraRunTableEscherichiacoliPRJNA292663.csv ecoli

# Check results
echo "=== KANSAS SAMPLES SUMMARY ==="
echo "Campylobacter: $(wc -l < campylobacter_kansas_srr_list.txt) samples"
echo "Salmonella: $(wc -l < salmonella_kansas_srr_list.txt) samples" 
echo "E. coli: $(wc -l < ecoli_kansas_srr_list.txt) samples"
```

### Expected Results
- Campylobacter: ~70 Kansas samples
- Salmonella: ~206 Kansas samples  
- E. coli: ~710 Kansas samples
- Total: ~986 Kansas samples

### Output Files
- `{pathogen}_kansas_samples.csv` - Full metadata for Kansas samples
- `{pathogen}_kansas_srr_list.txt` - SRR accession numbers for download

---

## Step 2: Sequence Data Download

### Overview
Download raw sequencing data from NCBI SRA using the Kansas sample accession lists generated in Step 1.

### 2.1 Load SRA Toolkit

```bash
# Check if SRA toolkit is available
module avail SRA-Toolkit

# Load the SRA toolkit module
module load SRA-Toolkit/3.0.3-gompi-2022a

# Verify tools are available
which prefetch
which fastq-dump
```

### 2.2 Create Download Scripts

Create download scripts for each pathogen that will download SRA files and convert them to paired FASTQ format:

```bash
# Create download script for Campylobacter
cat > download_kansas_campylobacter.sh << 'EOF'
#!/bin/bash
#SBATCH --job-name=download_campy_ks
#SBATCH --partition=batch.q
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --time=12:00:00
#SBATCH --mem=16G
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=your.email@institution.edu

module load SRA-Toolkit/3.0.3-gompi-2022a

cd /path/to/ks_samples_project/campylobacter_ks

echo "Downloading Kansas Campylobacter samples..."
while read srr; do
    echo "Downloading $srr"
    prefetch $srr
    fastq-dump --split-files --gzip $srr
    rm -rf $srr  # Remove SRA file after conversion to save space
done < ../metadata/campylobacter_kansas_srr_list.txt

echo "Download complete!"
EOF

# Create download script for Salmonella  
cat > download_kansas_salmonella.sh << 'EOF'
#!/bin/bash
#SBATCH --job-name=download_salm_ks
#SBATCH --partition=batch.q
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --time=12:00:00
#SBATCH --mem=16G
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=your.email@institution.edu

module load SRA-Toolkit/3.0.3-gompi-2022a

cd /path/to/ks_samples_project/salmonella_ks

echo "Downloading Kansas Salmonella samples..."
while read srr; do
    echo "Downloading $srr"
    prefetch $srr
    fastq-dump --split-files --gzip $srr
    rm -rf $srr  # Remove SRA file after conversion to save space
done < ../metadata/salmonella_kansas_srr_list.txt

echo "Download complete!"
EOF

# Create download script for E. coli
cat > download_kansas_ecoli.sh << 'EOF'
#!/bin/bash
#SBATCH --job-name=download_ecoli_ks
#SBATCH --partition=batch.q
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --time=12:00:00
#SBATCH --mem=16G
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=your.email@institution.edu

module load SRA-Toolkit/3.0.3-gompi-2022a

cd /path/to/ks_samples_project/ecoli_ks

echo "Downloading Kansas E. coli samples..."
while read srr; do
    echo "Downloading $srr"
    prefetch $srr
    fastq-dump --split-files --gzip $srr
    rm -rf $srr  # Remove SRA file after conversion to save space
done < ../metadata/ecoli_kansas_srr_list.txt

echo "Download complete!"
EOF

# Make scripts executable
chmod +x download_kansas_*.sh
```

### 2.3 Submit Download Jobs

```bash
# Navigate to metadata directory where scripts are located
cd metadata

# Submit all download jobs
sbatch download_kansas_campylobacter.sh
sbatch download_kansas_salmonella.sh
sbatch download_kansas_ecoli.sh

# Monitor job status
squeue -u your_username
```

### 2.4 Monitor Download Progress

```bash
# Check download progress for each pathogen
echo "Campylobacter downloads:"
ls campylobacter_ks/*_1.fastq.gz | wc -l

echo "Salmonella downloads:"  
ls salmonella_ks/*_1.fastq.gz | wc -l

echo "E. coli downloads:"
ls ecoli_ks/*_1.fastq.gz | wc -l
```

### Expected Output
- Each sample produces paired-end FASTQ files: `SRR#######_1.fastq.gz` and `SRR#######_2.fastq.gz`
- Downloads typically take 6-12 hours depending on data size and network speed
- File sizes range from ~100MB to 1GB per sample pair

### Data Organization
After downloads complete, your directory structure should look like:
```
ks_samples_project/
├── campylobacter_ks/
│   ├── SRR#######_1.fastq.gz
│   ├── SRR#######_2.fastq.gz
│   └── ...
├── salmonella_ks/
│   ├── SRR#######_1.fastq.gz  
│   ├── SRR#######_2.fastq.gz
│   └── ...
├── ecoli_ks/
│   ├── SRR#######_1.fastq.gz
│   ├── SRR#######_2.fastq.gz
│   └── ...
└── metadata/
    ├── *_kansas_samples.csv
    └── *_kansas_srr_list.txt
```

---

## Step 3: Quality Control and Assembly

### Overview
Assemble bacterial genomes from paired-end FASTQ files using SPAdes assembler. This step processes completed downloads and produces draft genome assemblies for downstream analysis.

### 3.1 Check Available Software

```bash
# Check if SPAdes is available
module avail SPAdes

# Load SPAdes module
module load SPAdes

# Verify installation
spades.py --version
```

### 3.2 Create Assembly Directory Structure

```bash
# Create assembly directories
mkdir -p assemblies/{campylobacter,salmonella,ecoli}
```

### 3.3 Count Completed Downloads

Before starting assemblies, verify which samples have completed downloading:

```bash
# Count completed paired samples for each pathogen
echo "Counting completed paired samples:"
echo "Campylobacter: $(ls campylobacter_ks/*_1.fastq.gz | wc -l) samples"
echo "Salmonella: $(ls salmonella_ks/*_1.fastq.gz | wc -l) samples"  
echo "E. coli: $(ls ecoli_ks/*_1.fastq.gz | wc -l) samples"
```

### 3.4 Create Assembly Script

Create a comprehensive assembly script that processes all completed samples:

```bash
cat > assemble_kansas_samples.sh << 'EOF'
#!/bin/bash
#SBATCH --job-name=kansas_assembly
#SBATCH --partition=batch.q
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --time=48:00:00
#SBATCH --mem=60G
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=your.email@institution.edu

module load SPAdes

cd /path/to/ks_samples_project

# Function to assemble a sample
assemble_sample() {
    local sample_id=$1
    local input_dir=$2
    local output_dir=$3
    
    echo "$(date): Starting assembly of $sample_id"
    
    spades.py \
        -1 ${input_dir}/${sample_id}_1.fastq.gz \
        -2 ${input_dir}/${sample_id}_2.fastq.gz \
        -o ${output_dir}/${sample_id}_spades \
        --threads 4 \
        --memory 12 \
        --careful \
        --only-assembler
        
    echo "$(date): Completed $sample_id"
}

# Process all three pathogens
for pathogen in campylobacter salmonella ecoli; do
    echo "$(date): Starting ${pathogen} assemblies..."
    
    for file in ${pathogen}_ks/*_1.fastq.gz; do
        if [ -f "$file" ]; then
            sample_id=$(basename "$file" "_1.fastq.gz")
            # Check if R2 exists and assembly doesn't already exist
            if [ -f "${pathogen}_ks/${sample_id}_2.fastq.gz" ] && [ ! -d "assemblies/${pathogen}/${sample_id}_spades" ]; then
                assemble_sample $sample_id ${pathogen}_ks assemblies/${pathogen} &
                # Limit to 4 concurrent jobs (16 cores / 4 threads each)
                (($(jobs -r | wc -l) >= 4)) && wait
            fi
        fi
    done
    
    wait  # Wait for current pathogen to complete
    echo "$(date): Completed all ${pathogen} assemblies"
done

echo "$(date): All assemblies completed!"
EOF

chmod +x assemble_kansas_samples.sh
```

### 3.5 Submit Assembly Job

```bash
# Submit the assembly job
sbatch assemble_kansas_samples.sh

# Monitor job status
squeue -u your_username

# Monitor assembly progress
tail -f slurm-JOBID.out
```

### 3.6 Alternative: Dedicated Node Assembly

If you have access to a dedicated compute node, you can use higher resource allocation:

```bash
cat > assemble_kansas_dedicated.sh << 'EOF'
#!/bin/bash
#SBATCH --job-name=kansas_assembly_dedicated
#SBATCH --nodelist=your_dedicated_node
#SBATCH --ntasks-per-node=32
#SBATCH --time=48:00:00
#SBATCH --mem=120G
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=your.email@institution.edu

module load SPAdes

cd /path/to/ks_samples_project

# Same assembly function but with more parallel jobs
assemble_sample() {
    local sample_id=$1
    local input_dir=$2
    local output_dir=$3
    
    echo "$(date): Starting assembly of $sample_id"
    
    spades.py \
        -1 ${input_dir}/${sample_id}_1.fastq.gz \
        -2 ${input_dir}/${sample_id}_2.fastq.gz \
        -o ${output_dir}/${sample_id}_spades \
        --threads 4 \
        --memory 12 \
        --careful \
        --only-assembler
        
    echo "$(date): Completed $sample_id"
}

# Process with higher parallelism (8 concurrent jobs)
for pathogen in campylobacter salmonella ecoli; do
    echo "$(date): Starting ${pathogen} assemblies..."
    
    for file in ${pathogen}_ks/*_1.fastq.gz; do
        if [ -f "$file" ]; then
            sample_id=$(basename "$file" "_1.fastq.gz")
            if [ -f "${pathogen}_ks/${sample_id}_2.fastq.gz" ] && [ ! -d "assemblies/${pathogen}/${sample_id}_spades" ]; then
                assemble_sample $sample_id ${pathogen}_ks assemblies/${pathogen} &
                (($(jobs -r | wc -l) >= 8)) && wait
            fi
        fi
    done
    
    wait
    echo "$(date): Completed all ${pathogen} assemblies"
done

echo "$(date): All assemblies completed!"
EOF
```

### 3.7 Monitor Assembly Progress

```bash
# Check assembly completion status
echo "Assembly progress:"
echo "Campylobacter: $(ls -d assemblies/campylobacter/*/scaffolds.fasta 2>/dev/null | wc -l) completed"
echo "Salmonella: $(ls -d assemblies/salmonella/*/scaffolds.fasta 2>/dev/null | wc -l) completed"  
echo "E. coli: $(ls -d assemblies/ecoli/*/scaffolds.fasta 2>/dev/null | wc -l) completed"

# Check assembly quality (example for one sample)
head -5 assemblies/campylobacter/SRR*/scaffolds.fasta
```

### Expected Results
- Each assembly produces a directory: `assemblies/{pathogen}/{sample}_spades/`
- Key output files:
  - `scaffolds.fasta` - Final assembled genome
  - `contigs.fasta` - Assembled contigs before scaffolding
  - `assembly_graph.fastg` - Assembly graph
  - `spades.log` - Assembly log file
- Assembly time: 15-45 minutes per sample for bacterial genomes
- Typical assembly stats:
  - Contigs: 20-200 per genome
  - N50: 50-500kb depending on genome quality
  - Total length: ~2-6Mb depending on species

### Quality Assessment (Optional)

Quick assembly quality check using basic statistics:

```bash
# Create a simple assembly stats script
cat > assembly_stats.sh << 'EOF'
#!/bin/bash

echo "Sample,Contigs,Total_Length,N50_estimate"
for fasta in assemblies/*/*/scaffolds.fasta; do
    if [ -f "$fasta" ]; then
        sample=$(basename $(dirname "$fasta"))
        contigs=$(grep -c ">" "$fasta")
        total_len=$(grep -v ">" "$fasta" | tr -d '\n' | wc -c)
        echo "$sample,$contigs,$total_len,TBD"
    fi
done
EOF

chmod +x assembly_stats.sh
./assembly_stats.sh > assembly_summary.csv
```

---

## Next Steps
The following steps will be added as the pipeline progresses:
- Step 4: AMR Analysis
- Step 5: Phylogenetic Analysis
