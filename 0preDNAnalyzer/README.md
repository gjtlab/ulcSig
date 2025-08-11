# preDNAnalyzer

A Snakemake workflow for DNA-seq data preprocessing.

## Overview

This workflow performs the following steps:
1. Quality control and trimming (fastp)
2. Read alignment (BWA-MEM2)
3. BAM processing (sambamba)

## Requirements

- Snakemake ≥7.0
- Conda/Mamba

## Usage

### Basic usage

```bash
# Create and activate conda environment
conda env create -f environment.yml
conda activate preDNAnalyzer

# Run the workflow
snakemake -c N  # Replace N with number of available cores
```

### Important Notes

1. Thread Usage
   - Always specify the number of cores using `snakemake -c N`
   - If not specified, all tools will run with single thread, which may significantly impact performance
   - Example thread allocation:
     * With `-c 32`:
       - BWA-MEM2: up to 16 threads
       - Sambamba merge: up to 8 threads
       - Fastp: up to 4 threads
     * Without `-c`: all tools use 1 thread

2. Input Configuration
   - Configure your input data in `configs/config.yaml`
   - Specify reference genome and other resources

## License

This project is licensed under the MIT License - see the LICENSE file for details.