#!/usr/bin/env python3
# parse_config.py

import configparser
import argparse
import sys

def main():
    parser = argparse.ArgumentParser(description="Parse species and information configs")
    parser.add_argument("--species_config", default="/workplace/database/database_species_taxid_config.txt", help="Path to species_config.txt")
    parser.add_argument("--information_config", default="/workplace/database/database_information_config.txt", help="Path to information_config.txt")
    parser.add_argument("--species", default="human_19", help="Species name (e.g., human_19)")
    parser.add_argument("--build", help="Build name (e.g., hg38), overrides species mapping")
    parser.add_argument("--gtf", help="Direct GTF file path, overrides config")
    parser.add_argument("--genome_fasta", help="Direct genome FASTA path (optional)")
    args = parser.parse_args()

    # -------- 读取 species_config ----------
    sp_config = configparser.ConfigParser()
    sp_config.read(args.species_config)
    if 'species' not in sp_config:
        sys.exit("ERROR: [species] section not found in species_config")
    
    # 如果指定了 --build，直接用它作为 section；否则查 species
    if args.build:
        species_index = args.build.upper()
    else:
        species_key = args.species.lower()
        if species_key not in sp_config['species']:
            sys.exit(f"ERROR: Species '{args.species}' not found in species_config")
        species_index = sp_config['species'][species_key].split('\t')[1]  # 取第二个字段
    
    # -------- 读取 information_config ----------
    info_config = configparser.ConfigParser()
    info_config.read(args.information_config)
    if species_index not in info_config:
        sys.exit(f"ERROR: Section '{species_index}' not found in information_config")
    
    # 提取所需键值
    try:
        bowtie2_index = info_config[species_index]['bowtie2_index']
        gtf = info_config[species_index]['gtf']
        gene_annotation = info_config[species_index]['gene_annotation']
        rrna_index = info_config[species_index]['rRNA_index']
    except KeyError as e:
        sys.exit(f"ERROR: Missing key {e} in section '{species_index}'")
    
    # 如果命令行提供了 --gtf，则覆盖
    if args.gtf:
        gtf = args.gtf

    # 输出键值对（一行一个）
    print(f"bowtie2_index: {bowtie2_index}")
    print(f"gtf: {gtf}")
    print(f"gene_annotation: {gene_annotation}")
    print(f"rRNA_index: {rrna_index}")

    # 可选：输出 GENOME_FASTA 如果提供
    if args.genome_fasta:
        print(f"genome_fasta: {args.genome_fasta}")

    # 也可以输出 BUILD 供参考（如果需要）
    print(f"build: {species_index}")

if __name__ == "__main__":
    main()