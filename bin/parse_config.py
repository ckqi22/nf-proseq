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
    
    # 提取所需键值（gtf 必需，其余可选：缺失输出空串，避免部分 section 缺键导致退出）
    gtf = info_config[species_index].get('gtf', '')
    if not gtf and not args.gtf:
        sys.exit(f"ERROR: Missing key 'gtf' in section '{species_index}'")

    bowtie2_index  = info_config[species_index].get('bowtie2_index', '')
    gene_annotation = info_config[species_index].get('gene_annotation', '')
    rrna_index     = info_config[species_index].get('rRNA_index', '')
    genome_fasta   = info_config[species_index].get('genome_fasta', '')

    # 命令行覆盖（用户指定路径优先于数据库）
    if args.gtf:
        gtf = args.gtf
    if args.genome_fasta:
        genome_fasta = args.genome_fasta

    # 输出键值对（一行一个）
    print(f"bowtie2_index: {bowtie2_index}")
    print(f"gtf: {gtf}")
    print(f"gene_annotation: {gene_annotation}")
    print(f"rRNA_index: {rrna_index}")
    print(f"genome_fasta: {genome_fasta}")

    # 也可以输出 BUILD 供参考（如果需要）
    print(f"build: {species_index}")

    # TODO(spike-in): 未来在此解析 spike 拼接（主基因组 + spike fasta/gtf 合并），
    # 输出 combined fasta/gtf 路径，供 prepare_genome.nf 的 CONCAT 扩展点使用。

if __name__ == "__main__":
    main()