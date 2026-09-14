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
    parser.add_argument("--output", help="Write config to this file (key: value lines)")
    parser.add_argument("--spike_genome", help="Spike-in genome section name in information_config (e.g., dm6)")
    parser.add_argument("--spike_fasta", help="Direct spike FASTA path (overrides --spike_genome lookup)")
    parser.add_argument("--spike_gtf", help="Direct spike GTF path (optional, reserved)")
    parser.add_argument("--spike_index", help="Pre-built combined (main+spike) bowtie2 index prefix; skips on-the-fly concat+build")
    parser.add_argument("--spike_chroms", help="Spike chromosome name list file (required with --spike_index)")
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

    # -------- 解析 spike-in 参考（可选；无 spike 则输出空串）---------
    # spike_genome 是 information_config 的 section 名（与 --build 同一语义，直接查段）；
    # 不经过 species_config 的 species→section 映射（spike 通常是跨物种参考，如 dm6）。
    spike_fasta = ''
    spike_gtf   = ''
    # 直接给 spike_fasta 路径时优先，完全跳过 spike_genome 段查找；
    # 只有没给 fasta 时才用 spike_genome 去 information_config 查段。
    if args.spike_fasta:
        spike_fasta = args.spike_fasta
    elif args.spike_genome:
        spike_section = None
        if args.spike_genome in info_config:
            spike_section = args.spike_genome
        else:
            for sec in info_config.sections():
                if sec.lower() == args.spike_genome.lower():
                    spike_section = sec
                    break
        if spike_section is None:
            sys.exit(f"ERROR: Spike section '{args.spike_genome}' not found in information_config "
                     f"(set spike_genome to a valid section, or give spike_fasta directly)")
        spike_fasta = info_config[spike_section].get('genome_fasta', '')
        spike_gtf   = info_config[spike_section].get('gtf', '')
    if args.spike_gtf:
        spike_gtf = args.spike_gtf

    # 预构建合并索引（可选）：spike_index 给「主+spike 合并」索引前缀，跳过现建 concat+build；
    # spike_chroms 给 spike 染色体名单文件（配合 spike_index，不再从 spike_fasta 派生）。
    spike_index  = args.spike_index or ''
    spike_chroms = args.spike_chroms or ''
    if spike_index and not spike_chroms:
        sys.exit("ERROR: --spike_index requires --spike_chroms (spike chromosome name list file)")

    # 输出键值对（一行一个）
    lines = [
        f"build: {species_index}",
        f"bowtie2_index: {bowtie2_index}",
        f"gtf: {gtf}",
        f"gene_annotation: {gene_annotation}",
        f"rRNA_index: {rrna_index}",
        f"genome_fasta: {genome_fasta}",
        f"spike_fasta: {spike_fasta}",
        f"spike_gtf: {spike_gtf}",
        f"spike_index: {spike_index}",
        f"spike_chroms: {spike_chroms}"
    ]
    for line in lines:
        print(line)

    # 可选：写入文件，供 pipeline 直接消费
    if args.output:
        with open(args.output, 'w') as f:
            for line in lines:
                f.write(line + '\n')

if __name__ == "__main__":
    main()