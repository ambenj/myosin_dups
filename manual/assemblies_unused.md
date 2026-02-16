# Long read assembly analysis (unused)

This section describes methods used to attempt assembly of PacBio CLR and Oxford Nanopore datasets. Due to lower read accuracy, we were unable to achieve confident high quality assemblies for these samples. Local assemblies were also attempted in an effort to simplify assembly but were also unsuccessful. The read depths of GIFU and and LITC samples from Ishikawa et al. were too low for assembly.

## Gynogenetic genome (*GasAcuGyn*)
### Check of available genome assembly
Run liftoff to transfer annotations gynogen assembly:
```bash
conda activate myosin_annotate
sbatch scripts/assemblies/liftoff_NCBI.sh /labs/kingsley/ambenj/myosin_dups/analysis/assemblies/Gynogen/Gynogen_pchrom_assembly_all.fasta
```

Align Gynogen pacbio reads to self genome assembly:
```bash
conda activate myosin_assembly
sbatch scripts/assemblies/pbmm2_align.sh \
/labs/kingsley/ambenj/myosin_dups/analysis/assemblies/Gynogen/Gynogen_pchrom_assembly_all.fasta \
/labs/kingsley/ambenj/myosin_dups/raw/sra/gynogen_pacbio.fofn \
/labs/kingsley/ambenj/myosin_dups/analysis/assemblies/Gynogen/Gynogen_pchrom_assembly_all.bam \
SUBREAD
```
### Local *MYH3C* region assembly
Try mapping reads to BLAU and performing local *MYH3C* region assembly:
```bash
# Map gynogenetic reads to BLAU
conda activate myosin_assembly
cd /labs/kingsley/ambenj/myosin_dups

# Index BLAU genome with subread setting (hap1 + hap2 although just hap1 would be sufficient)
sbatch scripts/assemblies/pbmm2_index.sh /labs/kingsley/ambenj/myosin_dups/analysis/assemblies/BLAU22_12/hifiasm/BLAU22_12.bp.hap12.p_ctg.fa SUBREAD

# Align Gynogen pacbio reads to BLAU genome
sbatch scripts/assemblies/pbmm2_align.sh \
/labs/kingsley/ambenj/myosin_dups/analysis/assemblies/BLAU22_12/hifiasm/BLAU22_12.bp.hap12.p_ctg.SUBREAD.mmi \
/labs/kingsley/ambenj/myosin_dups/raw/sra/gynogen_pacbio.fofn \
/labs/kingsley/ambenj/myosin_dups/analysis/assemblies/Gynogen/Gynogen.BLAU22_12.mapped.bam \
SUBREAD

# Extract reads from MYH region
sbatch scripts/assemblies/extract_reads_ROI.sh h1tg000036l:2640861-3140861 \
/labs/kingsley/ambenj/myosin_dups/analysis/assemblies/Gynogen/Gynogen.BLAU22_12.mapped.bam \
/labs/kingsley/ambenj/myosin_dups/analysis/assemblies/Gynogen/Gynogen.BLAU22_12.mapped.chrXIX_myh500kb.fa

# Start de novo assembly
sbatch scripts/assemblies/canu_pacbio.sh /labs/kingsley/ambenj/myosin_dups/analysis/assemblies/Gynogen/Gynogen.BLAU22_12.mapped.chrXIX_myh500kb.fa Gynogen_BLAU22_12.mapped.chrXIX_myh500kb 500k
```
## Lake Storvatnet ONT libraries
### *De novo* assembly
```bash
conda activate assembly

# De novo assemble all Lake Storvatnet, Norway spineless reads
sbatch scripts/assemblies/canu_nanopore.sh /labs/kingsley/ambenj/myosin_dups/raw/sra/SRR31134576.fastq /labs/kingsley/ambenj/myosin_dups/analysis/assemblies/Storvatnet_spineless/canu Storvatnet_spineless 500m

# De novo assemble all Lake Storvatnet, Norway spined reads
sbatch scripts/assemblies/canu_nanopore.sh /labs/kingsley/ambenj/myosin_dups/raw/sra/SRR31134577.fastq /labs/kingsley/ambenj/myosin_dups/analysis/assemblies/Storvatnet_spined/canu Storvatnet_spined 500m

# BLAT assemblies:
sbatch scripts/assemblies/blat.sh /labs/kingsley/ambenj/myosin_dups/analysis/assemblies/Storvatnet_spined/canu/Storvatnet_spined/Storvatnet_spined.contigs.fasta

sbatch scripts/assemblies/blat.sh /labs/kingsley/ambenj/myosin_dups/analysis/assemblies/Storvatnet_spineless/canu/Storvatnet_spineless/Storvatnet_spineless.contigs.fasta
```

### Local MYH3C region assembly
Try mapping reads to BLAU and performing local *MYH3C* region assembly:
```bash
# Map reads to BLAU
sbatch scripts/assemblies/minimap.sh /labs/kingsley/ambenj/myosin_dups/analysis/assemblies/BLAU22_12/hifiasm/BLAU22_12.bp.hap12.p_ctg.fa /labs/kingsley/ambenj/myosin_dups/raw/sra/SRR31134576.fastq /labs/kingsley/ambenj/myosin_dups/analysis/assemblies/Storvatnet_spineless Storvatnet_spineless.BLAUmapped

sbatch scripts/assemblies/minimap.sh /labs/kingsley/ambenj/myosin_dups/analysis/assemblies/BLAU22_12/hifiasm/BLAU22_12.bp.hap12.p_ctg.fa /labs/kingsley/ambenj/myosin_dups/raw/sra/SRR31134577.fastq /labs/kingsley/ambenj/myosin_dups/analysis/assemblies/Storvatnet_spined Storvatnet_spined.BLAUmapped

# Extract reads from MYH region
sbatch scripts/assemblies/extract_reads_ROI.sh h1tg000036l:2640861-3140861 \
/labs/kingsley/ambenj/myosin_dups/analysis/assemblies/Storvatnet_spineless/Storvatnet_spineless.BLAUmapped.sorted.bam \
/labs/kingsley/ambenj/myosin_dups/analysis/assemblies/Storvatnet_spineless/Storvatnet_spineless.BLAUmapped.chrXIX_myh500kb.fa

sbatch scripts/assemblies/extract_reads_ROI.sh h1tg000036l:2640861-3140861 \
/labs/kingsley/ambenj/myosin_dups/analysis/assemblies/Storvatnet_spined/Storvatnet_spined.BLAUmapped.sorted.bam \
/labs/kingsley/ambenj/myosin_dups/analysis/assemblies/Storvatnet_spined/Storvatnet_spined.BLAUmapped.chrXIX_myh500kb.fa

# Start local de novo assembly for spineless sample
conda activate assembly
sbatch scripts/assemblies/canu_nanopore.sh /labs/kingsley/ambenj/myosin_dups/analysis/assemblies/Storvatnet_spineless/Storvatnet_spineless.BLAUmapped.chrXIX_myh500kb.fa /labs/kingsley/ambenj/myosin_dups/analysis/assemblies/Storvatnet_spineless/canu Storvatnet_spineless.BLAUmapped.chrXIX_myh500kb 500k

# Start local de novo assembly for spined sample
sbatch scripts/assemblies/canu_nanopore.sh /labs/kingsley/ambenj/myosin_dups/analysis/assemblies/Storvatnet_spined/Storvatnet_spined.BLAUmapped.chrXIX_myh500kb.fa /labs/kingsley/ambenj/myosin_dups/analysis/assemblies/Storvatnet_spined/canu Storvatnet_spined.BLAUmapped.chrXIX_myh500kb 500k
```

## Paxton Benthic PacBio CLR
### *De novo* assembly
```bash
# De novo assemble all PacBio reads for PAXB
sbatch scripts/assemblies/canu_pacbio.sh /labs/kingsley/ambenj/myosin_dups/analysis/assemblies/PAXB/SRR10527428.fastq PAXB 500m

# BLAT assembly
sbatch scripts/assemblies/blat.sh /labs/kingsley/ambenj/myosin_dups/analysis/assemblies/PAXB/canu_PAXB/PAXB.contigs.fasta
```

### Local *MYH3C* region assembly
Try mapping reads to stickleback v.5 and performing local *MYH3C* assembly
```bash
# Index freshwater reference genome with subread setting to use for multiple alignments
sbatch scripts/assemblies/pbmm2_index.sh /labs/kingsley/ambenj/ref_genomes/Nath2020_sticklebackv5/stickleback_v5_assembly.fa SUBREAD

# Align PAXB pacbio reads to stickleback_v5 assembly
sbatch scripts/assemblies/pbmm2_align.sh \
/labs/kingsley/ambenj/ref_genomes/Nath2020_sticklebackv5/stickleback_v5_assembly.SUBREAD.mmi \
/labs/kingsley/ambenj/myosin_dups/raw/sra/SRR10527428.fastq \
/labs/kingsley/ambenj/myosin_dups/analysis/assemblies/PAXB/PAXB.stickleback_v5mapped.bam \
SUBREAD

# Extract PAXB reads from 500kb region around MYH on chrXIX and chrY
sbatch scripts/assemblies/extract_reads_ROI.sh chrXIX:2416650-2916650 \
/labs/kingsley/ambenj/myosin_dups/analysis/assemblies/PAXB/PAXB.stickleback_v5mapped.bam \
/labs/kingsley/ambenj/myosin_dups/analysis/assemblies/PAXB/PAXB.stickleback_v5mapped_chrXIX_myh500kb.fa

# Start PAXB de novo assemblies of 500kb region around myh on chrXIX and chrY
for c in chrXIX chrY; do sbatch scripts/assemblies/canu_pacbio.sh /labs/kingsley/ambenj/myosin_dups/analysis/assemblies/PAXB/PAXB.stickleback_v5mapped_${c}_myh500kb.fa PAXB_${c}_myh500kb 500k; done
```
