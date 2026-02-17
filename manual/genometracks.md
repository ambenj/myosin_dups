# Genome tracks figure

*Methods for lifting over annotations and visualizing genome tracks*

**Relevant Figure**: Figure 1A

## Liftover annotations to *stickleback v.5* reference
Liftover EcoPeaks from *gasAcu1-4* to *stickleback v.5*:
```
sbatch scripts/genome_tracks/liftover.sh
```

## Plot genome tracks
```bash
conda activate pyGenomeTracks
pyGenomeTracks --tracks scripts/genome_tracks/MYH_tracks.ini --region chrXIX:2600500-2840000 --outFileName analysis/assemblies/stickleback_v5/pyGenomeTracks/MYH_tracks.pdf --trackLabelFraction 0 --plotWidth 16
```
*Note: figure was manually modified after to add fragments of C2 and SYT19*