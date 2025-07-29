import pysam
import pandas as pd
import argparse
import os

def main():
    parser = argparse.ArgumentParser(description="Extract spanning reads for multiple regions from a BAM file.")
    parser.add_argument("-b", "--bam", required=True, help="Input BAM file (must be indexed)")
    parser.add_argument("-r", "--regions", required=True, help="Tab file with columns: name, chr, start, stop")
    parser.add_argument("-o", "--outdir", required=True, help="Directory to store output BAM files")

    args = parser.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    # Derive sample name from BAM filename (without path or extension)
    bam_filename = os.path.basename(args.bam)
    sample_prefix = os.path.splitext(bam_filename)[0]

    # Load input regions
    df = pd.read_csv(args.regions, sep="\t")
    required_cols = {'name', 'chr', 'start', 'stop'}
    if not required_cols.issubset(df.columns):
        raise ValueError(f"Input region file must contain columns: {required_cols}")

    # Open BAM
    bam = pysam.AlignmentFile(args.bam, "rb")

    results = []

    for _, row in df.iterrows():
        print(row)
        region_name = row['name']
        chrom = row['chr']
        start = int(row['start'])
        stop = int(row['stop'])

        output_bam_path = os.path.join(args.outdir, f"{sample_prefix}.{region_name}.bam")
        out_bam = pysam.AlignmentFile(output_bam_path, "wb", template=bam)

        count = 0
        for read in bam.fetch(chrom, start, stop):
            if read.is_unmapped:
                continue
            if read.reference_start <= start and read.reference_end >= stop:
                out_bam.write(read)
                count += 1

        out_bam.close()

        # Index the output BAM file
        pysam.index(output_bam_path) 
               
        results.append({
            "name": region_name,
            "chr": chrom,
            "start": start,
            "stop": stop,
            "read_count": count,
            "sample": sample_prefix
        })

    bam.close()

    

    # Write summary table
    output_table_path = os.path.join(args.outdir, f"{sample_prefix}.counts.txt")
    summary_df = pd.DataFrame(results)
    summary_df.to_csv(output_table_path, sep="\t", index=False)
    print(f"✅ Summary written to {output_table_path}")

if __name__ == "__main__":
    main()
