import sys

# Takes as input the list of filtered read group header data

with open(sys.argv[1]) as f1:
    for line in f1:
        ld = line.strip().split()
        sample = ld[0].split('.bam')[0].split('/')[1]
        data = ld[2:]
        rg = {}
        for i in data:
            split_data = i.split(':')
            field = split_data[0]
            value = ':'.join(split_data[1:])
            rg[field] = value
        command = "module load picard/3.1.1; picard MarkDuplicates -I 02_realign/" + sample + ".realignGA5_C4masked.sort.bam -O /dev/stdout -M 03_markdups_addRG/" + sample + ".metrics.txt | picard AddOrReplaceReadGroups -I /dev/stdin -O 03_markdups_addRG/" + sample + ".realignGA5_C4masked.sort.picard.bam --RGID " + rg["ID"] + " --RGPL " + rg["PL"] + " --RGPU " + rg["PU"] + " --RGLB " + rg["LB"] + " --RGCN " + rg["CN"] + " --RGSM '" + rg["SM"] + "' --VALIDATION_STRINGENCY LENIENT"
	print 'sbatch --account=kingsley --job-name=mkdup_addRGs --mem=24G --time=24:00:00 --wrap="' + command + '"'
