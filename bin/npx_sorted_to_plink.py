import os, sys, argparse, gzip

parser = argparse.ArgumentParser(prog=f"npx_sorted_to_wide", description="Convert sorted NPX files to wide formats")
parser.add_argument('--npx', type=str, required=True, help='TSV file containing NPX values in long format')
parser.add_argument('--samples', type=str, required=True, help='TSV file containing sample information')
parser.add_argument('--assays', type=str, required=True, help='TSV file containing assay information')
parser.add_argument('--out', type=str, required=True, help='Output file prefix (.npx.wide.tsv.gz, .missing.long.tsv.gz)')
parser.add_argument('--skip-mismatching-samples', action='store_true', help='Skip mismatched samples')

args = parser.parse_args()

def flexopen(filename, mode='rt'):
    if filename.endswith('.gz'):
        return gzip.open(filename, mode, encoding='utf-8' if 't' in mode else None)
    else:
        return open(filename, mode, encoding='utf-8' if 't' in mode else None)
    
samples = []
sample2idx = {}
assays = []
assay2idx = {}

print("Loading sample information")
with flexopen(args.samples) as f:
    for line in f:
        toks = line.strip().split("\t")
        (sampleid, wellid, plateid) = toks[0:3]
        if sampleid in sample2idx:
            print(f"Duplicate sample ID: {sampleid}")
            sys.exit(1)
        sample2idx[sampleid] = len(samples)
        samples.append([sampleid, wellid, plateid])

print("Loading assay information")
with flexopen(args.assays) as f:
    for line in f:
        (olinkid, uniprot, genename, other) = line.strip().split("\t")
        if olinkid in assay2idx:
            print(f"Duplicate assay ID: {olinkid}")
            sys.exit(1)
        assay2idx[olinkid] = len(assays)
        assays.append([olinkid, uniprot, genename])

nmiss = 0
nwrote = 0
(prev_sample, prev_assay) = ("", "")
with gzip.open(args.out + ".missing.long.tsv.gz", "wt") as wmiss:
    with gzip.open(args.out + ".npx.wide.tsv.gz", "wt") as wwide:
        ## write header
        wwide.write("FID\tIID")
        for assay in assays:
            wwide.write("\t" + assay[0] + ":" + assay[1] + ":" + assay[2])
        wwide.write("\n")
        with flexopen(args.npx) as f:
            isample = -1
            iassay = 0
            for line in f:
                (sampleid, olinkid, npx) = line.strip().split("\t")
                if args.skip_mismatching_samples and sampleid not in sample2idx:
                    continue
                if ( prev_sample > sampleid ) or ( prev_sample == sampleid and prev_assay > olinkid ):
                    print(f"Input file is not sorted: {prev_sample} {prev_assay} > {sampleid} {olinkid}")
                    sys.exit(1)
                if sampleid != samples[isample][0]:
                    newisample = sample2idx[sampleid]
                    if newisample != isample + 1:
                        print(f"Missing {newisample-isample-1} samples in a row... expecting {samples[newisample][0]} but observed {sampleid}. Please double check")
                        sys.exit(1)
                    if isample >= 0:
                        while iassay < len(assays):
                            #if sampleid == "TOP120262" and assays[iassay][0] == "OID45299":
                            #    print(f"Writing NA for {prev_sample} - olinkid={olinkid}")
                            wwide.write("\tNA")
                            wmiss.write(f"{prev_sample}\t{assays[iassay][0]}\n")
                            iassay += 1
                            nmiss += 1
                        wwide.write("\n")
                    wwide.write(sampleid + "\t" + sampleid)
                    isample = newisample
                    iassay = 0     
                if olinkid == assays[iassay][0]: ## expected assay ID is observed
                    wwide.write("\t" + npx)
                    iassay += 1
                    nwrote += 1
                else:
                    newiassay = assay2idx[olinkid]
                    while iassay < newiassay:
                        #if sampleid == "TOP120262" and assays[iassay][0] == "OID45299":
                        #    print(f"**Writing NA for {assays[iassay][0]} - olinkid={olinkid}")
                        wwide.write("\tNA")
                        wmiss.write(f"{sampleid}\t{assays[iassay][0]}\n")
                        iassay += 1
                        nmiss += 1
                    wwide.write("\t" + npx)
                    iassay += 1
                    nwrote += 1
                
                (prev_sample, prev_assay) = (sampleid, olinkid)
                if nwrote % 1000000 == 0:
                    print(f"Processed {nwrote} NPX values, {nmiss} ({nmiss/nwrote*100:.3f}%) missing values, at {sampleid} {olinkid}", file=sys.stderr)
            while iassay < len(assays):
                wwide.write("\tNA")
                wmiss.write(f"{sampleid}\t{assays[iassay][0]}\n")
                iassay += 1
                nmiss += 1
            wwide.write("\n")
print(f"Processed {nwrote} NPX values, {nmiss} missing values", file=sys.stderr)
