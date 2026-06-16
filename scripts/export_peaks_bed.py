import gzip
from liftover import ChainFile
import os

print("1. Loading liftover mm8 -> mm9 chain file...")
converter = ChainFile('../data/references/mm8ToMm9.over.chain.gz')

def export_peaks_to_bed(txt_file, out_bed, name):
    print(f"Exporting {txt_file} to {out_bed}...")
    
    count = 0
    with gzip.open(txt_file, 'rt') as f_in, open(out_bed, 'w') as f_out:
        for line in f_in:
            parts = line.strip().split()
            if len(parts) >= 1:
                loc = parts[0]
                if ':' not in loc or '-' not in loc:
                    continue
                chrom, coords = loc.split(':')
                start_str, end_str = coords.split('-')
                
                if not chrom.startswith('chr'):
                    chrom = 'chr' + chrom
                
                try:
                    start_mm8 = int(start_str)
                    end_mm8 = int(end_str)
                    
                    start_res = converter[chrom][start_mm8]
                    end_res = converter[chrom][end_mm8]
                    
                    if start_res and end_res:
                        new_chrom = start_res[0][0]
                        new_start = start_res[0][1]
                        new_end = end_res[0][1]
                        
                        if new_chrom == end_res[0][0]:
                            # Output in BED format
                            # chr start end name score strand
                            # Score is arbitrary (e.g. 1000)
                            f_out.write(f"{new_chrom}\t{new_start}\t{new_end}\t{name}_{count}\t1000\t+\n")
                            count += 1
                except ValueError:
                    continue

    print(f"  -> Exported {count} peaks to {out_bed}")

if __name__ == '__main__':
    export_peaks_to_bed('../data/raw/GSM288347_ES_Sox2.txt.gz', '../data/processed/SOX2_peaks_mm9.bed', 'SOX2')
    export_peaks_to_bed('../data/raw/GSM288345_ES_Nanog.txt.gz', '../data/processed/NANOG_peaks_mm9.bed', 'NANOG')
    print("Done!")
