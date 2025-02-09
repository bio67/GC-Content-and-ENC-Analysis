from Bio import SeqIO
import csv

# Define synonymous codons dictionary
SYNONYMOUS_CODONS = {
    'Phe': ['TTT', 'TTC'], 'Leu': ['TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG'],
    'Ile': ['ATT', 'ATC', 'ATA'], 'Val': ['GTT', 'GTC', 'GTA', 'GTG'],
    'Ser': ['TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC'], 'Pro': ['CCT', 'CCC', 'CCA', 'CCG'],
    'Thr': ['ACT', 'ACC', 'ACA', 'ACG'], 'Ala': ['GCT', 'GCC', 'GCA', 'GCG'],
    'Tyr': ['TAT', 'TAC'], 'His': ['CAT', 'CAC'], 'Gln': ['CAA', 'CAG'],
    'Asn': ['AAT', 'AAC'], 'Lys': ['AAA', 'AAG'], 'Asp': ['GAT', 'GAC'],
    'Glu': ['GAA', 'GAG'], 'Cys': ['TGT', 'TGC'], 'Trp': ['TGG'],
    'Arg': ['CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'], 'Gly': ['GGT', 'GGC', 'GGA', 'GGG']
}

def calculate_gc_content(sequence):
    """Calculate GC1, GC2, GC3, GC3s, and GC12 values"""
    codons = [sequence[i:i+3] for i in range(0, len(sequence) - len(sequence) % 3, 3)]
    gc1, gc2, gc3 = 0, 0, 0
    gc3s_count, synonymous_count = 0, 0

    for codon in codons:
        if len(codon) < 3:  # Skip incomplete codons
            continue
        # GC content at the 1st and 2nd positions
        if codon[0] in "GC":
            gc1 += 1
        if codon[1] in "GC":
            gc2 += 1
        # GC content at the 3rd position
        if codon[2] in "GC":
            gc3 += 1
        # Check if the codon is synonymous
        for synonymous_codons in SYNONYMOUS_CODONS.values():
            if codon.upper() in synonymous_codons:
                synonymous_count += 1
                if codon[2] in "GC":
                    gc3s_count += 1

    total_codons = len(codons)
    gc1_percentage = (gc1 / total_codons) * 100 if total_codons > 0 else 0
    gc2_percentage = (gc2 / total_codons) * 100 if total_codons > 0 else 0
    gc3_percentage = (gc3 / total_codons) * 100 if total_codons > 0 else 0
    gc3s_percentage = (gc3s_count / synonymous_count) * 100 if synonymous_count > 0 else 0
    gc12_percentage = (gc1_percentage + gc2_percentage) / 2  # Calculate GC12

    return gc1_percentage, gc2_percentage, gc3_percentage, gc3s_percentage, gc12_percentage


def process_fasta_and_calculate_gc(input_fasta, output_csv):
    """Read FASTA file, calculate GC1, GC2, GC3, GC3s, GC12, and save to CSV file"""
    with open(output_csv, 'w', newline='') as csvfile:
        fieldnames = ['Gene', 'GC1 (%)', 'GC2 (%)', 'GC3 (%)', 'GC3s (%)', 'GC12 (%)']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()

        for record in SeqIO.parse(input_fasta, "fasta"):
            sequence = str(record.seq).upper()
            gc1, gc2, gc3, gc3s, gc12 = calculate_gc_content(sequence)
            writer.writerow({
                'Gene': record.id,
                'GC1 (%)': round(gc1, 2),
                'GC2 (%)': round(gc2, 2),
                'GC3 (%)': round(gc3, 2),
                'GC3s (%)': round(gc3s, 2),
                'GC12 (%)': round(gc12, 2)
            })


# Input and output files
input_fasta = "name.fasta"  # Replace with your FASTA file path
output_csv = "gc_results.csv"  # Output results file path

# Run the script
process_fasta_and_calculate_gc(input_fasta, output_csv)
print(f"GC1, GC2, GC3, GC3s, and GC12 results saved to {output_csv}")