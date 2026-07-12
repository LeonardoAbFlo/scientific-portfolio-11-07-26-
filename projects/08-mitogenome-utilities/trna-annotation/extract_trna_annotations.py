from Bio import SeqIO
import pandas as pd

gb = "barcode81_mtDNA_contig_2.gb"
rec = SeqIO.read(gb, "genbank")

rows = []
with open("tRNAs.fasta","w") as outfa:
    for f in rec.features:
        if f.type == "tRNA":
            seq = f.location.extract(rec.seq)
            start_1based = int(f.location.start) + 1
            end_1based   = int(f.location.end)   # Biopython end is exclusive; this prints inclusive end in 1-based coords
            strand = "+" if f.location.strand == 1 else "-"
            name = (f.qualifiers.get("gene") or f.qualifiers.get("product") or ["tRNA"])[0]
            anticodon = (f.qualifiers.get("anticodon") or [""])[0]  # often like "pos:123..125,aa:Leu"
            ident = f"{name}|{start_1based}-{end_1based}({strand})"
            outfa.write(f">{ident}\n{str(seq)}\n")
            rows.append([name, start_1based, end_1based, strand, anticodon, len(seq)])

pd.DataFrame(rows, columns=["tRNA","start","end","strand","anticodon","length"]).to_csv(
    "tRNAs_coordinates.tsv", sep="\t", index=False
)

# after running this py you acan use this command: arwen -sequence tRNAs.fasta -mt -o barcode83.out
