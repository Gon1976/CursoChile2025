from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import pandas as pd

# Archivos de entrada
blast_file = "CURSO/resultados_blast.txt"
genome_file = "CURSO/TcDm25_genome.fasta"

# Leer genoma en diccionario
genome = SeqIO.to_dict(SeqIO.parse(genome_file, "fasta"))

# Leer resultados de BLAST (formato tabular)
cols = ["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
        "qstart", "qend", "sstart", "send", "evalue", "bitscore"]
blast_df = pd.read_csv(blast_file, sep="\t", names=cols)

# Lista para secuencias de ADN y proteínas
dna_records = []
protein_records = []

for i, row in blast_df.iterrows():
    chrom = row["sseqid"]
    start = int(row["sstart"])
    end = int(row["send"])
    qname = row["qseqid"]

    # Asegurar coordenadas crecientes
    seq_start = min(start, end)
    seq_end = max(start, end)
    strand = "+" if start < end else "-"

    # Extraer secuencia
    seq = genome[chrom].seq[seq_start - 1 : seq_end]  # índices desde 0
    if strand == "-":
        seq = seq.reverse_complement()

    # Generar ID
    seq_id = f"{chrom}_{seq_start}-{seq_end}_{qname}"

    # Crear record ADN
    dna_rec = SeqRecord(seq, id=seq_id, description="")
    dna_records.append(dna_rec)

    # Traducir a proteína (primera fase de lectura, sin buscar ORF completa)
    protein_seq = seq.translate(to_stop=True)
    prot_rec = SeqRecord(protein_seq, id=seq_id, description="")
    protein_records.append(prot_rec)

# Guardar multifastas
SeqIO.write(dna_records, "CURSO/hits_dna.fasta", "fasta")
SeqIO.write(protein_records, "CURSO/hits_proteins.fasta", "fasta")

print("✅ Multifastas generados:")
print(" - CURSO/hits_dna.fasta")
print(" - CURSO/hits_proteins.fasta")
