from Bio import SeqIO
from Bio.Seq import Seq

# Archivos de entrada y salida
blast_output = "CURSO/resultados_blast.txt"
genome_fasta = "CURSO/TcDm25_genome.fasta"
proteins_fasta = "CURSO/proteinasCurso.fasta"
output_dna = "CURSO/hits_dna.fasta"
output_protein = "CURSO/hits_proteins.fasta"

# Cargar genoma como diccionario
genome_dict = SeqIO.to_dict(SeqIO.parse(genome_fasta, "fasta"))

# Procesar archivo BLAST
dna_records = []
protein_records = []

with open(blast_output) as blast_file:
    for line in blast_file:
        parts = line.strip().split('\t')
        if len(parts) < 12:
            continue  # Saltar líneas incompletas

        qseqid, sseqid, pident, length, mismatch, gapopen, \
        qstart, qend, sstart, send, evalue, bitscore = parts[:12]

        sstart = int(sstart)
        send = int(send)

        # Determinar dirección de la hebra
        strand = '+'
        if sstart > send:
            sstart, send = send, sstart  # Asegurarse de que start < end
            strand = '-'

        # Extraer secuencia
        if sseqid not in genome_dict:
            print(f"WARNING: {sseqid} no encontrado en el genoma.")
            continue

        seq_record = genome_dict[sseqid]
        subseq = seq_record.seq[sstart-1:send]  # BLAST usa 1-based

        # Si está en hebra complementaria: reverse-complement
        if strand == '-':
            subseq = subseq.reverse_complement()

        # Generar header
        header = f"{sseqid}_{sstart}-{send}_{qseqid}_strand{strand}"

        # Guardar secuencia de ADN
        dna_records.append(SeqIO.SeqRecord(subseq, id=header, description=""))

        # Traducir a proteína (sin stop codones incompletos)
        protein_seq = subseq.translate(to_stop=True)
        protein_records.append(SeqIO.SeqRecord(protein_seq, id=header, description=""))

# Escribir multifasta de ADN
with open(output_dna, "w") as out_dna:
    SeqIO.write(dna_records, out_dna, "fasta")

# Escribir multifasta de proteínas
with open(output_protein, "w") as out_prot:
    SeqIO.write(protein_records, out_prot, "fasta")

print("✅ Extracción y traducción completadas.")
