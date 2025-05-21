import os
import sys
from Bio.SeqRecord import SeqRecord
from io_utils import write_tf_fastas

def extract_sequences(genome_record: SeqRecord, peaks_dict: dict, output_dir: str):
    """
    Extrae secuencias del genoma para cada pico asociado a un factor de transcripción (TF)
    y guarda las secuencias en archivos FASTA separados por TF.

    Parámetros:
        genome_record (SeqRecord): Objeto SeqRecord que contiene la secuencia genómica completa.
        peaks_dict (dict): Diccionario con picos agrupados por factor de transcripción.
                           Formato esperado: {tf_name: [(start, end, peak_id), ...]}
                           Las coordenadas start y end son 1-based.
        output_dir (str): Ruta al directorio donde se guardarán los archivos FASTA generados.

    Proceso:
        - Valida que las coordenadas de cada pico estén dentro de los límites de la secuencia genómica.
        - Omite picos con coordenadas inválidas y emite advertencias.
        - Para cada TF, escribe un archivo FASTA con las secuencias correspondientes a sus picos.
        - Muestra un resumen del número de secuencias guardadas por TF y totales.
    """
    genome_seq = genome_record.seq
    genome_length = len(genome_seq)

    print(f"\nProcesando {len(peaks_dict)} factores de transcripción...")

    summary = {}

    for tf_name, peaks in peaks_dict.items():
        sequences = []
        for start, end, peak_id in peaks: 
            # Convertir coordenadas (1-based → 0-based para Python)
            start_idx = start - 1
            end_idx = end   # indice final exclusivo

            # Validar coordenadas
            if start_idx < 0 or end_idx > genome_length:
                print(f"Advertencia: Pico {peak_id} ({start}-{end}) fuera de los límites. Se omite.", file=sys.stderr)
                continue
            
            sequence = genome_seq[start_idx:end_idx]
            sequences.append((peak_id, start, end, str(sequence)))

        if sequences:
            try:
                write_tf_fastas(tf_name, sequences, output_dir)
                summary[tf_name] = len(sequences)
                print(f"Guardadas {len(sequences)} secuencias para {tf_name} → {tf_name}.fa")
            except IOError as e:
                print(f"ERROR: No se pudo escribir en {tf_name}.fa: {str(e)}", file=sys.stderr)

    total_tfs = len(summary)
    total_seqs = sum(summary.values())

    print(f"\nExtracción completada: {total_seqs} secuencias de {total_tfs} factores de transcripción guardadas.")

    # para debugging
    example_tfs = list(summary.items())[:10]
    print("Ejemplo de conteo por TF:")
    for tf, count in example_tfs:
        print(f"  {tf}: {count} secuencias")

    if total_tfs > 10:
        print(f"  ... y {total_tfs - 10} TFs adicionales.")
