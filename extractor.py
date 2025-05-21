import os
import sys
from io_utils import write_tf_fastas

def extract_sequences(genome_seq, peaks_dict, output_dir):
    """
    Extrae secuencias del genoma y guarda archivos FASTA por cada TF usando io_utils.
    
    Args:
        genome_seq (Seq): Secuencia completa del genoma (Bio.Seq.Seq).
        peaks_dict (dict): {tf_name: [(start, end, peak_id), ...], ...}
        output_dir (str): Directorio donde se guardan los FASTA.
    
    Proceso:
        1. Calcula la longitud total del genoma para validación de coordenadas
        2. Itera sobre cada TF y sus picos asociados
        3. Para cada pico:
           a. Convierte coordenadas 1-based a 0-based
           b. Valida que las coordenadas estén dentro de los límites del genoma
           c. Extrae la secuencia correspondiente
        4. Escribe un archivo FASTA por TF usando la función write_tf_fastas
    
    Notas:
        - Las coordenadas en el archivo de entrada son 1-based (formato estándar en genómica)
        - Las coordenadas para extracción son convertidas a 0-based (requerido por Python)
        - Los archivos de salida siguen el formato: {TF_name}.fa
    """
    # Longitud del genoma para validación
    genome_length = len(genome_seq)
    
    print(f"\nGenoma cargado (longitud: {genome_length} pb)")
    print(f"Procesando {len(peaks_dict)} factores de transcripcion...")
    
    # Procesamiento por TF
    for tf_name, peaks in peaks_dict.items():
        sequences_to_write = [] # Almacena tuplas (peak_id, start, end, sequence)

        peaks_processed = 0    # Contador de picos

        for start, end, peak_id in peaks:
            # Conversión a coordenadas 0-based para Python
            start_idx = start - 1
            end_idx = end
            
            # Validación de coordenadas
            if start_idx < 0 or end_idx > genome_length:
                import sys
                print(f"Advertencia: Pico {peak_id} con coordenadas ({start}-{end}) fuera de los limites. Se omite.", file=sys.stderr)
                continue
            
            # Extracción de secuencia
            seq = str(genome_seq[start_idx:end_idx]) 
            sequences_to_write.append((peak_id, start, end, seq))
            peaks_processed += 1

        # Escritura de archivo FASTA para el TF actual
        if sequences_to_write:
            write_tf_fastas(tf_name, sequences_to_write, output_dir)
            print(f"{tf_name}: {peaks_processed} secuencias extraidas")
        else:
            print(f"{tf_name}: Sin picos válidos - archivo no generado", file=sys.stderr)
        
    print("\n[COMPLETADO] Extraccion finalizada")