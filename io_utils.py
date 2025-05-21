import os

def write_tf_fastas(tf_name, sequences, output_dir):
    """
    Guarda las secuencias extraídas en un archivo FASTA, agrupadas por factor de transcripción (TF).

    Parámetros:
        tf_name (str): Nombre del factor de transcripción para usar como nombre de archivo.
        sequences (lista de tuplas): Lista de secuencias a escribir, cada una representada como una tupla:
                                   (peak_id, start, end, seq), donde:
                                   - peak_id (str): Identificador único del pico.
                                   - start (int): Posición inicial (1-based) del pico en el genoma.
                                   - end (int): Posición final (1-based, inclusiva) del pico.
                                   - seq (str): Secuencia de nucleótidos correspondiente al pico.
        output_dir (str): Ruta al directorio donde se guardará el archivo FASTA.

    Comportamiento:
        - Crea un nombre de archivo seguro a partir de tf_name, removiendo caracteres no permitidos.
        - Escribe un archivo FASTA con encabezados formateados como '>peak_id|start-end'.
        - Captura y reporta errores de escritura en la salida de error estándar.

    Excepciones:
        - No lanza excepciones, pero imprime errores de IO si no puede escribir el archivo.
    """
    # Crear nombre seguro para el archivo
    safe_name = "".join(c for c in tf_name if c.isalnum() or c in ('_', '-')).rstrip()
    output_path = os.path.join(output_dir, f"{safe_name}.fa")
    
    try:
        with open(output_path, 'w') as out_file:
            for peak_id, start, end, seq in sequences:
                out_file.write(f">{peak_id}|{start}-{end}\n{seq}\n")
    except IOError as e:
        import sys
        print(f"ERROR: No se pudo escribir en {output_path}: {str(e)}", file=sys.stderr)
