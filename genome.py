from Bio import SeqIO
import sys
import os

def load_genome(genome_file_path):
    """
    Carga la secuencia del genoma desde un archivo FASTA.
    
    Args:
        genome_file_path (str): Ruta al archivo FASTA del genoma.
        
    Returns:
        tuple: (SeqRecord, genome_length) si tiene éxito.
    
    Sale con mensaje de error si:
        - El archivo no existe
        - El archivo está vacío/no es válido
        - Ocurren errores inesperados al analizar
    """
    # Verificar si el archivo existe
    if not os.path.exists(genome_file_path):
        print(f"ERROR: Archivo de genoma no encontrado: {genome_file_path}", file=sys.stderr)
        sys.exit(1)

    try:
        # Intentar analizar el archivo FASTA
        with open(genome_file_path) as f:
            genome_record = next(SeqIO.parse(f, 'fasta'))
            genome_length = len(genome_record.seq)
            return genome_record, genome_length

    except StopIteration:  # Se lanza si el FASTA está vacío
        print(f"ERROR: El archivo de genoma está vacío: {genome_file_path}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:  # Captura otros errores de análisis
        print(f"ERROR: Fallo al analizar el archivo de genoma: {str(e)}", file=sys.stderr)
        sys.exit(1)