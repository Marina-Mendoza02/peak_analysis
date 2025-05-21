from Bio import SeqIO
import sys

def load_genome(genome_file_path):
    """
    Carga la secuencia del genoma desde un archivo FASTA.
    
    Args:
        genome_file_path (str): Ruta al archivo FASTA del genoma.
        
    Returns:
        SeqRecord: Registro del genoma cargado.
    """
    try:
        # Abrir y parsear el FASTA, se espera un solo registro (genoma)
        genome_record = next(SeqIO.parse(genome_file_path, 'fasta'))

        # Obtener longitud del genoma para verificación
        genome_length = len(genome_record.seq)
    
        # Devolver registro para uso posterior
        return genome_record, genome_length

    except Exception as e:
        print(f"ERROR al cargar genoma: {str(e)}", file=sys.stderr)
        raise # para manejo en main.py
