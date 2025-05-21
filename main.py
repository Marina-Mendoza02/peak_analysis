import argparse
import os
import sys
from peaks import parse_peaks
from genome import load_genome
from extractor import extract_sequences as write_sequences

# Ejemplo de uso
# Con archivo de picos más corto: python main.py -p data/union_peaks_file_short.tsv -g data/E_coli_K12_MG1655_U00096.3.fasta -o results/
# Archivo de picos completo: python main.py -p data/union_peaks_file.tsv -g data/E_coli_K12_MG1655_U00096.3.fasta -o results/


def parse_args():
    """
    Configura y analiza los argumentos de línea de comandos para el script.
    
    Define tres argumentos obligatorios:
    - --peaks (-p): Archivo TSV con datos de picos ChIP-Seq
    - --genome (-g): Archivo FASTA con el genoma de referencia
    - --outdir (-o): Directorio de salida para los resultados
    
    Returns:
        argparse.Namespace: Objeto con los argumentos parseados
    
    Nota:
        Los argumentos originalmente eran posicionales, pero se cambiaron a argumentos
        opcionales con flags requeridas para mejorar la claridad y flexibilidad del script.
    """
    parser = argparse.ArgumentParser(
        description='Extrae secuencias de union a TFs usando datos de ChIP-Seq',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser = argparse.ArgumentParser(
        description='Extrae secuencias de unión a factores de transcripción (TFs) usando datos de ChIP-Seq',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        '-p', '--peaks', required=True,
        help='Archivo TSV con información de picos de ChIP-Seq (formato: TF, start, end, ID)'
    )
    parser.add_argument(
        '-g', '--genome', required=True,
        help='Archivo FASTA con la secuencia del genoma de referencia'
    )
    parser.add_argument(
        '-o', '--outdir', required=True,
        help='Directorio donde se guardarán los archivos FASTA por TF'
    )
    return parser.parse_args()

def main():
    """
    Función principal que orquesta el proceso completo:
    1. Parsea los argumentos de línea de comandos
    2. Valida los archivos de entrada
    3. Procesa los picos de ChIP-Seq
    4. Carga el genoma de referencia
    5. Extrae y guarda las secuencias de unión
    """
    # Paso 1: Configuración inicial
    args = parse_args()

    # Convertir rutas a absolutas
    peak_file = os.path.abspath(args.peaks)
    genome_file = os.path.abspath(args.genome)
    output_dir = os.path.abspath(args.outdir)

    # Paso 2: Validación de entradas
    if not os.path.exists(peak_file):
        print(f"ERROR: Archivo de picos no encontrado: {peak_file}", file=sys.stderr)
        sys.exit(1)
    if not os.path.exists(genome_file):
        print(f"ERROR: Archivo del genoma no encontrado: {genome_file}", file=sys.stderr)
        sys.exit(1)

    # Crear directorio de salida si no existe
    os.makedirs(output_dir, exist_ok=True)

    print("\n" + "="*60)
    print(f"{'EXTRACCION DE SITIOS DE UNION':^60}")
    print("="*60)
    print(f"Archivo de picos: {peak_file}")
    print(f"Archivo del genoma: {genome_file}")
    print(f"Directorio de salida: {output_dir}")
    print("="*60)

    # Paso 3: Procesar archivo de picos
    print("\nPaso 1/2: Procesando archivo de picos...")
    peaks_dict = parse_peaks(peak_file)
    if not peaks_dict:
        print("ERROR: No se encontraron picos validos en el archivo", file=sys.stderr)
        sys.exit(1)
    print(f"Encontrados {len(peaks_dict)} factores de transcripcion con picos validos")

    # Paso 4: Cargar genoma
    print("\nPaso 2/2: Extrayendo secuencias...")
    genome_record = load_genome(genome_file)
    genome_seq = genome_record.seq
    genome_length = len(genome_seq)

    print(f"Genoma cargado (longitud: {genome_length} pb)")

    # Paso 5. Extraer secuencias
    write_sequences(genome_seq, peaks_dict, output_dir)

    # Mensaje final
    print("\n" + "="*60)
    print(f"{'PROCESO COMPLETADO CON EXITO':^60}")
    print("="*60)

if __name__ == '__main__':
    main()
