import argparse
import os
import sys
from peaks import parse_peaks
from genome import load_genome
from io_utils import write_tf_sequences as write_sequences


def parse_args():
    """
    Define y analiza los argumentos de línea de comandos.

    Retorna:
        argparse.Namespace: Objeto con los argumentos proporcionados por el usuario.
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
    args = parse_args()

    # Convertir rutas a absolutas
    peak_file = os.path.abspath(args.peaks)
    genome_file = os.path.abspath(args.genome)
    output_dir = os.path.abspath(args.outdir)

    # Validar archivos y directorio de salida
    if not os.path.exists(peak_file):
        print(f"ERROR: Archivo de picos no encontrado: {peak_file}", file=sys.stderr)
        sys.exit(1)
    if not os.path.exists(genome_file):
        print(f"ERROR: Archivo del genoma no encontrado: {genome_file}", file=sys.stderr)
        sys.exit(1)
    os.makedirs(output_dir, exist_ok=True)

    print("\n" + "="*60)
    print(f"{'EXTRACCION DE SITIOS DE UNION':^60}")
    print("="*60)
    print(f"Archivo de picos: {peak_file}")
    print(f"Archivo del genoma: {genome_file}")
    print(f"Directorio de salida: {output_dir}")
    print("="*60)

    # Paso 1: Procesar archivo de picos
    print("\nPaso 1/2: Procesando archivo de picos...")
    peaks_dict = parse_peaks(peak_file)
    if not peaks_dict:
        print("ERROR: No se encontraron picos validos en el archivo", file=sys.stderr)
        sys.exit(1)
    print(f"Encontrados {len(peaks_dict)} factores de transcripcion con picos validos")

    # Paso 2: Cargar genoma y extraer secuencias
    print("\nPaso 2/2: Extrayendo secuencias...")
    genome_seq, genome_length = load_genome(genome_file)
    print(f"Genoma cargado (longitud: {genome_length} pb)")

    write_sequences(peaks_dict, genome_seq, output_dir)

    # Confirmación final
    print("\n" + "="*60)
    print(f"{'PROCESO COMPLETADO CON EXITO':^60}")
    print("="*60)

if __name__ == '__main__':
    main()
