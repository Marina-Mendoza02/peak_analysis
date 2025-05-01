#!/usr/bin/env python3

import argparse
import os
from collections import defaultdict
from Bio import SeqIO # Para el manejo de secuencias

def main():
    # Analizar de argumentos de linea de comandos
    parser = argparse.ArgumentParser(
        description = 'Extraer secuencias de union a TFs desde el genoma usando datos de picos'
        formatter_class = argparse.ArgumentDefaultsHelpFormatter)
    
    parser.add_argument('peak_file', help='Archivo TSV con informacion de picos')
    parser.add_argument('genome_file', help='Archivo FASTA con la secuencia del genoma')
    parser.add_argument('output_dir', help='Directorio para guardar los archivos FASTA de cada TF')

    args = parser.parse_args()

    # Validar entradas
    if not os.path.exists(args.peak_file):
        raise FileNotFoundError(f"No se encontro el archivo de picos: {args.peak_file}")
    if not os.path.exists(args.genome_file):
        raise FileNotFoundError(f"No se encontro el archivo del genoma: {args.genome_file}")
    
    # Crear directorio de salida si no existe
    os.makedirs(args.output_dir, exist_ok=True)

    print(f'Procesando picos desde {args.peak_file}')
    print(f'Usando genoma de {args.genome_file}')
    print(f'La salida se guardara en {args.output_dir}')

if __name__=='__main__':
    main()