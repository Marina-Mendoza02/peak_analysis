#!/usr/bin/env python3

import argparse
import os
from collections import defaultdict
from Bio import SeqIO # Para el manejo de secuencias

def parse_peaks(peak_file_path):
    """
    Analiza el archivo de picos y devuelve un diccionario con el formato:
    {tf_name: [(start, end, peak_id), ...]}
    """ 

    peaks_by_tf = defaultdict(list)

    with open(peak_file_path) as f:
        # Leer encabezado
        header = f.readline().strip().split('\t')

        # Verificar que existan las columnas necesarias
        required_columns = {'TF_name', 'Peak_start', 'Peak_end', 'Dataset_Ids', 'Peak_number'}
        if not required_columns.issubset(set(header)):
            missing = required_columns - set(header)
            raise ValueError(f"Faltan columnas requeridas en el archivo de picos: {missing}")
        
        # Procesar cada linea del archivo de picos
        for line_num, line in enumerate(f,2): # comenzar desde la linea 2
            fields = line.strip().split('\t')
            if len(fields) != len(header):
                print(f"Advertencia: La linea {line_num} tiene {len(fields)} campos, se esperaban {len(header)}. Se omite.")
                continue
            
            peak_data = dict(zip(header, fields))

            try:
                tf_name = peak_data['TF_name']
                start = int(peak_data['Peak_start'])
                end = int(peak_data['Peak_end'])
                peak_id = f"{peak_data['Dataset_Ids']}_{peak_data['Peak_number']}"

                # Validar coordenadas
                if start > end:
                    print(f"Advertencia: La linea {line_num} tiene start > end ({start} > {end}). Se omite.")
                    continue
                
                peaks_by_tf[tf_name].append((start, end, peak_id))
            except (ValueError, KeyError) as e:
                print(f"Advertencia: Error al procesar la linea {line_num}: {str(e)}. Se omite.")
                continue

    return peaks_by_tf

    

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

