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
        # Leer encabezado y sus columnas
        header = f.readline().strip().split('\t')

        # Verificar que existan las columnas necesarias
        required_columns = {'TF_name', 'Peak_start', 'Peak_end', 'Dataset_Ids', 'Peak_number'}
        if not required_columns.issubset(set(header)):
            missing = required_columns - set(header)
            print(f"ERROR: Faltan columnas requeridas en el archivo de picos: {missing}", file=sys.stderr)
            sys.exit(1) #interrumpe el proceso si faltan columnas
        

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

def extract_sequences(genome_file_path, peaks_dict, output_dir):
    """
    Extrae secuencias del genoma y guarda archivos FASTA por cada TF
    """
    # Cargar el genoma (asumiendo un solo contig/cromosoma)
    try:
        genome_record = next(SeqIO.parse(genome_file_path, 'fasta'))
    except Exception as e:
        print(f"ERROR: No se pudo leer el archivo FASTA del genoma: {str(e)}", file=sys.stderr)
        sys.exit(1)

    genome_seq = genome_record.seq
    genome_length = len(genome_seq)
    
    print(f"Genoma cargado con longitud de {genome_length} pb")
    print(f"Procesando {len(peaks_dict)} factores de transcripcion...\n")
    
    # Procesar cada factor de transcripción
    for tf_name, peaks in peaks_dict.items():
        safe_tf_name = "".join(c for c in tf_name if c.isalnum() or c in ('_','-')).rstrip()
        output_path = os.path.join(output_dir, f"{tf_name.replace('/', '_')}.fa")

        sequences_extracted = 0
        
        try:
            with open(output_path, 'w') as out_file:
                for start, end, peak_id in peaks:
                    # Ajustar coordenadas (1-based a 0-based en Python)
                    start_idx = start - 1
                    end_idx = end  # slicing en Python excluye el índice final
                
                    # Validar coordenadas
                    if start_idx < 0 or end_idx > genome_length:
                        print(f"Advertencia: Pico {peak_id} con coordenadas ({start}-{end}) "
                          f"fuera de los límites del genoma (1-{genome_length}). Se omite.")
                        continue
                
                    # Extraer secuencia
                    sequence = genome_seq[start_idx:end_idx]
                
                    # Escribir en formato FASTA
                    out_file.write(f">{peak_id}|{start}-{end}\n{sequence}\n")
                    sequences_extracted += 1
        
            print(f"{sequences_extracted} secuencias guardadas para {tf_name} en {output_path}")

            except IOError as e:
                print(f"ERROR: No se pudo escribir {output_path}: {str(e)}", file=sys.stderr)
                continue

def main():
    # Analizar de argumentos de linea de comandos
    parser = argparse.ArgumentParser(
        description = 'Extraer secuencias de union a TFs desde el genoma usando datos de picos',
        formatter_class = argparse.ArgumentDefaultsHelpFormatter)
    
    parser.add_argument('peak_file', help='Archivo TSV con informacion de picos')
    parser.add_argument('genome_file', help='Archivo FASTA con la secuencia del genoma')
    parser.add_argument('output_dir', help='Directorio para guardar los archivos FASTA de cada TF')


    args = parser.parse_args()

    # Convertir a rutas absolutas
    try:
        args.peak_file = os.path.abspath(args.peak_file)
        args.genome_file = os.path.abspath(args.genome_file)
        args.output_dir = os.path.abspath(args.output_dir)


        # Validar entradas
        if not os.path.exists(args.peak_file):
            raise FileNotFoundError(f"Peak file not found: {args.peak_file}")
        if not os.path.exists(args.genome_file):
            raise FileNotFoundError(f"Genome file not found: {args.genome_file}")
        
        os.makedirs(args.output_dir, exist_ok=True)

    except Exception as e:
        print(f"Error de inicializacion: {str(e)}", file=sys.stderr)
        sys.exit(1)

    print("\n" + "="*60)
    print(f"{'TF Binding Site Extraction':^60}")
    print("="*60)
    print(f"Input peak file: {args.peak_file}")
    print(f"Input genome file: {args.genome_file}")
    print(f"Output directory: {args.output_dir}")
    print("="*60 + "\n")

    # Procesar datos
    try:
        # ========== MAIN FUNCTION CALLS ==========
        print("Paso 1/2: Analizando archivo de picos...")
        peaks_dict = parse_peaks(args.peak_file)
        
        if not peaks_dict:
            raise ValueError("No valid peaks found in input file")
        
        print(f"Se encontraron {len(peaks_dict)} TFs con picos validos")
        
        print("\nPaso 2/2: Extrayendo secuencias...")
        extract_sequences(args.genome_file, peaks_dict, args.output_dir)
        # ========================================

        print("\n" + "="*60)
        print(f"{'Procesamiento completado exitosamente.':^60}")
        print("="*60)

    except Exception as e:
        print(f"\nError: {str(e)}", file=sys.stderr)
        sys.exit(1)

if __name__ == '__main__':
    main()