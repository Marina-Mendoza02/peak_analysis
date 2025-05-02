#!/usr/bin/env python3
# Script para extraer secuencias de union a factores de transcripcion

import argparse
import os        # Para manejo de rutas y sistema de archivos
import sys       # Para salida de errores y codigos de salida
from collections import defaultdict
from Bio import SeqIO   # Para leer archivos FASTA

def parse_peaks(peak_file_path):
    """
    Analiza el archivo de picos y devuelve un diccionario con el formato:
    {tf_name: [(start, end, peak_id), ...]}
    """ 
    # Diccionario para agrupar picos por factor de transcripcion
    peaks_by_tf = defaultdict(list)

    with open(peak_file_path) as f:
        # Leer el encabezado del archivo TSV
        header = f.readline().strip().split('\t')
        
        # Columnas obligatorias que debe tener el archivo
        required_columns = {'TF_name', 'Peak_start', 'Peak_end', 'Dataset_Ids', 'Peak_number'}
        
        # Verificar que todas las columnas requeridas estén presentes
        if not required_columns.issubset(set(header)):
            missing = required_columns - set(header)
            print(f"ERROR: Faltan columnas requeridas: {missing}", file=sys.stderr)
            sys.exit(1)  # Terminar el programa si faltan columnas

        # Procesar cada línea del archivo (comenzando desde la línea 2)
        for line_num, line in enumerate(f, 2):
            fields = line.strip().split('\t')
            
            # Verificar que el número de campos coincida con la cabecera
            if len(fields) != len(header):
                print(f"Advertencia: Linea {line_num} tiene {len(fields)} campos (se esperaban {len(header)}). Se omite.", file=sys.stderr)
                continue
            
            try:
                # Crear diccionario con los datos del pico
                peak_data = dict(zip(header, fields))
                
                # Extraer información clave
                tf_name = peak_data['TF_name']  # Nombre del factor de transcripción
                start = int(float(peak_data['Peak_start']))  # Inicio del pico
                end = int(float(peak_data['Peak_end']))      # Fin del pico
                
                # Validar que las coordenadas sean correctas
                if start > end:
                    print(f"Advertencia: Linea {line_num} tiene coordenadas invalidas ({start} > {end}). Se omite.", file=sys.stderr)
                    continue
                
                # Crear ID unico para el pico
                peak_id = f"{peak_data['Dataset_Ids']}_{peak_data['Peak_number']}"
                
                # Almacenar en el diccionario
                peaks_by_tf[tf_name].append((start, end, peak_id))
                
            except (ValueError, KeyError) as e:
                print(f"Advertencia: Error procesando linea {line_num}: {str(e)}. Se omite.", file=sys.stderr)
                continue

    return peaks_by_tf

def extract_sequences(genome_file_path, peaks_dict, output_dir):
    """
    Extrae secuencias del genoma y guarda archivos FASTA por cada TF
    """
    try:
        # Cargar el genoma (se asume un solo cromosoma/contig)
        genome_record = next(SeqIO.parse(genome_file_path, 'fasta'))
        genome_seq = genome_record.seq  # Secuencia completa del genoma
        genome_length = len(genome_seq)  # Longitud del genoma
        
        print(f"\nGenoma cargado (longitud: {genome_length} pb)")
        print(f"Procesando {len(peaks_dict)} factores de transcripcion...")

        # Procesar cada factor de transcripcion
        for tf_name, peaks in peaks_dict.items():
            # Crear nombre de archivo seguro (reemplazando caracteres especiales)
            safe_name = "".join(c for c in tf_name if c.isalnum() or c in ('_', '-')).rstrip()
            output_path = os.path.join(output_dir, f"{safe_name}.fa")  # Ruta completa del archivo de salida
            sequences_extracted = 0  # Contador de secuencias
            
            try:
                # Abrir archivo FASTA de salida
                with open(output_path, 'w') as out_file:
                    # Procesar cada pico para este TF
                    for start, end, peak_id in peaks:
                        # Convertir coordenadas (1-based → 0-based para Python)
                        start_idx = start - 1
                        end_idx = end  # indice final exclusivo
                        
                        # Validar coordenadas
                        if start_idx < 0 or end_idx > genome_length:
                            print(f"Advertencia: Pico {peak_id} con coordenadas ({start}-{end}) fuera de los limites. Se omite.", file=sys.stderr)
                            continue
                            
                        # Extraer la secuencia
                        sequence = genome_seq[start_idx:end_idx]
                        
                        # Escribir en formato FASTA
                        out_file.write(f">{peak_id}|{start}-{end}\n{sequence}\n")
                        sequences_extracted += 1
                        
                # Reportar éxito para este TF
                print(f"Guardadas {sequences_extracted} secuencias para {tf_name} en {os.path.basename(output_path)}")
                
            except IOError as e:
                print(f"ERROR: No se pudo escribir en {output_path}: {str(e)}", file=sys.stderr)
                continue
                
    except Exception as e:
        print(f"ERROR CRITICO: {str(e)}", file=sys.stderr)
        sys.exit(1)

def main():
    # Configurar el parser de argumentos de linea de comandos
    parser = argparse.ArgumentParser(
        description='Extraer secuencias de union a TFs usando datos de ChIP-Seq',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    
    # Definir argumentos requeridos
    parser.add_argument('peak_file', help='Archivo TSV con informacion de picos ChIP-Seq')
    parser.add_argument('genome_file', help='Archivo FASTA con la secuencia del genoma')
    parser.add_argument('output_dir', help='Directorio para los archivos FASTA de salida')

    # Procesar argumentos
    args = parser.parse_args()

    try:
        # Convertir rutas a absolutas
        args.peak_file = os.path.abspath(args.peak_file)
        args.genome_file = os.path.abspath(args.genome_file)
        args.output_dir = os.path.abspath(args.output_dir)

        # Validar que los archivos de entrada existan
        if not os.path.exists(args.peak_file):
            raise FileNotFoundError(f"Archivo de picos no encontrado: {args.peak_file}")
        if not os.path.exists(args.genome_file):
            raise FileNotFoundError(f"Archivo del genoma no encontrado: {args.genome_file}")
        
        # Crear directorio de salida si no existe
        os.makedirs(args.output_dir, exist_ok=True)

        # Mostrar información de configuración
        print("\n" + "="*60)
        print(f"{'EXTRACCION DE SITIOS DE UNION':^60}")
        print("="*60)
        print(f"Archivo de picos: {args.peak_file}")
        print(f"Archivo del genoma: {args.genome_file}")
        print(f"Directorio de salida: {args.output_dir}")
        print("="*60)

        # Paso 1: Procesar archivo de picos
        print("\nPaso 1/2: Procesando archivo de picos...")
        peaks_dict = parse_peaks(args.peak_file)
        
        if not peaks_dict:
            raise ValueError("No se encontraron picos validos en el archivo")
        
        print(f"Encontrados {len(peaks_dict)} factores de transcripcion con picos validos")
        
        # Paso 2: Extraer secuencias
        print("\nPaso 2/2: Extrayendo secuencias...")
        extract_sequences(args.genome_file, peaks_dict, args.output_dir)

        # Mensaje de finalizacion exitosa
        print("\n" + "="*60)
        print(f"{'PROCESO COMPLETADO CON EXITO':^60}")
        print("="*60)

    except Exception as e:
        print(f"\nERROR: {str(e)}", file=sys.stderr)
        sys.exit(1)

if __name__ == '__main__':
    main()