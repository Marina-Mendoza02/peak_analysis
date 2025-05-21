from collections import defaultdict
import sys       # Para salida de errores y codigos de salida

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