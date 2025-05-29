# README - Extracción de Secuencias de Unión a Factores de Transcripción

## Descripción

Este proyecto permite extraer secuencias genómicas asociadas a sitios de unión de factores de transcripción (TFs) en *Escherichia coli* a partir de datos de experimentos ChIP-Seq. Genera archivos FASTA individuales para cada factor de transcripción, listos para su uso en análisis posteriores.

## Requisitos

- Python 3.6 o superior
- Biopython
- Sistema operativo Unix/Linux

## Instalación

1. Clonar el repositorio:
```bash
git clone https://github.com/Marina-Mendoza02/peak_analysis.git
cd peak_analysis
```

2. Instalar dependencias:
```bash
pip install biopython
```

## Uso

### Ejecución básica

```bash
python main.py -p archivo_picos.tsv -g genoma.fasta -o directorio_salida/
```

### Parámetros

| Argumento | Descripción | Requerido |
|-----------|-------------|-----------|
| `-p`, `--peaks` | Ruta al archivo TSV con datos de picos ChIP-Seq | Sí |
| `-g`, `--genome` | Ruta al archivo FASTA del genoma de referencia | Sí |
| `-o`, `--outdir` | Directorio donde se guardarán los archivos FASTA | Sí |

### Formato del archivo de picos

El archivo TSV debe contener al menos estas columnas:
- `TF_name`: Nombre del factor de transcripción
- `Peak_start`: Posición inicial del pico (1-based)
- `Peak_end`: Posición final del pico (1-based)
- `Dataset_Ids`: Identificador del dataset
- `Peak_number`: Número identificador del pico

## Estructura del proyecto

```
peak_analysis/
├── data/                     # Datos de entrada
│   ├── E_coli_genome.fasta   # Genoma de referencia
│   └── peaks_data.tsv        # Datos de picos ChIP-Seq
│
├── src/                      # Código fuente
│   ├── main.py               # Punto de entrada principal
│   ├── genome.py             # Manejo del genoma
│   ├── peaks.py              # Procesamiento de picos
│   ├── extractor.py          # Extracción de secuencias
│   └── io_utils.py           # Utilidades de entrada/salida
│
├── results/                  # Resultados generados por el programa
│
├── docs/                     # Documentación del proyecto
│   ├── test_cases.md         # Casos de prueba manuales
│   └── project_overview.md   # Descripción general del proyecto
│
├── README.md                 # Este archivo
├── LICENSE                   # Licencia del proyecto
└── .gitignore                # Exclusiones de seguimiento en Git

```

## Ejemplo de salida

Para cada factor de transcripción se genera un archivo `.fa` con el formato:

```
>ID_pico|inicio-fin
SECUENCIA_EXTRACTA
```

Ejemplo (`AraC.fa`):
```
>DS001_1|12345-12400
ATGACATCATGACATCATGACATCATGACATCATGACATCATGACATCATGA
>DS001_2|45600-45650
CGTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA
```

## Manejo de errores

El programa detecta y reporta:
- Archivos de entrada faltantes o corruptos
- Coordenadas inválidas (fuera del rango del genoma)
- Problemas de permisos en el directorio de salida
- Errores de formato en los archivos de entrada

Los mensajes de error se muestran en la salida estándar de error (stderr).

## Licencia

Este proyecto está bajo la licencia [MIT](LICENSE).

## Contacto

Para preguntas o comentarios, contactar a:
Marina Mendoza Suárez <marinams@lcg.unam.mx>