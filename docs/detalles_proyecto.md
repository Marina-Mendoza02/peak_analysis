# Proyecto de Automatización para la Identificación de Sitios de Unión de Factores de Transcripción en E. coli en experimentos de ChIP-Seq

Fecha: [30/05/2025]

Participantes: 

- [Marina Mendoza Suárez]  <email: marinams@lcg.unam.mx > 

## Descripción del Problema

El proyecto busca automatizar la extracción y el análisis de secuencias genómicas donde los factores de transcripción se unen en _Escherichia coli_. Se cuenta con un archivo que contiene información sobre los picos de unión, y con otro archivo que posee la secuencia completa del genoma. El objetivo es generar archivos FASTA específicos para cada factor de transcripción (TF), agrupando las secuencias de los picos de unión correspondientes. Estas secuencias servirán como entrada para análisis posteriores de motivos de unión utilizando el software 'MEME'.

## Especificación de Requisitos

### Requisitos Funcionales:

#### A. Extracción de Secuencias FASTA:
    
1.  **Entrada de datos:**
    
    -  Implementado en main.py mediante argparse.ArgumentParser
    - Argumentos obligatorios:

        -p/--peaks: Ruta al archivo TSV de picos

        -g/--genome: Ruta al archivo FASTA del genoma

        -o/--outdir: Directorio de salida para FASTA

2. **Validación de entradas**
    - Verificación de existencia de archivos (main.py líneas 52-58)
    - Creación automática del directorio de salida (main.py línea 60)

3. **Procesamiento de picos**
    - Implementado en peaks.py con la función parse_peaks
    - Valida:
        - Presencia de columnas requeridas (TF_name, Peak_start, Peak_end, etc.)
        - Consistencia entre número de campos y cabecera
        - Coordenadas válidas (start ≤ end)

4. **Carga del genoma**
    - Implementado en genome.py con load_genome
    - Verifica:
        - Existencia del archivo
        - Formato FASTA válido
        - Secuencia no vacía

5. **Extracción de secuencias**
    - Implementado en extractor.py con extract_sequences
    - Características:
        - Conversión de coordenadas 1-based a 0-based
        - Validación de límites del genoma
        - Extracción eficiente de subsecuencias

6. **Generación de FASTA**
    - Implementado en io_utils.py con write_tf_fastas
    - Formato de salida:
        - Un archivo por TF (ej: AraC.fa)
        - Encabezados con formato: >peak_id|start-end
        - Nombres de archivo seguros (solo caracteres alfanuméricos y _-)  

### **Requisitos No Funcionales:**

-   **Portabilidad y Usabilidad:**
    
    -   Compatible con sistemas Unix/Linux.
        - Implementado en:
            - Uso de rutas relativas/absolutas con os.path (main.py, genome.py, io_utils.py)
            - Manejo de permisos y creación de directorios con os.makedirs(exist_ok=True) (main.py)
            - Sin dependencias de sistema operativo específico.

    -   El sistema debe ser ejecutable desde la línea de comandos.
        - Implementado en main.py
        - Con validación automática de argumentos:
            - Rutas convertidas a absolutas (os.path.abspath)
            - Verificación de existencia de archivos (os.path.exists)
            - Creación de directorio de salida si no existe (os.makedirs).

    -   Todos los datos de entrada a los programas deben pasarse via argumentos.
        - No hay rutas, nombres de archivo o parámetros fijos en el código.
        - Todo se pasa mediante argparse (-p, -g, -o).

    -   Si se implementa código debe usarse python o scripts shell.
        - Solo requiere Biopython para manejo de FASTA (from Bio import SeqIO).
        - Código autocontenido (no necesita wrappers en Bash).
    
-   **Calidad y Mantenimiento:**
    
    - Estructura de proyecto compatible con Git:
        - Módulos independientes para fácil control de versiones
        - Archivos de requisitos (requirements.txt) para dependencias

    - Documentación clara y comentarios efectivos deben acompañar todo el proyecto.
        - Docstrings detallados en todas las funciones.
        - Comentarios explicativos en partes críticas.

    -   Validación robusta:
        - Validación de entradas
            - Archivo de picos
                - Verificación de columnas obligatorias (TF_name, Peak_start, Peak_end, etc.).
                - Manejo de líneas mal formateadas (peaks.py).
            - Archivo FASTA
                - Detección de archivos vacíos o corruptos (genome.py).
            - Coordenadas
                - Filtrado de picos fuera de los límites del genoma (extractor.py).
        - Manejo de errores
            - Mensajes descriptivos en stderr.


### C. Descripción de Datos de Entrada y Salida 

#### Formato del Archivo de Picos

Este archivo contiene información crucial sobre las regiones de unión de los 144 factores de transcripción (TFs) en _Escherichia coli_. Los datos están organizados en columnas que permiten identificar detalles específicos sobre la unión de los TFs a lo largo del genoma. El formato del archivo y la descripción de cada columna se detallan a continuación:

-   **Dataset_Ids:**
    
    -   _Descripción:_ Identificadores únicos para cada conjunto de datos. Estas IDs indican diferentes experimentos o condiciones bajo las cuales se determinaron los sitios de unión para los TFs.
    -   _Ejemplo:_ "DS001","DS002", etc.
-   **TF_name:**
    
    -   _Descripción:_ El nombre del factor de transcripción que se une al genoma en la región especificada.
    -   _Ejemplo:_ "AraC", "LacI", etc.
-   **Peak_start:**
    
    -   _Descripción:_ La posición inicial en el genoma donde comienza el pico de unión. Se refiere a la ubicación del primer nucleótido del pico.
    -   _Ejemplo:_ 345676, 123456, etc.
-   **Peak_end:**
    
    -   _Descripción:_ La posición final en el genoma donde termina el pico de unión. Se refiere a la ubicación del último nucleótido del pico.
    -   _Ejemplo:_ 345786, 123556, etc.
-   **Peak_center:**
    
    -   _Descripción:_ Posición central del pico de unión, calculada como el promedio o posición entre el `Peak_start` y `Peak_end`.
    -   _Ejemplo:_ 345731, 123501, etc.
-   **Peak_number:**
    
    -   _Descripción:_ Número secuencial utilizado para identificar picos dentro de un conjunto de datos. Esto es útil para referencias internas.
    -   _Ejemplo:_ 1, 2, 3, etc.
-   **Max_Fold_Enrichment:**
    
    -   _Descripción:_ Valor que representa el máximo enriquecimiento observado en el sitio de unión del pico.
    -   _Ejemplo:_ 15.4, 22.3, etc.
-   **Max_Norm_Fold_Enrichment:**
    
    -   _Descripción:_ Valor de máximo enriquecimiento normalizado, ajustado por un factor de control para comparaciones equitativas entre experimentos.
    -   _Ejemplo:_ 12.0, 20.1, etc.
-   **Proximal_genes:**
    
    -   _Descripción:_ Lista de genes cercanos al pico de unión, proporcionando contexto para el análisis funcional.
    -   _Ejemplo:_ "geneA, geneB", "geneX, geneY", etc.
-   **Center_position_type:**
    
    -   _Descripción:_ Denota la ubicación genómica del pico central, como intergénica, intrónica, etc.
    -   _Ejemplo:_ "intergénica", "intrónica", etc.

Debe contener (como mínimo) estas columnas:
- TF_name: Nombre del factor de transcripción

- Peak_start: Inicio del pico (entero, 1-based)

- Peak_end: Fin del pico (entero, 1-based)

- Dataset_Ids: Identificador del dataset

- Peak_number: Número de pico único

#### Salida esperada

1. Archivos FASTA en directorio de salida:

- Formato: {TF_name}.fa
- Ejemplo: AraC.fa, LacI.fa

2. Estructura de archivo FASTA:
```
>DS001_1|12345-12400
AGCTAGCTAGCTAGCTAGCT...
>DS001_2|45600-45680
TTAGCTAGCTAGCTAGCTAG...
```

## Análisis y Diseño

### Arquitectura modular
El sistema se implementó en Python con los siguientes módulos: 
1. main.py:
- Punto de entrada principal
- Manejo de argumentos de línea de comandos
- Coordinación del flujo de trabajo

2. peaks.py:
- Lectura y parseo del archivo TSV de picos
- Validación de datos y estructura
- Agrupación de picos por TF

3. genome.py:
- Carga y validación del archivo FASTA del genoma
- Manejo de errores en formato y contenido

4. extractor.py:
- Extracción de secuencias del genoma
- Validación de coordenadas
- Conversión de formatos de coordenadas

5. io_utils.py:
- Escritura de archivos FASTA
- Formateo de nombres y encabezados
- Manejo seguro de archivos

#### Módulo 1: Extractor y creador de secuencias FASTA

**Objetivo:** Extraer las secuencias genómicas correspondientes a los picos de unión de los factores de transcripción y generar archivos FASTA individuales para cada `TF_name`.

**Flujo de Trabajo:**

1. El usuario proporciona los archivos de entrada

2. El sistema valida y procesa los datos
    -  Módulo 1: Procesamiento de entradas (main.py, peaks.py, genome.py)
        - Valida y parsea los archivos de entrada
        - Agrupa picos por factor de transcripción

3. Se extraen las secuencias genómicas
    - Módulo 2: Extracción de secuencias (extractor.py)
        - Conversión de coordenadas genómicas
        - Validación de rangos
        - Extracción de subsecuencias

4. Se generan los archivos FASTA de salida
    - Módulo 3: Generación de FASTA (io_utils.py)
        - Escritura formateada de archivos
        - Manejo seguro de nombres

5. Detalles Técnicos:

    - Se muestran las conversiones de coordenadas
    - Validaciones críticas (formato, rangos)
    - Operaciones de E/S con manejo de errores


**Pseudocódigo principal**

```
1. Leer argumentos de línea de comandos
2. Validar existencia de archivos de entrada
3. Crear directorio de salida si no existe
4. Parsear archivo de picos:
   a. Leer encabezado y validar columnas
   b. Procesar cada línea:
      i. Extraer TF_name, Peak_start, Peak_end
      ii. Validar coordenadas
      iii. Crear ID único
      iv. Agrupar por TF
5. Cargar genoma:
   a. Verificar formato FASTA
   b. Obtener secuencia completa
6. Para cada TF:
   a. Para cada pico:
      i. Convertir coordenadas a 0-based
      ii. Validar rango
      iii. Extraer secuencia
   b. Escribir archivo FASTA
```

**Pruebas y validación**
Dado que el proyecto no incluye pruebas automatizadas ejecutables (como unit tests con unittest o pytest), se implementó un sistema robusto de validación en tiempo de ejecución que garantiza la calidad del software:

1. Validación de entradas:
- Archivos corruptos o faltantes
- Estructura incorrecta de TSV
- Coordenadas inválidas

2. Pruebas de flujo:
- Ejecución con datos de prueba
- Verificación de archivos de salida
- Comprobación de formatos FASTA

3. Manejo de casos límite:
- Picos en bordes del genoma
- Nombres de TF con caracteres especiales
- Directorios de salida existentes/nuevos

**Casos de prueba recomendados**
1. Caso normal:
- Archivo TSV válido con múltiples TFs
- Genoma completo de E. coli
- Verificar:
    - Número correcto de archivos FASTA
    - Integridad de secuencias extraídas

2. Coordenadas Inválidas:
    - Picos con start > end
    - Coordenadas fuera del rango del genoma
    - Verificar:
        - Omisión de picos inválidos
        - Mensajes de error apropiados
3. Formato incorrecto:
- TSV con columnas faltantes
- FASTA mal formado
- Verificar:
    - Salida con código de error
    - Mensajes descriptivos


### Diagrama de Caso de Uso (PlantUML) para Visualizar el Proceso:

```
%% Diagrama de Casos de Uso - Extracción de Secuencias de Unión a TFs
%% Basado en la implementación actual del código

graph TD
    %% Actores
    Usuario[🧑 Usuario] --> Sistema[🖥️ Sistema de Extracción]

    %% Sistema Principal
    Sistema --> M1[[Módulo 1: Procesamiento de Entradas]]
    Sistema --> M2[[Módulo 2: Extracción de Secuencias]]
    Sistema --> M3[[Módulo 3: Generación de FASTA]]

    %% Subcasos de Módulo 1
    M1 --> UC1["📂 Leer archivo de picos (peaks.py)"]
    UC1 --> UC1a["✔ Validar estructura TSV"]
    UC1 --> UC1b["🧮 Agrupar picos por TF"]
    UC1 --> UC1c["⚙️ Generar IDs únicos"]

    M1 --> UC2["🧬 Cargar genoma (genome.py)"]
    UC2 --> UC2a["🔍 Verificar formato FASTA"]
    UC2 --> UC2b["📏 Validar longitud del genoma"]

    %% Subcasos de Módulo 2
    M2 --> UC3["📍 Convertir coordenadas (extractor.py)"]
    UC3 --> UC3a["1-based → 0-based"]
    UC3 --> UC3b["🚫 Filtrar picos fuera de rango"]

    M2 --> UC4["🧬 Extraer secuencias"]
    UC4 --> UC4a["✂️ Segmentar genoma"]
    UC4 --> UC4b["📦 Almacenar en buffer"]

    %% Subcasos de Módulo 3
    M3 --> UC5["🖋️ Escribir archivos (io_utils.py)"]
    UC5 --> UC5a["🛡️ Sanitizar nombres"]
    UC5 --> UC5b["📝 Formatear encabezados"]
    UC5 --> UC5c["💾 Guardar en disco"]

    %% Flujos
    Usuario -->|"-p peaks.tsv\n-g genome.fasta\n-o output/"| M1
    M1 -->|Diccionario de picos| M2
    M2 -->|Secuencias agrupadas| M3
    M3 -->|"Archivos {TF}.fa"| Usuario

    %% Notas técnicas
    note[Notas técnicas:\n1. Coordenadas 1-based en entrada\n2. Validación de límites del genoma\n3. Formato FASTA estandarizado]
```

