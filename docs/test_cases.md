### Caso de prueba 1: [Validación de coordenadas inválidas]

- Descripción: Verificar que el sistema registre un error o notifique al usuario cuando los valores de Peak_start y Peak_end no se encuentren dentro del rango permitido del genoma.
- Datos de entrada: 
    - Archivo de picos con un registro que tenga Peak_start o Peak_end mayores que la longitud del genoma o negativos.
    - Archivo FASTA del genoma.
- Resultado esperado: 
    - El sistema debe detectar la inconsistencia en las coordenadas.
    - Se debe registrar un mensaje de error (o mostrarlo en pantalla) indicando que las coordenadas son inválidas.
    - El proceso debe omitir el registro con coordenadas erróneas y continuar con los registros válidos, o detener la ejecución según la configuración.

### Caso de prueba 2: [Detección de archivos inexistentes o ilegibles]

- Descripción: Verificar que el sistema notifique apropiadamente si uno o ambos archivos de entrada (archivo de picos o FASTA) no existen o no se pueden leer.
- Datos de entrada: 
    - Probar con rutas de archivo incorrectas o archivos con permisos restringidos.
- Resultado esperado: 
    - El sistema debe mostrar un mensaje de error indicando la ausencia o imposibilidad de lectura, y detener la ejecución o gestionar el error según lo especificado.

### Caso de prueba 3: [Generación correcta del script de automatización y logs]

- Descripción: Verificar que el script run_meme.sh se genere correctamente, contenga las líneas de comando necesarias para ejecutar meme en cada archivo FASTA, que sea ejecutable y que se genere un log con la información de la ejecución.
- Datos de entrada: 
    - Un directorio con varios archivos FASTA válidos
- Resultado esperado: 
    - Se genera el archivo run_meme.sh en el directorio de trabajo.
    - El script contiene una línea de comando para cada archivo FASTA, con los parámetros predefinidos.
    - Se genera un archivo de log.
    - El archivo es ejecutable y, al ejecutarlo, se inicia el análisis con meme en cada archivo.