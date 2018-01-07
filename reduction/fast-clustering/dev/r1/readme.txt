run_multicore.py

Es el script orquestador.

Parámetros:
- inputDir = Directorio con las proteínas
- outputDir = Directorio de salida
- nCpus = número de procesadores disponibles para el cálculo

Se encarga de dividir el trabajo para cada procesador. Además crea un directorio con enlaces simbólicos para los archivos de entrada, creando subdirectorios sobre los que aplicará el clustering.
Creará en outputDir tantos subdirectorios como procesadores se hayan indicado en nCpus. 


apply_cluster_by_blocks

Parámetros:
- sub_dir: Subdirectorio con las secuencias que le corresponden a cada procesador.
- maxFiles: Tamaño del bloque que se podrá procesar. Corresponde al limite para que el sistema no colapse por memoria.

Acá se establece el porcentaje (de 0 a 1) que se desea utilizar para la reducción


cluster.r
La primera línea debe modificarse para indicar la ruta donde R se encuentra instalado… por ahora hay que hacer este cambio manual


Muscle debe ser instalado. El error:
You do not have muscle installed/working locally on your machine