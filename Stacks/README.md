
# Pipeline bioinformático con Stacks para *Octopus mimus*

Este es un pipeline para el procesamiento de datos ddRADseq de *Octopus mimus*, desde las lecturas crudas hasta las estadísticas poblacionales, usando Stacks. 

## Resumen del pipeline

1. **process_radtags**: demultiplexa y limpia las lecturas crudas
2. **denovo_map.pl**: ensambla de novo y genera compendios de loci
3. **populations**: estadísticas poblacionales y exportación de datos

<img src="imagenes/Stax_pipeline.png" width="700">
Stacks pipeline (Catchen et al., 2013)

---

## ¿Cómo entro al servidor?

En este caso al server del CIBNOR
1. Abro Hillstone secure connect para conectarme con VPN (si no me encuentro presencialmente en el CIBNOR) y en mi caso hago auto-connect.

<img src="imagenes/vpn_connect.png" width="200">

2. Me conecto al server usando Putty, escribo la IP del server (200.23.162.240), puerto 22.
3. Escribo mi usuario y contraseña.

<img src="imagenes/Server_connect.png" width="300">

## Almacenamiento de las lecturas crudas

Dentro del servidor (CIBNOR) las lecturas crudas se deben almacenar en la carpeta de Datos, entonces tanto el raw_data (cada pool) y el barcodefile, así como el output con las lecturas demultiplexadas y limpias, se almacenan en Datos (/Datos/user/raw_pools).

<img src="imagenes/raw_pools.png" width="400">

## Formatos de los archivos: barcodes & popmap

**BARCODEFILE:**
La estructura del barcodefile es **BARCODE** | **INDEX** | **SAMPLE_NAME**

Debe ser un archivo sin extensión o de texto o separado por tabulaciones.

<img src="imagenes/barcodefile_format.png" width="500">

**Notas a considerar**

☙ — Los barcodes e indexes que se utilizan durante la library prep provienen del manual de ddRADseq de Peterson *et al.* (2012). En el Laboratorio de Molecular y Genética del Cibnor estos NO se modificaron, se utilizan siempre los originales. 

En mi caso, para los pools 2 y 3 (individuos de Ecuador y Perú), utilicé los indexes 2 (CGATGT) y 3 (TTAGGC).

<img src="imagenes/PCRindex_primers.png" width="400">

☙ — Se tiene que checar que todas las carpetas que se ocupan como datos input deben estar localizadas en el mismo escalón o rama de árbol de directorios.

<img src="imagenes/InputDataLocation.png" width="500">

## 1. process_radtags

Demultiplexa las lecturas por individuo a partir de los archivos crudos 
del secuenciador Novogene y elimina lecturas de baja calidad.

Primeramente se llama al ambiente de stacks.

```bash
conda activate stacks
```

```bash
nohup process_radtags -P -p ./raw_pools -b ./barcodes/barcodes_Pool2y3.txt -o ./demultiplexed -c -q -r
-s 25 --inline_index --renz_1 ecoRI --renz_2 mspI &> process_log &
```

**Parámetros:**
- `-P`: datos paired-end (lecturas R1 + R2 emparejadas)
- `-p -/raw_pools`: la ruta hacia la carpeta donde están mis lecturas crudas del pool
- `-b barcodes_Pool2y3.txt`: ruta hacia la carpeta de archivo de barcodes
- `-o ./demultiplexed`: carpeta de salida para los archivos FastQ separados
- `-c`: clean: remueve aquellas lecturas con bases indeterminadas (N)
- `-q`: filtro  de calidad (Phred)
- `-r`: rescata barcodes/sitios de corte con errores menores
- `-s 25`: umbral de Phred para el parámetro de -q
- `--inline_index`: índice dentro de la lectura, no como archivo aparte
- `--renz_1 ecoRI`: enzima corte raro (sitio GAATTC)
- `--renz_2 mspI`: enzima corte frecuente (CCGG)

**Notas**
- `nohup`: permite dejar el análisis en el background. Si uno se va a tomar un café y cierra la PC, el server sigue computando el run.
- `process_log`: genera las notificaciones del análisis.

Para checar los procesos actuales del servidor que están corriendo:

```bash
top
```

o

```bash
htop
```

<img src="imagenes/topcommand.png" width="600">

Para observar si *process_radtags* se congeló y no computó, o si simplemente se quiere saber el pasó en el que va o, en su caso, el desenlace de la corrida, entonces se utiliza el comando:

```bash
more process_log
```

<img src="imagenes/moreProcess_radtags_log.png" width="600">

**Consideraciones (OJO)**

Los nombres de los archivos puede ser que no sean los adecuados para que Stacks los ejecute. Los archivos crudos deben tener la extensión de .fastq.gz y NO .fq.gz (como se muestra en la siguiente imagen).

<img src="imagenes/ExtensionFastq.png" width="600">

**Process_radtags output**

Al finalizar el run, se obtuvieron más de 99M de lecturas para cada pool. El tiempo de computo fue de aproximadamente 1 hora x pool. A continuación se muestran los archivos demultiplexados del pool 3:

<img src="imagenes/demultiplexed_pool3.png" width="600">

Una vez concluido el *process_radtags*, se obtuvieron los números siguientes:

<img src="imagenes/process_radtags_finished.png" width="500">

---

## 2. denovo_map.pl

Ensambla los loci de novo (sin genoma de referencia) y construye el 
catálogo compartido entre individuos.



En un inicio se llevó a cabo tres análisis exploratorios de *denovo_map*. Con base en las lecturas por individuo obtenidas durante el demultiplexado (en la etapa de *process_radtags*), se realizaron tres cortes de filtrado de individuos, ya que el número de lecturas entre las 96 muestras fue heterogéneo. Los cortes fueron para 1) retener aquellos individuos que presentaron igual o mayor a un millón de lecturas (1M), 2) retener aquellos con igual o mayor a 750 millones de lecturas (750K) y 3) igual o mayor a 500 millones (500K).

**1M READS**
```bash
nohup denovo_map.pl --samples ./demultiplexed --popmap ./barcodes/popmap_1M.txt -o ./stacks/R1M -m 5 -M 2 -n 4 -T 10 &> denovo_1M_log &
```

**750K READS**
```bash
nohup denovo_map.pl --samples ./demultiplexed --popmap ./barcodes/popmap_750k.txt -o ./stacks/R750K -m 5 -M 2 -n 4 -T 20 &> denovo_750K_log &
```

**500K READS**
```bash
nohup denovo_map.pl --samples ./demultiplexed --popmap ./barcodes/popmap_500k.txt -o ./stacks/R500K -m 5 -M 2 -n 4 -T 20 &> denovo_500K_log &
```

**Parámetros clave:**
- `-m 5`: mínimo número de lecturas idénticas para formar un stack
- `-M 2`: número máximo de mismatches permitidos entre stacks de un mismo individuo
- `-n 4`: número máximo de mismatches permitidos entre individuos al construir el catálogo

**Notas**

En una primera instancia, los tres cortes se corrieron al mismo tiempo. Sin embargo, los runs de 1M y 750K abortaron el análisis por falta de memoria del server. 750K si generó la corrida completa. Para evitar que las tres cortes compitieran por recursos de memoria, se corrió una por una. Correr *denovo.map* tardó en mi caso aproximadamente >9 horas x corrida, esto depende en parte del número de muchos factores, como el número de reads por individuo, la parametrización que se utilice, el número de threads. Para los análisis de 750K y 1M, utilicé 20 threads porque ningún otro usuario del servidor estaba utilizándolo. Importante revisar antes de enviar cada run.

**Output**

EJEMPLO. Se generaron los archivos resultantes del run de 1M reads (/Datos/smunguia/stacks/R1M):

<img src="imagenes/Output_1Mreads_archivos.png" width="600">

Revisé el archivo *population.log*, el cual compila diferentes estadísticas, como el número total de loci retenidos, el número de sitios variantes; y algunas estadísticas poblacionales por localidad (e.g. diversidad nucleotídica, sitios polimórficos, alelos privados, etc.).

<img src="imagenes/PopLOG_Locality_R1M.png" width="600">

Asimismo, para las otras dos cortes de 500K y 750K reads, se generaron los archivos output correspondientes (no se muestran en este readme).

---

###**Segunda aproximación: correr *denovo_map.pl* sin filtros de número de lecturas y en el módulo de *populations* ir depurando la base de datos**

Para esto, se generó un nuevo archivo Popmap incluyendo a todos los individuos; los valores de la parametrización no se modificaron (-m, -M, -n). Se ejecutó el siguiente comando:

```bash
nohup denovo_map.pl --samples ./demultiplexed --popmap ./barcodes/PopMap_sinFiltros.tsv -o ./stacks/RAllinOne -m 5 -M 2 -n 4 -T 10 &> denovo_all_log &
```

<img src="imagenes/RAllinOne_m5M2n4_files.png" width="600">



**Consideraciones**

Posteriormente, cuando obtenga las secuencias del total de pools para *O. mimus* y con el objetivo de optimizar el ensamblaje de novo, evaluaré sistemáticamente diferentes combinaciones de parámetro usando RADstackshelpR (DeRaad, 2021; https://github.com/DevonDeRaad/RADstackshelpR). Por ahora me quedo con la parametrización estándar.

Es importante mencionar que los valores de **-m**, **-M** y **-n**, son modificables respecto a los datos obtenidos de la secuenciacion y su procesamiento con el pipeline de stacks. Por eso es que para diferentes especies, estos números son distintos:

<img src="imagenes/Denovo_parametros_ejemplos.png" width="600">


Corrida unificada "RAllinOne" — catálogo construido con las 10 
localidades combinadas.

![Resumen del catálogo](imagenes/denovo_map_summary.png)

---

## 3. populations

Genera estadísticas poblacionales, filtra loci según parámetros de 
representación, y exporta formatos de salida (VCF, structure, etc.).

```bash
populations -P ./stacks_output/ -M ./popmap_RAllinOne.tsv \
  -r 0.8 --vcf --genepop -t 8
```

**Notas:**
- `-r 0.8`: el locus debe estar presente en el 80% de los individuos 
  por población (filtro r80)
- Salidas usadas después para DAPC (adegenet) y análisis de missing 
  data (VCFtools)

![Estadísticas de populations](imagenes/populations_stats.png)

---

## Problemas conocidos / decisiones tomadas

- Individuos **LI17** y **LI3** mostraron >40% de missing data — 
  evaluados para exclusión.
- Localidades de bajo n (Los Órganos, Lobitos) — en evaluación para 
  fusionar con localidades vecinas antes de re-correr denovo_map.
- RAM: evitar correr cstacks en simultáneo con múltiples hilos si el 
  servidor tiene memoria limitada (crashes observados).

## Referencias

- Paris, J.R., Stevens, J.R., Catchen, J.M. (2017). Lost in parameter 
  space: a road map for stacks. *Methods in Ecology and Evolution*.
