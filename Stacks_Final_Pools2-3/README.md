# Stacks *de novo*

1. **process_radtags**
2. **denovo_map.pl**
3. **populations**

---

```bash
200.23.162.240
```

## 1. process_radtags

```bash
nohup process_radtags -P -p ./raw_pools -b ./barcodes/barcodes_Pool2y3.txt -o ./demultiplexed -c -q -r
-s 25 --inline_index --renz_1 ecoRI --renz_2 mspI &> process_log &
```

**Process_radtags output**

Pool 2 y 3: >99M de lecturas/pool.

| Categoría                     | Lecturas    | Porcentaje |
|-------------------------------|-------------|------------|
| Total de secuencias           | 398,933,454 | —          |
| Descartadas: barcode not found| 1,347,854   | 0.3%       |
| Descartadas: low quality read | 14,357,013  | 3.6%       |
| Descartadas: RAD cutsite not found | 2,322,453 | 0.6%    |
| Lecturas retenidas            | 380,906,134 | 95.5%      |

Lecturas Retenidas: Pool2 + Pool3

<img src="../Stacks_Prueba_1/imagenes/ecdf_reads_processRadtagsTotal.png" width="500">

*Curva de distribución acumulada de los reads retenidos por individuo post-process_radtags.* ~33 individuos obtuvieron menos de 1 millón de lecturas. 

<img src="../Stacks_Prueba_1/imagenes/barras_por_muestra.png" width="1000">
Reads retenidos por individuo. 


## Muestras excluidas (<900,000 reads)

| Barcode        | Muestra | Reads retenidos |
|----------------|---------|----------------:|
| CAACC-CGATGT   | EP13    |         589,759 |
| ACTGG-CGATGT   | EP25    |         208,992 |
| ACTTC-CGATGT   | EP26    |         351,010 |
| ATACG-CGATGT   | EP27    |         147,145 |
| ATGAG-CGATGT   | EP29    |         182,698 |
| ATTAC-CGATGT   | EP31    |         668,853 |
| CGTAC-CGATGT   | EP38    |         294,980 |
| CGTCG-CGATGT   | EP39    |         139,032 |
| CTGTC-CGATGT   | LI11    |          91,804 |
| GCCGT-CGATGT   | LI18    |         503,469 |
| GGCTC-CGATGT   | LI26    |         558,496 |
| GTAGT-CGATGT   | E4      |         413,270 |
| TACCG-CGATGT   | E13     |          38,346 |
| TCAGT-CGATGT   | SR8     |         864,141 |
| TCCGG-CGATGT   | SR10    |          46,839 |
| TCTGC-CGATGT   | SR11    |          68,910 |
| TGGAA-CGATGT   | SR13    |         140,476 |
| TTACC-CGATGT   | SR14    |         726,035 |
| ACTTC-TTAGGC   | IA17    |         557,048 |
| CATAT-TTAGGC   | IA28    |         145,897 |
| GGATA-TTAGGC   | BS6     |         103,150 |
| GTCGA-TTAGGC   | BS14    |         245,262 |
| TACCG-TTAGGC   | LOB1    |          21,999 |
| TACGT-TTAGGC   | LOB2    |         116,318 |
| TATAC-TTAGGC   | LO2     |          39,255 |
| TCACG-TTAGGC   | LO3     |          50,848 |
| TCAGT-TTAGGC   | LO4     |         490,205 |
| TCCGG-TTAGGC   | LO5     |          86,851 |
| TCTGC-TTAGGC   | LO6     |         721,876 |
| TTACC-TTAGGC   | LO10    |         365,500 |

**Total excluidas:** 30 muestras

Se excluyeron de análisis posteriores y se construyó el Popmap con 64 individuos (`popmap_1M_POPMODULE.txt`, `7 localidades`). Los individuos `E12` y `PS3`, ambos con ~900K lecturas, se ignoraron del análisis (no existen en el catálogo R1M).

---

## 2. denovo_map.pl

- Catálogo: R1M (denovo_1M_log)
- Corte 1M reads
- m5M2n4


```bash
nohup denovo_map.pl --samples ./demultiplexed --popmap ./barcodes/popmap_1M.txt -o ./stacks/R1M -m 5 -M 2 -n 4 -T 10 &> denovo_1M_log &
```

**Parámetros**
- `-m 5`: mínimo número de lecturas idénticas para formar un stack
- `-M 2`: número máximo de mismatches permitidos entre stacks de un mismo individuo
- `-n 4`: número máximo de mismatches permitidos entre individuos al construir el catálogo


Los valores de **-m**, **-M** y **-n** son modificables. No hay valores universales:

<img src="../Stacks_Prueba_1/imagenes/Denovo_parametros_ejemplos.png" width="650">


---

## 3. populations

- 7 localidades (fusión SR-E y PS-LOB-LO).
- popmap `popmap_1M_POPMODULE.txt`, `-p 5`


```bash
populations -P ./stacks/R1M --popmap ./barcodes/popmap_1M_POPMODULE.txt -O ./populations/1M/p5/ -p 5 -r 0.80 -t 5 --min-maf 0.05 --write-single-snp --genepop --vcf --fasta-loci --fasta-samples
```

**Parámetros:**
- `-r 0.8`: el locus debe estar presente en el 80% de los individuos por población (filtro r80)
- `-p 5`: número de poblaciones en las que un locus debe estar presente para conservarse en el análisis.
- `--min-maf 0.05`: filtro de frecuencia alélica menor mínima. Este umbral es estándar.
  
 **OUTPUT**:
 
- Loci retenidos: `12,757`
- Sitios variantes: `8,282`


---

## Análisis del missing data (VCFtools)


**MISSING DATA POR LOCUS**

```bash
vcftools --vcf populations.snps.vcf --missing-site --out missing_site_first01
```

```bash
sort -k6 -n -r missing_site_first01.lmiss | head -20
```

**Resultado**: muchos loci con alto missing data >20%.

Para revisar el número de loci retenidos por cada umbral, se genera un loop:

```bash
for t in 0.5 0.6 0.7 0.8 0.9 0.95 1.0; do
  n=$(awk -v t=$t 'NR>1 && (1-$6) >= t' missing_site_first01.lmiss | wc -l)
  echo -e "${t}\t${n}"
done
```

**Resultados**: a un umbral de 0.5, se retienen 8,260 loci, a 0.7 = 4,706 loci retenidos. Se seleccionó el umbral de `0.8` reteniendo `1,782`, permitiendo el 20% de missing data

```bash
vcftools --vcf populations.snps.vcf --max-missing 0.8 --recode --recode-INFO-all --out vcf_filtrado_locus_08
```

**MISSING DATA POR INDIVIDUO**

Se aplicó VCF sobre el archivo ya filtrado `vcf_filtrado_locus_08.recode.vcf`

```bash
vcftools --vcf vcf_filtrado_locus_08.recode.vcf --missing-indv --out missing_indv_post_locus_08
```

**Resultados**

| Individuo | % Missing Data |
|-----------|----------------|
| LI17 | 51% |
| LI3 | 43% |
| E14 | 40% |
| LI13 | 40% |
| LI15 | 36% |
| SR7 | 35% |
| EP5 | 33% |
| AN4 | 33% |
| EP9 | 31% |
| SR6 | 31% |
| LI24 | 31% |
| LI12 | 31% |
| LI22 | 31% |
| LI20 | 30% |
| E11 | 25% |
| EP34 | 25% |
| BS2 | 23% |
| E15 | 23% |
| EP21 | 22% |
| IA31 | 18% |


- Se eliminaron las poblaciones completas de LI y E (fusión SR + E).
- Corte: >20% de missing data.
- Número de individuos eliminados: 19


---

## 3. populations final

- 5 localidades: EP, IA, AN, PS, BS
- popmap `Popmap_5loc_p5_final.tsv`, `-p 5`
- N (individuos) = 44

```bash
populations -P ./stacks/R1M --popmap ./barcodes/Popmap_5loc_p5_final.tsv \
  -O ./populations/1M/p5_5loc_definitivo/ -p 5 -r 0.80 -t 5 --min-maf 0.05 \
  --write-single-snp --genepop --vcf --fasta-loci --fasta-samples
```

Loci retenidos: `10,161`
Total de sitios	1,484,111
Sitios variantes (SNPs) retenidos: `6,404`


*Estadísticas poblacionales*

| Población | Muestras/locus | π (pi) | Sitios (todos/variantes/polimórficos) | Alelos privados |
|-----------|----------------|--------|----------------------------------------|------------------|
| EP | 11 | 0.253 | 1,484,102 / 6,404 / 5,588 | 154 |
| PS | 8 | 0.236 | 1,484,096 / 6,404 / 5,235 | 1 |
| IA | 11 | 0.230 | 1,484,092 / 6,404 / 5,391 | 5 |
| BS | 7 | 0.231 | 1,484,081 / 6,404 / 4,700 | 2 |
| AN | 4 | 0.230 | 1,484,064 / 6,404 / 3,912 | 0 |

Verificación de missing data por individuo. Tres individuos con >20% (EP23, 24%; IA25, 21%; IA11, 20%). No se descartaron.

---

## Guideline

1. *populations*  (./populations/1M/p5/)
2. missing-site sobre populations.snps.vcf, decisión de umbral
3. --max-missing, VCF filtrado por locus
4. missing-indv sobre el VCF filtrado x locus
5. Exclusión de individuos con missingness alto y nuevo POPMAP
6. *populations* con POPMAP nuevo
7. verificación con missing-site / missing-indv sobre el output nuevo 




