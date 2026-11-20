## Optimización de parámetros con RadStackshelpR

[RADstackshelpR](https://devonderaad.github.io/RADstackshelpR/index.html) que se basa en [Paris et al. (2017)](https://doi.org/10.1111/2041-210X.12775): al shufflear valores de `m`, `M` & `n`, maximizand el número de loci polimórficos presentes en ≥ 80% de los indv.

**Popmap + parámetro:** `RF53_m{m}M{M}n{n}` (RF53: popmap final de 53 individuos). *esto fue post-eliminación de muestras ruidosas*

---

## Primer run con las 96 muestras (`PopMap_sinFiltros.tsv`) *R80 en RadstackshelpR descartó la mayoría porque no hubo un consenso*

- `denovo_map.pl` con m = 3-6, M = 2, n = 1, popmap con los 96 individuos

*looping everything from now on*:

```bash
for m in 3 4 5 6; do d=./stacks/RAllinOne_m${m}M2n1; mkdir -p $d; nohup nice -n 10 denovo_map.pl --samples ./demultiplexed --popmap ./barcodes/PopMap_sinFiltros.tsv -o $d -m $m -M 2 -n 1 -T 15 &> ./denovo_all_log_m${m}M2n1 & done
```

`populations` **sin ningún filtro** (el VCF que espera RADstackshelpR), desde `stacks/`:

```bash
for m in 3 4 5 6; do (cd RAllinOne_m${m}M2n1 && nohup populations -P ./ --popmap ../../barcodes/PopMap_sinFiltros.tsv --vcf -t 8 &> ../populations_log_m${m}M2n1 &); done
```

Descarga desde Cib-Server, renombrando x corrida:

```bash
for m in 3 4 5 6; do rsync -ah --progress --partial cib-server:/Datos/meow/stacks/RAllinOne_m${m}M2n1/populations.snps.vcf ~/Documents/Populations/${m}21.vcf; done
```

Pues primeramente en R se corrió la función de optimize que se basa en el R80 y en mi caso, colapsó
```r
library(RADstackshelpR)
m.out <- optimize_m(
  m3 = "~/Documents/Populations/321.vcf", m4 = "~/Documents/Populations/421.vcf",
  m5 = "~/Documents/Populations/521.vcf", m6 = "~/Documents/Populations/621.vcf")
vis_depth(output = m.out); vis_snps(output = m.out, stacks_param = "m"); vis_loci(output = m.out, stacks_param = "m")
m.out$loci.R80; m.out$snp.R80
```

**Resultado (96 muestras):** el R80 colapsa. =(

| m | Loci R80 | SNPs R80 |
|---|---|---|
| 3 | 9 | 28 |
| 4 | 6 | 17 |
| 5 | 6 | 20 |
| 6 | 7 | 18 |

La profundidad por muestra era adecuada (~20x), pero las distribuciones de SNPs y loci por muestra eran **bimodales** alcanzando almost los 0 SNPs retenidos.

**Por qué?:** `optimize_m` conserva un SNP si le faltan como máximo el 20% de los genotipos (`rowSums(is.na(gt)) / ncol(gt) <= 0.2`). Las muestras con pocos datos aportan casi solo faltantes y hunden el conteo. La documentación del paquete describe este caso y recomienda eliminar las muestras de pocos datos y repetir Stacks.

---

## Para elimininar muestras ruidosas, diagnosticamos qué pasa

**Loci por muestra en `ustacks`**:

```bash
cd /Datos/meow/stacks/RAllinOne_m3M2n1
for f in *.tags.tsv.gz; do echo -e "${f%.tags.tsv.gz}\t$(zcat $f | grep -c consensus)"; done | sort -k2,2n | head -40
```

Las más pobres iban de 264 loci (SR10) a ~47,000

**SNPs por muestra desde el VCF** (R):

```r
n.snp <- colSums(!is.na(gt))
pop <- read.table("~/Documents/Populations/PopMap_sinFiltros.tsv", sep = "\t", col.names = c("sample", "pop"))
tab <- merge(data.frame(sample = names(n.snp), n.snp), pop, by = "sample")
tab$rank <- rank(-tab$n.snp)
tab <- tab[order(tab$rank), ]
```

**Curva de R80 añadiendo muestras de mejor a peor** (m3):

```r
ord <- names(sort(colSums(!is.na(gt)), decreasing = TRUE))
res <- data.frame()
for (k in seq(24, 96, by = 4)) {
  gk <- gt[, ord[1:k]]
  pass <- rowSums(is.na(gk)) / ncol(gk) <= 0.2
  res <- rbind(res, data.frame(k = k, snps = sum(pass), loci = length(unique(vcf.r@fix[pass, 1]))))
}
plot(res$k, res$loci, type = "b", xlab = "muestras incluidas", ylab = "loci R80")
```

| Muestras | Loci R80 |
|---|---|
| 24 | 55,871 |
| 48 | 40,288 |
| 60 | 27,841 |
| 72 | 10,181 |
| 96 | 9 |

Es un gradiente: caída fuerte hacia 60-64 muestras, colapso después de ~80.

**Aproximación rápida (recortando el VCF de 96 a las 48 de `Popmap_6loc_p6.tsv`)**. Esto ayudó para orientarme porque el catálogo se construyó con las 96:

```r
keep <- read.table("~/Documents/Populations/Popmap_6loc_p6.tsv", sep = "\t")[, 1]
for (m in 3:6) {
  v <- vcfR::read.vcfR(sprintf("~/Documents/Populations/%d21.vcf", m), verbose = FALSE)
  g <- v@gt[, -1]; g <- g[, colnames(g) %in% keep]
  pass <- rowSums(is.na(g)) / ncol(g) <= 0.2
  print(c(m = m, snps = sum(pass), loci = length(unique(v@fix[pass, 1]))))
  rm(v, g, pass); gc()
}
```

| m | Muestras | SNPs R80 | Loci R80 |
|---|---|---|---|
| 3 | 48 | 227,635 | 38,625 |
| 4 | 48 | 206,335 | 34,871 |
| 5 | 48 | 182,674 | 30,932 |
| 6 | 48 | 158,620 | 27,076 |

**Controles de que elegir las 48 con la corrida m5M2n4 no sesga la comparación:**

```r
sn <- function(f) { v <- vcfR::read.vcfR(f, verbose = FALSE); colSums(!is.na(v@gt[, -1])) }
cor(sn("~/Documents/Populations/321.vcf"), sn("~/Documents/Populations/621.vcf"), method = "spearman")  # 0.984
length(intersect(ord[1:48], keep))   # 45 de las 48 mejores por SNPs
```
*hice un check de que se alinearon los chakras bien*, y entonces se probó lo siguiente:

---

## Selección del popmap final

Recuperación de muestras de calidad intermedia. Popmaps ampliados por rango de SNPs (a diferentes cortes: ≤ 56, 60 y 66), con las fusiones previas de SR con E y LO/LOB con PS:

```r
pm48 <- read.table("~/Documents/Populations/Popmap_6loc_p6.tsv", sep = "\t", col.names = c("sample", "pop"))
orig <- read.table("~/Documents/Populations/PopMap_sinFiltros.tsv", sep = "\t", col.names = c("sample", "pop"))
orig$pop <- ifelse(orig$pop == "SR", "E", ifelse(orig$pop %in% c("LO", "LOB"), "PS", orig$pop))
for (k in c(56, 60, 66)) {
  extra <- setdiff(tab$sample[tab$rank <= k], pm48$sample)
  out <- rbind(pm48, orig[orig$sample %in% extra, ])
  write.table(out, sprintf("~/Documents/Populations/popmap_top%d_fus.tsv", k),
              sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
}
```

```bash
scp ~/Documents/Populations/popmap_top*_fus.tsv cib-server:/Datos/meow/barcodes/
cd /Datos/meow/stacks
for k in 56 60 66; do mkdir -p test_top${k}_fus; nohup populations -P RAllinOne_m3M2n1 --popmap ../barcodes/popmap_top${k}_fus.tsv -O test_top${k}_fus -r 0.8 -p 7 -t 4 &> test_top${k}_fus/stdout.log & done
mkdir -p test_p48; nohup populations -P RAllinOne_m3M2n1 --popmap ../barcodes/Popmap_6loc_p6.tsv -O test_p48 -r 0.8 -p 6 -t 4 &> test_p48/stdout.log &
grep -H "^Kept" test_p48/populations.log test_top*_fus/populations.log
```

| Popmap | Poblaciones / `-p` | Loci conservados |
|---|---|---|
| 48 (control) | 6 / `-p 6` | 11,829 |
| top56_fus | 7 / `-p 7` | 6,992 |
| top60_fus | 7 / `-p 7` | 3,560 |
| top66_fus | 7 / `-p 7` | 2,775 |

E+SR (3-5 individuos) hundía los loci: `-r 0.8` exige que el locus esté en casi todos. Se repitió sin E, con 6 poblaciones y `-p 6`:

```r
for (k in c(56, 60, 66)) {
  f <- read.table(sprintf("~/Documents/Populations/popmap_top%d_fus.tsv", k), sep = "\t", col.names = c("sample", "pop"))
  f <- f[f$pop != "E", ]
  write.table(f, sprintf("~/Documents/Populations/popmap_top%d_sinE.tsv", k), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
}
```

| Popmap (6 poblaciones, `-p 6`) | Muestras | Loci | vs. control |
|---|---|---|---|
| 48 (control) | 48 | 11,829 | |
| top56_sinE | 53 | 14,324 | +21 % |
| top60_sinE | 56 | 12,147 | +3 % |
| top66_sinE | 61 | 4,737 | −60 % |

**Popmap final:** `popmap_top56_sinE.tsv`, **N = 53** (AN 5, BS 8, EP 12, IA 13, LI 5, PS 10). E+SR queda fuera de la optimización de parámetros

> **Nota:** respecto al popmap de 48, este incorpora 5 individuos (IA31, BS2, LI15, AN4 y PS3), elegidos por SNPs genotipados por muestra en el VCF sin filtrar y para dar más individuos a AN y LI. El `vcftools --missing-indv` del ensamblaje final debe confirmar si se quedan.

---

## Repetición del ensamblaje con las 53 muestras: m = 3-6

```bash
for m in 3 4 5 6; do d=./stacks/RF53_m${m}M2n1; mkdir -p $d; nohup nice -n 10 denovo_map.pl --samples ./demultiplexed --popmap ./barcodes/popmap_top56_sinE.tsv -o $d -m $m -M 2 -n 1 -T 15 &> ./denovo_RF53_m${m}M2n1 & done
```

Comprobar que terminaron (`denovo_map.pl is done.`):

```bash
pgrep -af denovo_map
tail -n 2 denovo_RF53_m*M2n1
grep -iE "error|killed|aborted|segmentation|failed" denovo_RF53_m*M2n1
```

`populations` sin filtros:

```bash
cd /Datos/meow/stacks
for m in 3 4 5 6; do (cd RF53_m${m}M2n1 && nohup populations -P ./ --popmap ../../barcodes/popmap_top56_sinE.tsv --vcf -t 8 &> ../populations_RF53_m${m}M2n1 &); done
```
Descargo VCFs:

```bash
for m in 3 4 5 6; do rsync -ah --progress --partial cib-server:/Datos/meow/stacks/RF53_m${m}M2n1/populations.snps.vcf ~/Documents/Populations/RF53_${m}21.vcf; done
```
Corro RadstackshelpR:
```r
m.out <- optimize_m(
  m3 = "~/Documents/Populations/RF53_321.vcf", m4 = "~/Documents/Populations/RF53_421.vcf",
  m5 = "~/Documents/Populations/RF53_521.vcf", m6 = "~/Documents/Populations/RF53_621.vcf")
m.out$loci.R80; m.out$snp.R80
```

| m | Loci R80 | Cambio | SNPs R80 | Cambio |
|---|---|---|---|---|
| **3** | **32,745** | | 177,277 | |
| 4 | 28,542 | −12.8 % | 154,959 | −12.6 % |
| 5 | 24,250 | −15.0 % | 131,479 | −15.2 % |
| 6 | 19,946 | −17.7 % | 107,495 | −18.2 % |

Con 53 muestras, un SNP pasa el R80 si está en ≥ 43. **`m = 3` es el óptimo**.

---

## Optimización de `M` (m = 3 fijo, n = M)

```bash
for M in 1 2 3 4 5 6 7 8; do d=./stacks/RF53_m3_M${M}n${M}; mkdir -p $d; nohup nice -n 10 denovo_map.pl --samples ./demultiplexed --popmap ./barcodes/popmap_top56_sinE.tsv -o $d -m 3 -M $M -n $M -T 7 &> ./denovo_RF53_m3_M${M}n${M} & done
```

```bash
cd /Datos/smunguia
pgrep -af denovo_map | grep -oE "RF53_m3_M[0-9]n[0-9]"
for M in 1 2 3 4 5 6 7 8; do echo "M$M: $(tail -n 1 denovo_RF53_m3_M${M}n${M} | tr -d '\n')"; done
cd stacks
for M in 1 2 3 4 5 6 7 8; do (cd RF53_m3_M${M}n${M} && nohup populations -P ./ --popmap ../../barcodes/popmap_top56_sinE.tsv --vcf -t 7 &> ../populations_RF53_m3_M${M}n${M} &); done
grep -H "^Kept" RF53_m3_M*n*/populations.log
```

loci de 1,168,787 (M1) a 869,688 (M8); sitios variantes de 577,393 a 1,075,422.

```bash
for M in 1 2 3 4 5 6 7 8; do rsync -ah --progress --partial cib-server:/Datos/meow/stacks/RF53_m3_M${M}n${M}/populations.snps.vcf ~/Documents/Populations/RF53_M${M}.vcf; done
```

```r
d <- "~/Documents/Populations/"
M.out <- optimize_bigM(
  M1 = paste0(d, "RF53_M1.vcf"), M2 = paste0(d, "RF53_M2.vcf"), M3 = paste0(d, "RF53_M3.vcf"), M4 = paste0(d, "RF53_M4.vcf"),
  M5 = paste0(d, "RF53_M5.vcf"), M6 = paste0(d, "RF53_M6.vcf"), M7 = paste0(d, "RF53_M7.vcf"), M8 = paste0(d, "RF53_M8.vcf"))
M.out$loci.R80; M.out$snp.R80
vis_snps(output = M.out, stacks_param = "M"); vis_loci(output = M.out, stacks_param = "M")
```

| M | Loci R80 | Cambio | SNPs R80 |
|---|---|---|---|
| 1 | 30,850 | | 151,074 |
| 2 | 34,862 | +13.0 % | 202,655 |
| 3 | 36,126 | +3.6 % | 222,001 |
| 4 | 36,440 | +0.9 % | 229,757 |
| **5** | **36,552** | +0.3 % | 234,324 |
| 6 | 36,488 | −0.2 % | 235,893 |
| 7 | 36,412 | −0.2 % | 237,004 |
| 8 | 36,301 | −0.3 % | 237,620 |

**`M = 5`** es el máximo y el valor que marca el paquete, pero entre M = 4 y 8 la diferencia es de ~0.3 %. De M1 a M3 la mejora es real y a partir de M5 los loci bajan levemente. 

---

## Optimización de `n` (m = 3, M = 5) 

`n` = M−1, M, M+1 = 4, 5, 6. `n = 5` ya existe (`RF53_m3_M5n5`):

```bash
cd /Datos/meow
for n in 4 6; do d=./stacks/RF53_m3_M5n${n}; mkdir -p $d; nohup nice -n 10 denovo_map.pl --samples ./demultiplexed --popmap ./barcodes/popmap_top56_sinE.tsv -o $d -m 3 -M 5 -n $n -T 15 &> ./denovo_RF53_m3_M5n${n} & done
```

```bash
cd /Datos/meow/stacks
for n in 4 6; do (cd RF53_m3_M5n${n} && nohup populations -P ./ --popmap ../../barcodes/popmap_top56_sinE.tsv --vcf -t 8 &> ../populations_RF53_m3_M5n${n} &); done
```

```bash
for n in 4 6; do rsync -ah --progress --partial cib-server:/Datos/meow/stacks/RF53_m3_M5n${n}/populations.snps.vcf ~/Documents/Populations/RF53_M5n${n}.vcf; done
```

| n   | Poli-loci R80 |
|:----|--------------:|
| M-1 | 36,402        |
| =M  | 36,552        |
| M+1 | 36,638        |

## Final grid with `m` | `M` | `n` 

![Optimización de parámetros](Stacks_Prueba_1/imagenes/Radstackshelp_resultP2-P3.png)


---

## Guideline: resumen

1. `denovo_map.pl` con las 96, se corrió `populations` sin filtros y con `optimize_m` el R80 colapsa (9 loci)
2. Diagnóstico por muestra: un bloque de muestras recupera muy pocos loci.
3. Elección del popmap: `popmap_top56_sinE.tsv` (53 individuos, 6 poblaciones).
4. Repetición del `denovo_map.pl` con el popmap final, m = 3-6 y `optimize_m` seleccinó **m = 3**.
5. `denovo_map.pl` con M = 1-8 (n = M). El `optimize_bigM` eligió **M = 5**.
6. `n` = 4, 5, 6: `optimize_n` (en curso, van >6 horas).

