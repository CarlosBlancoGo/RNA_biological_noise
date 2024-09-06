#!/bin/bash

#SBATCH --job-name=coverage_length_distribution
#SBATCH --time=1-00:00:00
#SBATCH --qos=short
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=70gb
#SBATCH --output=coverage_length_distribution.out
#SBATCH --error=coverage_length_distribution.err

# Cargar el módulo de R
module load R/4.3.3

# Definir variables de entrada y salida
input_file="/storage/gge/Carlos/aging/data_nanopore/not_filtered_nor_subsampled_bams_analysis/coverage_df_nanopore.csv"  # Reemplaza con la ruta de tu fichero
output_barplot="ONT_percent_coverage.png"
output_violinplot="ONT_lengths.png"

# Verificar si el archivo existe antes de ejecutar el código R
if [[ ! -f "$input_file" ]]; then
  echo "Error: El archivo $input_file no existe."
  exit 1
fi

# Llamar a R para procesar el fichero y generar gráficos
Rscript --vanilla -e "
# Librerías necesarias
library(ggplot2)
library(dplyr)

# Leer el fichero CSV
data <- read.csv('$input_file', header = TRUE)

# Comprobar si las columnas necesarias existen
if(!all(c('percent_coverage', 'length', 'sample') %in% colnames(data))) {
    stop('Las columnas percent_coverage, length o sample no están presentes en el fichero.')
}

samples_order <- c('B31', 'B32', 'B33', 'B34', 'B35', 'B151', 'B152', 'B153', 'B154')
data\$group <- ifelse(data\$sample %in% c('B151', 'B152', 'B153', 'B154'), 'old', 'young')
data\$sample <- factor(data\$sample, levels = samples_order)
data\$group <- factor(data\$group, levels = c('young', 'old'))

# Crear el barplot de percent_coverage por muestra
plot_coverage <- ggplot(data, aes(x = sample, y = percent_coverage, fill = group, color = group)) +
  geom_violin(alpha = 0.3) +
  geom_boxplot(alpha = 0.5, width = 0.5, outlier.shape = NA) + 
  scale_fill_manual(values = c('young' = '#15918A', 'old' = '#FDC659')) +
  scale_color_manual(values = c('young' = '#15918A', 'old' = '#FDC659')) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 20),
    axis.text.y = element_text(size = 20),
    axis.title.x = element_text(size = 25),
    axis.title.y = element_text(size = 25),
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 18),
    legend.key.size = unit(0.5, 'cm'),
    panel.background = element_rect(fill = 'transparent', colour = NA), 
    plot.background = element_rect(fill = 'transparent', colour = NA),
    legend.background = element_rect(fill = 'transparent', colour = NA),
    plot.margin = margin(t = 10, r = 10, b = 10, l = 10)) +
  labs(x = 'Sample', y = '% of coverage')

# Guardar el barplot
ggsave('$output_barplot', plot = plot_coverage, width = 10, height = 7)

# Crear el violin plot de length por muestra
plot_length <- ggplot(data, aes(x = sample, y = length, fill = group, color = group)) +
  geom_violin(alpha = 0.3) +
  geom_boxplot(alpha = 0.5, width = 0.5, outlier.shape = NA) + 
  scale_fill_manual(values = c('young' = '#15918A', 'old' = '#FDC659')) +
  scale_color_manual(values = c('young' = '#15918A', 'old' = '#FDC659')) +
  coord_cartesian(ylim = c(0, 10000)) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 20),
    axis.text.y = element_text(size = 20),
    axis.title.x = element_text(size = 25),
    axis.title.y = element_text(size = 25),
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 18),
    legend.key.size = unit(0.5, 'cm'),
    panel.background = element_rect(fill = 'transparent', colour = NA), 
    plot.background = element_rect(fill = 'transparent', colour = NA),
    legend.background = element_rect(fill = 'transparent', colour = NA),
    plot.margin = margin(t = 10, r = 10, b = 10, l = 10)) +
  labs(x = 'Sample', y = 'Length')

# Guardar el violin plot
ggsave('$output_violinplot', plot = plot_length, width = 10, height = 7)
"

echo "Gráficos generados y guardados como PNG."
