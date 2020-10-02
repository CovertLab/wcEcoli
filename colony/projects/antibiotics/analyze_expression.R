library(tidyverse)
library(readxl)
library(ggridges)

TANIGUCHI_S6_PATH = file.path("assets", "TableS6.xls")
GENE = "yaaA"

# Plot Reference Data

taniguchi_s6 <- read_excel(TANIGUCHI_S6_PATH)
gamma_parameters <- taniguchi_s6[, c("Gene Name", "A_Protein",
                                     "B_Protein")]
gene_parameters_map <- gamma_parameters$"Gene Name" == GENE
a <- gamma_parameters[gene_parameters_map, "A_Protein"][[1]]
b <- gamma_parameters[gene_parameters_map, "B_Protein"][[1]]
taniguchi_plot <- ggplot(data=tibble(x = c(0, 10)), aes(x = x)) +
    stat_function(fun=dgamma, args=list(shape=a, rate=b)) +
    labs(title = "Distribution of Counts per Cell for yaaA",
         subtitle = "Uses Parameters from Table S6, Taniguchi 2010") + 
    xlab("Counts per Cell") + ylab("Frequency")
ggsave("taniguchi_plot.pdf", taniguchi_plot)

# Plot Expression Data
data <- read_csv("expression.csv")
protein <- c()
expression <- c()
for (value in data$"EG10040-MONOMER[p]"){
    protein <- c(protein, "Beta-Lactamase")
    expression <- c(expression, value)
}
for (value in data$"TRANS-CPLX-201[s]"){
    protein <- c(protein, "AcrAB-TolC")
    expression <- c(expression, value)
}
transformed <- tibble(protein=protein, expression=expression)
expression_plot <- ggplot(transformed, aes(x=expression, y=protein)) +
    geom_density_ridges(
        jittered_points = TRUE,
        position = position_points_jitter(width = 0.05, height = 0),
        point_shape = '|', point_size = 3, point_alpha = 1,
        alpha = 0.7,
    ) +
    labs(title = "Protein Concentration Distributions",
         subtitle = "From Experiment 20200924.162551") + 
    xlab("Protein Concentration (counts/fL)") + ylab("Protein")

ggsave("expression_plot.pdf", expression_plot)
