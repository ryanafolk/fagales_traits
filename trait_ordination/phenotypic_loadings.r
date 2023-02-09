pheno = read.table("coded_oneval_out.csv", header = TRUE, sep = ",", stringsAsFactors=FALSE)
pheno$X <- NULL


trait = read.csv("trait_phylomds_PC1.txt", sep = "\t", header = FALSE)
colnames(trait) <- c("TAXON", "trait.vector")

# Match matrix taxa to trait.vector taxa
pheno <- pheno[ as.character(pheno$TAXON) %in% trait$TAXON, ]

pheno$Habitat <- NULL
pheno$Elevation <- NULL
pheno$Name.source <- NULL
pheno$Special.characters <- NULL
pheno$Notes <- NULL
pheno$High.missing <- NULL
pheno$Presence.. <- NULL

pheno = merge(pheno, trait, by = "TAXON")

# Treat categorical data as factors for proper R^2
# as.numeric removes any non-numeric items
pheno$Woody <- as.factor(as.numeric(pheno$Woody))
pheno$Perennial <- as.factor(as.numeric(pheno$Perennial))
pheno$Leaf.shape <- as.factor(as.numeric(pheno$Leaf.shape))
pheno$Inflorescence.type <- as.factor(as.numeric(pheno$Inflorescence.type))
pheno$Sepal.number <- as.factor(as.numeric(pheno$Sepal.number))
pheno$Sepal.color <- as.factor(as.numeric(pheno$Sepal.color))
pheno$Sepal.shape <- as.factor(as.numeric(pheno$Sepal.shape))
pheno$Petal.number <- as.factor(as.numeric(pheno$Petal.number))
pheno$Petal.color <- as.factor(as.numeric(pheno$Petal.color))
pheno$Petal.shape <- as.factor(as.numeric(pheno$Petal.shape))
pheno$Stamen.number <- as.factor(as.numeric(pheno$Stamen.number))
pheno$Flowering.time <- as.factor(as.numeric(pheno$Flowering.time))

correlations <- data.frame(NA,NA,NA,NA)
names(correlations) <- c("variable", "R.sq", "coefficients", "coefficient.labels")

model = lm(trait.vector ~ Woody, data = pheno); correlations <- rbind(correlations, c("Woody", summary(model)$r.squared, paste(model$coefficients, collapse = ', '), paste(names(model$coefficients), collapse = ', ')))
model = lm(trait.vector ~ Perennial, data = pheno); correlations <- rbind(correlations, c("Perennial", summary(model)$r.squared, paste(model$coefficients, collapse = ', '), paste(names(model$coefficients), collapse = ', ')))
model = lm(trait.vector ~ Plant.height, data = pheno); correlations <- rbind(correlations, c("Plant.height", summary(model)$r.squared, paste(model$coefficients, collapse = ', '), paste(names(model$coefficients), collapse = ', ')))
model = lm(trait.vector ~ Petiole.length, data = pheno); correlations <- rbind(correlations, c("Petiole.length", summary(model)$r.squared, paste(model$coefficients, collapse = ', '), paste(names(model$coefficients), collapse = ', ')))
model = lm(trait.vector ~ Leaf.length, data = pheno); correlations <- rbind(correlations, c("Leaf.length", summary(model)$r.squared, paste(model$coefficients, collapse = ', '), paste(names(model$coefficients), collapse = ', ')))
model = lm(trait.vector ~ Leaf.width, data = pheno); correlations <- rbind(correlations, c("Leaf.width", summary(model)$r.squared, paste(model$coefficients, collapse = ', '), paste(names(model$coefficients), collapse = ', ')))
model = lm(trait.vector ~ Leaf.shape, data = pheno); correlations <- rbind(correlations, c("Leaf.shape", summary(model)$r.squared, paste(model$coefficients, collapse = ', '), paste(names(model$coefficients), collapse = ', ')))
model = lm(trait.vector ~ Inflorescence.type, data = pheno); correlations <- rbind(correlations, c("Inflorescence.type", summary(model)$r.squared, paste(model$coefficients, collapse = ', '), paste(names(model$coefficients), collapse = ', ')))
model = lm(trait.vector ~ Sepal.number, data = pheno); correlations <- rbind(correlations, c("Sepal.number", summary(model)$r.squared, paste(model$coefficients, collapse = ', '), paste(names(model$coefficients), collapse = ', ')))
model = lm(trait.vector ~ Sepal.length, data = pheno); correlations <- rbind(correlations, c("Sepal.length", summary(model)$r.squared, paste(model$coefficients, collapse = ', '), paste(names(model$coefficients), collapse = ', ')))
model = lm(trait.vector ~ Sepal.width, data = pheno); correlations <- rbind(correlations, c("Sepal.width", summary(model)$r.squared, paste(model$coefficients, collapse = ', '), paste(names(model$coefficients), collapse = ', ')))
model = lm(trait.vector ~ Sepal.color, data = pheno); correlations <- rbind(correlations, c("Sepal.color", summary(model)$r.squared, paste(model$coefficients, collapse = ', '), paste(names(model$coefficients), collapse = ', ')))
model = lm(trait.vector ~ Sepal.shape, data = pheno); correlations <- rbind(correlations, c("Sepal.shape", summary(model)$r.squared, paste(model$coefficients, collapse = ', '), paste(names(model$coefficients), collapse = ', ')))
model = lm(trait.vector ~ Petal.number, data = pheno); correlations <- rbind(correlations, c("Petal.number", summary(model)$r.squared, paste(model$coefficients, collapse = ', '), paste(names(model$coefficients), collapse = ', ')))
model = lm(trait.vector ~ Petal.length, data = pheno); correlations <- rbind(correlations, c("Petal.length", summary(model)$r.squared, paste(model$coefficients, collapse = ', '), paste(names(model$coefficients), collapse = ', ')))
model = lm(trait.vector ~ Petal.width, data = pheno); correlations <- rbind(correlations, c("Petal.width", summary(model)$r.squared, paste(model$coefficients, collapse = ', '), paste(names(model$coefficients), collapse = ', ')))
model = lm(trait.vector ~ Petal.color, data = pheno); correlations <- rbind(correlations, c("Petal.color", summary(model)$r.squared, paste(model$coefficients, collapse = ', '), paste(names(model$coefficients), collapse = ', ')))
model = lm(trait.vector ~ Petal.shape, data = pheno); correlations <- rbind(correlations, c("Petal.shape", summary(model)$r.squared, paste(model$coefficients, collapse = ', '), paste(names(model$coefficients), collapse = ', ')))
model = lm(trait.vector ~ Stamen.number, data = pheno); correlations <- rbind(correlations, c("Stamen.number", summary(model)$r.squared, paste(model$coefficients, collapse = ', '), paste(names(model$coefficients), collapse = ', ')))
model = lm(trait.vector ~ Stamen.length, data = pheno); correlations <- rbind(correlations, c("Stamen.length", summary(model)$r.squared, paste(model$coefficients, collapse = ', '), paste(names(model$coefficients), collapse = ', ')))
model = lm(trait.vector ~ Style.length, data = pheno); correlations <- rbind(correlations, c("Style.length", summary(model)$r.squared, paste(model$coefficients, collapse = ', '), paste(names(model$coefficients), collapse = ', ')))
model = lm(trait.vector ~ Seed.length, data = pheno); correlations <- rbind(correlations, c("Seed.length", summary(model)$r.squared, paste(model$coefficients, collapse = ', '), paste(names(model$coefficients), collapse = ', ')))
model = lm(trait.vector ~ Flowering.time, data = pheno); correlations <- rbind(correlations, c("Flowering.time", summary(model)$r.squared, paste(model$coefficients, collapse = ', '), paste(names(model$coefficients), collapse = ', ')))

correlations <- na.omit(correlations)
options(scipen = 999) # Do not display scientific notation

# Get relative "loadings" by normalizing correlations to sum to one
correlations$R.sq.corrected = as.numeric(correlations$R.sq)/(sum(as.numeric(correlations$R.sq)))

correlations
attach(correlations)
correlations.sorted <- correlations[order(R.sq.corrected),]
detach(correlations)

correlations.sorted

correlations.sorted.reduced = correlations.sorted
correlations.sorted.reduced$coefficients <- NULL
correlations.sorted.reduced$coefficient.labels <- NULL

write.table(correlations.sorted, "correlations.sorted.txt", quote = FALSE, row.names = FALSE, sep = "\t")
write.table(correlations.sorted.reduced, "correlations.sorted.reduced.txt", quote = FALSE, row.names = FALSE, sep = "\t")
