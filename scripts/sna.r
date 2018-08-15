#!/usr/bin/env Rscript

#
# Runs the social network analysis using Peter Hoff's GBME
#
suppressMessages(library("optparse"))
suppressMessages(library("xtable"))
suppressMessages(library("R.utils"))

cargs = commandArgs(trailingOnly = FALSE)
source_dir = dirname(sub("^--file=", "", cargs[grep("^--file=", cargs)]))
source(file.path(source_dir, "gbme.r"))

option_list = list(
  make_option(
    c("-f", "--file"),
    default = "",
    type = "character",
    help = "Similarity matrix",
    metavar = "character"
  ),
  make_option(
    c("-o", "--out_dir"),
    default = "",
    type = "character",
    help = "outdir",
    metavar = "character"
  ),
  make_option(
    c("-s", "--sna"),
    default = "sna-gbme.pdf",
    type = "character",
    help = "GBME SNA output filename",
    metavar = "character"
  ),
  make_option(
    c("-n", "--number"),
    default = as.integer(20000),
    type = "integer",
    help = "Number of GBME iterations",
    metavar = "integer"
  ),
  make_option(
    c("-a", "--alias"),
    default = "",
    type = "character",
    help = "Sample name alias file",
    metavar = "alias"
  )
);

opt_parser   = OptionParser(option_list = option_list);
opt          = parse_args(opt_parser);
matrix.file  = opt$file  
out.dir      = opt$out_dir
n_iter       = opt$number
alias.file   = opt$alias
sna.filename = opt$sna

Y = as.matrix(read.table(matrix.file, header = TRUE))

if (nchar(matrix.file) == 0) {
  stop("Missing matrix file argument")
}

if (!file.exists(matrix.file)) {
  stop(paste("Invalid --matrix", matrix.file))
}

if (nchar(out.dir) == 0) {
  out.dir = dirname(matrix.file)
}

if (!dir.exists(out.dir)) {
  printf("Creating outdir '%s'\n", out.dir)
  dir.create(out.dir)
}
out.dir = normalizePath(out.dir)
setwd(out.dir)

GBME_OUT = file.path(out.dir, "gbme.out")
if (file.exists(GBME_OUT)) {
  file.remove("Failed to run GBME!")
}

# Look for the "*.meta" files 
meta_dir = file.path(out.dir, "meta")

printf("Reading %s matrix\n", matrix.file)
Y = as.matrix(read.table(matrix.file, header = TRUE))
n = nrow(Y)
Xss = NULL
k = 0

if (dir.exists(meta_dir)) {
  meta_files = list.files(path = meta_dir, pattern = "*.meta$")
  print(meta_files)
  k = length(meta_files)
  
  printf("Found %s file%s in meta dir '%s'\n", k, 
         if (k == 1) { '' } else {'s'}, meta_dir)
  
  if (k > 0) {
    Xss = array(NA, dim = c(n,n,k))
    
    for (i in 1:k) {
      file = file.path(meta_dir, meta_files[i])
      d = as.matrix(read.table(file, header = TRUE))
      if (nrow(d) == n) {
        Xss[,,i] = d
      }
    }
  }
}

printf("Running GBME with %s scans\n", n_iter)
Z_path = file.path(out.dir, "Z")

if (is.null(Xss)) {
  gbme(Y = Y, fam = "gaussian", k = 2, direct = F, NS = n_iter, odens = 10, ofilename = GBME_OUT, zfilename = Z_path)
} else {
  tryCatch(gbme(Y = Y, Xss, fam = "gaussian", k = 2, direct = F, NS = n_iter, odens = 10, ofilename = GBME_OUT, zfilename = Z_path),
           error = function(c) {
               stop("Error running GBME with metadata -- TRY WITHOUT!")
           }
           )
}

printf("looking for %s", GBME_OUT)
if (!file.exists(GBME_OUT)) {
  stop("Failed to run GBME!")
}

printf("Using GBME_OUT '%s'\n", GBME_OUT)
OUT = read.table(GBME_OUT, header = T)

# posterior samples, dropping
# the first half of the chain
# to allow for burn in
PS <- OUT[OUT$scan > round(max(OUT$scan) / 2),-(1:3)]

#
# analysis of latent positions
#
Z <- read.table(Z_path)

#
# convert to an array
#
nss <- dim(OUT)[1]
n <- dim(Z)[1] / nss
k <- dim(Z)[2]
PZ <- array(dim = c(n,k,nss))
for (i in 1:nss) {
  PZ[,,i] <- as.matrix(Z[((i - 1) * n + 1):(i * n) ,])
}

PZ <- PZ[,,-(1:round(nss / 2))] #drop first half for burn in

#
# find posterior mean of Z %*% t(Z)
#
ZTZ <- matrix(0,n,n)
for (i in 1:dim(PZ)[3]) {
  ZTZ <- ZTZ + PZ[,,i] %*% t(PZ[,,i])
}
ZTZ <- ZTZ / dim(PZ)[3]

#
# a configuration that approximates posterior mean of ZTZ
#
tmp <- eigen(ZTZ)
Z.pm <- tmp$vec[, 1:k] %*% sqrt(diag(tmp$val[1:k]))

#now transform each sample Z to a common orientation
for (i in 1:dim(PZ)[3]) {
  PZ[,,i] <- proc.rr(PZ[,,i],Z.pm)
}

#
# a two dimensional plot of "mean" latent locations
# and marginal confidence regions
#
r <- atan2(Z.pm[,2], Z.pm[,1])
r <- r + abs(min(r))
r <- r / max(r)
g <- 1 - r
b <- (Z.pm[,2] ^ 2 + Z.pm[,1] ^ 2)
b <- b / max(b)

par(mfrow = c(1,1))
sna_path = file.path(out.dir, sna.filename)
printf('Creating "%s"\n', sna_path)
pdf(sna_path, width = 6, height = 6)
plot(
  Z.pm[,1],
  Z.pm[,2],
  xlab = "",
  ylab = "",
  type = "n",
  xlim = range(PZ[,1,]),
  ylim = range(PZ[,2,])
)

abline(h = 0,lty = 2);abline(v = 0,lty = 2)

for (i in 1:n) {
  points(PZ[i,1,], 
         PZ[i,2,], 
         pch = 46, 
         col = rgb(r[i], g[i], b[i]))
}

# add labels here
labels = rownames(Y)
if (length(alias.file) > 0 & file.exists(alias.file)) {
  printf("Using alias file '%s'\n", alias.file)
  aliases = read.table(alias.file, header=T, as.is=T)
  for (i in 1:length(labels)) {
    label = labels[i]
    alias = aliases[aliases$name == label, "alias"]
    if (length(alias) > 0) {
      printf("Alias '%s' -> '%s'\n", label, alias)
      labels[i] = alias
    }
  }
}

text(Z.pm[,1],Z.pm[,2], cex = 0.3, labels = labels)   
invisible(dev.off())

printf("Done, see output dir '%s'\n", out.dir)
