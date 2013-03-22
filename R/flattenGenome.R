## Creates a GRanges object that mimics (with more annotations) the GRanges
## object returned from .buildRange,TranscriptDb. The big exception is that
## non-coding exons are annotated as such.
##
## Also if an anno.source.key is provided, we  try to annotated each range with
## gene symbol (and entrez.id -- although that's essentially already there).

if (interactive()) {
  library(GenomicFeatures)
  library(data.table)
}

## The dependency on data.table is to ensure that this function is as zippy
## as I can make it.
flattenGenome <- function(x, anno.source.key='UCSC Table') {
  stopifnot(inherits(x, 'TranscriptDb'))

  genome <- unname(genome(seqinfo(x))[1L])

  ## ---------------------------------------------------------------------------
  ## Protein coding transcripts
  utr5s <- fiveUTRsByTranscript(x)
  u5 <- unname(unlist(utr5s))
  values(u5)$transcript <- rep(names(utr5s), elementLengths(utr5s))
  values(u5)$exon.anno <- "utr5"

  cdss <- cdsBy(x, "tx")
  cs <- unname(unlist(cdss))
  values(cs)$transcript <- rep(names(cdss), elementLengths(cdss))
  values(cs)$exon.anno <- "cds"

  utr3s <- threeUTRsByTranscript(x)
  u3 <- unname(unlist(utr3s))
  values(u3)$transcript <- rep(names(utr3s), elementLengths(utr3s))
  values(u3)$exon.anno <- "utr3"

  cols <- c("seqnames", "start", "end", "strand", "transcript", "exon.anno")
  et <- rbind(as.data.frame(u5, stringsAsFactors=FALSE)[, cols],
              as.data.frame(cs, stringsAsFactors=FALSE)[, cols],
              as.data.frame(u3, stringsAsFactors=FALSE)[, cols])
  et <- defactor(et)

  ## ---------------------------------------------------------------------------
  ## non-coding transcripts
  all.tx <- exonsBy(x, 'tx')
  rest.tx <- all.tx[!names(all.tx) %in% et$transcript]
  rest.gr <- unname(unlist(rest.tx))
  values(rest.gr)$transcript <- rep(names(rest.tx), elementLengths(rest.tx))
  values(rest.gr)$exon.anno <- 'ncRNA'
  rest <- defactor(as.data.frame(rest.gr)[, cols])

  et <- rbind(et, rest)

  ## ---------------------------------------------------------------------------
  ## Add expected columns
  g2t <- transcriptsBy(x, "gene")
  gids <- rep(names(g2t), elementLengths(g2t))
  g2t <- defactor(transform(as.data.frame(unname(unlist(g2t))), gene_id=gids))
  xref <- match(et$transcript, g2t$tx_id)

  et$symbol <- g2t$tx_name[xref]
  et$gene <- g2t$gene_id[xref]

  ## Rig up our own exon_id. Gviz puts this in the "exon" column
  ## Using data.table to make the next two steps as fast as I know how
  et <- data.table(et, key=c('seqnames', 'strand', 'start', 'end'))
  uexons <- unique(et)
  uexons$exon <- 1:nrow(uexons)
  uexons <- uexons[, c(key(uexons), 'exon'), with=FALSE]
  et <- uexons[et]

  ## calculate exon_rank and put it in $exon.
  ## Rank of exons is increasing from 5' -> 3'
  ## Note that the xx[, yy := zz] syntax is specific to data.table
  ## and updates the data.table in place (really quick)
  setkeyv(et, c('transcript', 'strand', 'start'))
  et[, rank := if (strand[1] == '+') 1L:.N else .N:1L, by=transcript]

  features <- c(cds='protein_coding', utr3='utr', utr5='utr', ncRNA='ncRNA')
  et$feature <- features[et$exon.anno]
  et$density <- 1
  et$id <- et$exon

  ## Add gene column (entrez.id) and gene_symbol via addGeneInfo
  if (is.character(anno.source.key)) {
    anno.source.key <- tolower(anno.source.key)
    md <- transform(metadata(x), name=tolower(name))
    anno.source <- md[md$name == anno.source.key, 'value']
    if (length(anno.source) == 1L) {
      uber <- addGeneInfo(et, 'symbol', anno.source, genome)
      if (nrow(uber) != nrow(et)) {
        warning("Inconsistent number of rows in 'uber' table", immediate.=TRUE)
      }
      et <- uber
      if (all(c('symbol', 'gene.symbol') %in% names(et))) {
        et$tx_id <- et$symbol
        et$symbol <- ifelse(is.na(et$gene.symbol), et$tx_id, et$gene.symbol)
      }
    }
  }

  setkeyv(et, c('transcript', 'start'))

  ans <- with(et, GRanges(seqnames, IRanges(start, end), strand))
  rm.cols <- c('seqnames', 'strand', 'start', 'end', 'width')
  df.cols <- setdiff(names(et), rm.cols)

  ## ---------------------------------------------------------------------------
  ## do.call(DataFrame, as.list(et)) is super slow, so I'm currently stitching
  ## this together manually (and I'm surprised it's so much faster this way).
  meta <- DataFrame(transcript=et$transcript, exon=et$exon, rank=et$rank,
                    symbol=et$symbol, gene=et$gene, id=et$id,
                    feature=et$feature, density=et$density)
  for (mc in df.cols) {
    if (!mc %in% names(meta)) {
      meta[[mc]] <- et[[mc]]
    }
  }
  ## ---------------------------------------------------------------------------
  values(ans) <- meta

  ## Add TranscriptDb seqinfo to the result for posterity
  if (!all(seqlevels(ans) %in% seqlevels(x))) {
    warning("seqlevels of result look suspicious", immediate.=TRUE)
  } else {
    si <- seqinfo(x)[seqlevels(x)[seqlevels(x) %in% seqlevels(ans)]]
    new2old <- match(seqlevels(si), seqlevels(ans))
    seqinfo(ans, new2old) <- si
  }

  ans
}

## as.data.frame(GRanges) converts strings to factors :-/
defactor <- function(x) {
  for (col in names(x)) {
    if (is.factor(x[[col]])) x[[col]] <- as.character(x[[col]])
  }
  x
}

##' Wrapper function to add more meta columns to the range cache
##' query.col is set to 'symbol' because Gviz is using the symbol column
##' to store refseq IDs -- this assumption probably needs to be fixed.
##'
##' @param x Assumed to be a data.table
##' @param key.col The column that holds the transcript id to use as query
##' against the org.*.db
##' @param genome The genome accession number
##' @param anno.source refGene? knownGene? This tells us
##' what to query in the org.*.db package.
addGeneInfo <- function(x, key.col='symbol', anno.source='refGene',
                        genome='hg19') {
  if (!key.col %in% names(x)) {
    warning("query.col '", key.col, "' not in `x`", immediate.=TRUE)
    return(x)
  }
  db.key.col <- annoSource2DbKeyCol(anno.source)
  if (is.null(db.key.col)) {
    warning("Unknown column to query org.*.db by ", db.key.col, immediate.=TRUE)
    return(x)
  }

  org.db <- genome2orgDb(genome)
  if (is.null(org.db)) {
    warning("Unsupported genome <--> org.db mapping for genome: ", genome,
            immediate.=TRUE)
    return(x)
  }
  if (!require(org.db, character.only=TRUE)) {
    warning(org.db, " library not installed, skipping ...", immediate.=TRUE)
    return(x)
  }

  keys <- unique(as.character(x[[key.col]]))
  keys <- keys[!is.na(keys)]

  db <- get(org.db)
  db.cols <- orgGeneInfoColumns(org.db)
  if (is.null(db.cols)) {
    warning("Undefined return columns for ", org.db, immediate.=TRUE)
    return(x)
  }
  if (!all(db.cols %in% cols(db))) {
    warning("Unexpected query column names for ", org.db, immediate.=TRUE)
    return(x)
  }

  ## TODO: Poke the BioC folks w/ an easy test case to show big differences
  ##       in speed between the very nice `select` interface vs. hand weaving
  ##       the results I need.
  ## Sadly the `select` stuff is super slow
  ## info <- select(db, keys, db.cols, db.key.col)
  ## info <- info[!duplicated(info[[db.key.col]]),,drop=FALSE]
  ## names(info) <- c(key.col, names(db.cols))

  info <- fetchDbInfo(org.db, db.key.col, keys)
  if (is.null(info)) {
    warning("Fetching addition gene info unexpectedly failed", immediate.=TRUE)
    return(x)
  }
  names(info)[1] <- key.col
  m <- merge(x, info, by=key.col, all.x=TRUE)
  new.cols <- setdiff(names(m), names(x))
  m[, c(names(x), new.cols), with=FALSE]
}

## get entrez ids and symbols
## database is already loaded
fetchDbInfo <- function(org.db, db.key.col, keys) {
  bn <- gsub("\\.db", "", org.db)

  ## entrez
  if (db.key.col == "UCSCKG") {
    map.name <- paste0(bn, db.key.col)
    entrez.map <- revmap(get(map.name))
  } else {
    if (!db.key.col %in% c("REFSEQ", "ENSEMBL")) {
      return(NULL)
    }
    entrez.map <- get(paste0(bn, db.key.col, "2EG"))
  }
  entrez.ids <- mget(keys, entrez.map, ifnotfound=NA)
  entrez.ids <- sapply(entrez.ids, '[[', 1L)
  entrez.ids <- entrez.ids[!is.na(entrez.ids)]

  ans <- data.frame(key=names(entrez.ids),
                    entrez.id=unname(entrez.ids),
                    stringsAsFactors=FALSE)
  ans <- subset(ans, !duplicated(key))

  ## symbols (yeast should be SGD, but ugh -- forget this for now)
  symbol.map <- tryCatch(get(paste0(bn, "SYMBOL")), error=function(e) NULL)
  if (!is.null(symbol.map)) {
    symbols <- mget(entrez.ids, symbol.map, ifnotfound=NA)
    ans$gene.symbol <- sapply(symbols, '[[', 1L)
  }

  ans
}

genome2orgDb <- function(genome) {
  org.dbs <- list(hg19="org.Hs.eg.db",
                  mm9="org.Mm.eg.db", mm10="org.Mm.eg.db",
                  dm3="org.Dm.eg.db")
  org.dbs[[genome]]
}

annoSource2DbKeyCol <- function(anno.source) {
  a2k <- list(knownGene="UCSCKG", refGene="REFSEQ")
  a2k[[anno.source]]
}

orgGeneInfoColumns <- function(org.db) {
  cols <- list(org.Hs.eg.db=c(entrez.id="ENTREZID", gene.symbol="SYMBOL"),
               org.Mm.eg.db=c(entrez.id="ENTREZID", gene.symbol="SYMBOL"),
               org.Dm.eg.db=c(entrez.id="ENTREZID", gene.symbol="SYMBOL"),
               org.Sc.sgd.db=c(entrez.id="ENTREZID", gene.symbol="SGD"))
  cols[[org.db]]
}

## `keys` is the character() of ids to use as the keys for the query
## `db.key.col` is the name of the column that `keys` are in, eg.
## 'REFSEQ', 'UCSCKG', etc.
getGeneInfo.org.Hs.eg.db <- function(keys, db.key.col="REFSEQ") {
  org.db <- get("org.Hs.eg.db")
  info <- select(org.db, , c("ENTREZID", "SYMBOL"), query.id)
  names(info) <- c(query.col, 'entrez.id', 'gene.symbol')
  info[!duplicated(info[[query.col]]),]
}

getGeneInfo.org.Dm.eg.db <- function(x) {
}

getGeneInfo.org.Mm.eg.db <- function(x) {
  org.db <- get("org.Mm.eg.db")
  info <- select(org.db, , c("ENTREZID", "SYMBOL"), query.id)
  names(info) <- c(query.col, 'entrez.id', 'gene.symbol')
  info[!duplicated(info[[query.col]]),]
}
