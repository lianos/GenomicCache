setClassUnion("MaybeTranscriptDb", c('TranscriptDb', 'NULL'))

setClass("GenomicFeaturesX", contains="VIRTUAL")

## Load up the transcripts, exonsBy, threeUTRsByTranscript,
## fiveUTRsByTranscript. Seems we want to query the *Ranges objects
## in memory
setClass("GenomicCache",
         representation(.txdb='MaybeTranscriptDb',
                        ## .transcripts="GRanges",
                        ## .exons="GRangesList", ## byTranscript
                        ## .utr5="GRangesList",  ## byTranscript
                        ## .utr3="GRangesList",
                        .cache="environment"))

## OO Way to deal with Gene objects
setClass("GFGene",
         representation(.id='character',
                        .entrez.id='character',
                        .symbol='character',
                        .chromosome='factor',
                        .strand='factor',
                        .exons="GRangesList",
                        .cds="GRangesList",
                        .utr5="GRangesList",
                        .utr3="GRangesList",
                        .transcript.names="character",
                        .genome='character',
                        .cache='environment'))

setClass("RefSeqGene", contains="GFGene")
setClass("EnsemblGene", contains="GFGene")
setClass("AceviewGene", contains="GFGene")
