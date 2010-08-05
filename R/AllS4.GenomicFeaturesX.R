setClass("GenomicFeaturesX", contains="VIRTUAL")


## Load up the transcripts, exonsBy, threeUTRsByTranscript,
## fiveUTRsByTranscript. Seems we want to query the *Ranges objects
## in memory
setClass("GenomiCache",
         representation(transcripts="GRanges",
                        exons="GRangesList", ## byTranscript
                        utr5="GRangesList",  ## byTranscript
                        utr3="GRangesList"))


## OO Way to deal with Gene objects
setClass("GFGene", representation(.id='character',
                                  .symbol='character',
                                  .chromosome='factor',
                                  .strand='factor',
                                  .transcripts="GRangesList",
                                  .transcript.names="character",
                                  .txBounds='IRanges',
                                  .cdsBounds'IRanges',
                                  .cache='environment'))
setClass("RefSeqGene", contains="GFGene")
setClass("EnsemblGene", contains="GFGene")
setClass("AceviewGene", contains="GFGene")
