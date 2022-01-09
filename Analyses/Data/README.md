# Data for reproducible analyses


### Original ASV table from DADA2 pipeline: ASV_all.Rdata
It included 13 ASV tables with original ASVs (not truncated). Taxonomy assignments are not provided. Data sets Ba and Bb are not combined. 


### Multiple ASV tables and metadata: ASV_meta.Rdata
It included the common and top ASV tables, taxonomy assignments for ASVs in common and top table, orginal ASVs with truncated ASV sequence, and also metadata for each data set. Data sets Ba and Bb are not combined. 

The metadata data frame includes the following columns:\
Study\
Abbr (Abbreviation of study name)\
SampleID\
SubjectID\
GAAC (gestational age at collection) \
GAAD (gestational age at delivery) \
mAge (maternal age)
mBMI (maternal BMI)
mRace (maternal race, A: asian; B: black; C: caucasian; O: other)\
PPROM (1 is yes; 0 is no)\
pPTB (prior PTB; 1 is yes; 0 is no)\
sDel (spontaneous delivery, 1 is yes; 0 is no)\
GC (gestational complication)\
Intervention \
Preterm (if preterm birth based on 37 weeks; 1 is yes; 0 is no) \

### Genus-level data with different transformation: genus_data.Rdata
It includes genus-level data with the following data transformation methods: \
proportional data (prop) \
centered log-ratio (clr) \
centered log-ratio using top 4 abundance genus as reference (clr4) \
natural logarithm  (log) \
rank (rank) \
scaled proportional data (props) \
scaled entered log-ratio (clrs) \
scaled centered log-ratio using top 4 abundance genus as reference (clr4s) \
scaled natural logarithm  (logs) \
scaled rank (ranks) \

### Relative abundance with differnt taxaonmy level data: taxa_data.Rdata
It included relative abundance information of the following level:
ASV, genus, family, order, class, phylum
Both common and top table are included.


