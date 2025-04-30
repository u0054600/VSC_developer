#suppressMessages(library(optparse))
#suppressMessages(library(ShortRead))
#suppressMessages(library(tidyverse))
#suppressMessages(library(furrr))
#suppressMessages(library(data.table))

# in v18 we added '-num_threads 15', in the blast command  and fixed a fa2 new line that was commented out.
# in v19 we fixed the patch on line 1100 
suppressMessages(library("optparse",  lib="/staging/leuven/stg_00080/vrc/Projects/Project_Rushuang/Rpack"))
suppressMessages(library("ShortRead",  lib="/staging/leuven/stg_00080/vrc/Projects/Project_Rushuang/Rpack"))
suppressMessages(library("ggplot2",  lib="/staging/leuven/stg_00080/vrc/Projects/Project_Rushuang/Rpack"))
suppressMessages(library("forcats",  lib="/staging/leuven/stg_00080/vrc/Projects/Project_Rushuang/Rpack"))
suppressMessages(library("tzdb",  lib="/staging/leuven/stg_00080/vrc/Projects/Project_Rushuang/Rpack"))
suppressMessages(library("forcats",  lib="/staging/leuven/stg_00080/vrc/Projects/Project_Rushuang/Rpack"))
suppressMessages(library("readr",  lib="/staging/leuven/stg_00080/vrc/Projects/Project_Rushuang/Rpack"))
suppressMessages(library("tidyr",  lib="/staging/leuven/stg_00080/vrc/Projects/Project_Rushuang/Rpack"))
suppressMessages(library("lubridate",  lib="/staging/leuven/stg_00080/vrc/Projects/Project_Rushuang/Rpack"))
suppressMessages(library("tidyverse",  lib="/staging/leuven/stg_00080/vrc/Projects/Project_Rushuang/Rpack"))
suppressMessages(library("data.table",  lib="/staging/leuven/stg_00080/vrc/Projects/Project_Rushuang/Rpack"))
suppressMessages(library("furrr",  lib="/staging/leuven/stg_00080/vrc/Projects/Project_Rushuang/Rpack"))
suppressMessages(library("vegan",  lib="/staging/leuven/stg_00080/vrc/Projects/Project_Rushuang/Rpack"))

 set.seed(123)


print.noquote("libraries loaded") 

option_list = list(
  make_option(c("--sample"), action="store", help = "sample name"),
  make_option(c("--data_path"), action="store", help = "path to fastq files"),
  make_option(c("--kraken_path"), action="store", help = "path to Kraken2 executable"),
  make_option(c("--kraken_db_path"), action="store", help = "path to Kraken2 reference database"),
  make_option(c("--seqkit_path"), action="store", help = "path to SeqKit"),
  make_option(c("--star_path"), action="store", help = "path to STAR"),
  make_option(c("--star_genome_dir"), action="store", help = "path to STAR genome index (for --genomeDir parameter)"),
  make_option(c("--model_org_taxids"), action="store", default = 'NA', help = "path to .txt file with model organism taxids"),
  make_option(c("--blast_path"), action="store", help = "path to blastn"),
  make_option(c("--blast_db_path"), action="store", help = "path to blast reference database"),
  make_option(c("--prism_path"), action="store", help = "path to PRISM directory"),
  make_option(c("--min_read_per"), action="store", default = 10^4, help = "minimum reads per X to analyze. Default is 1 per 10,000 (min_read_per=10^4)"),
  make_option(c("--min_uniq_frac"), action="store", default = 5, help = "minimum ratio of unique k-mer to Kraken reads. Default=5"),
  make_option(c("--paired"), action="store", default = T, help = "paired end (T) or single end (F) reads"),
  make_option(c("--max_sample"), action="store", default = 1000, help = "max reads to sample"),
  make_option(c("--barcode_only"), action="store", default = F, help = "does read 1 only contain a barcode (T) or not (F)"),
  make_option(c("--fq1_end"), action="store", default = "_1.fastq", help = "fastq file name end format after sample name (e.g. sample_1.fastq, or sample_R1.fastq)"),
  make_option(c("--fq2_end"), action="store", default = "_2.fastq", help = "fastq file name end format after sample name (e.g. sample_1.fastq, or sample_R1.fastq)")
)
opt = parse_args(OptionParser(option_list = option_list))

sample = opt$sample
data_path = opt$data_path
kraken_path = opt$kraken_path
kraken_db_path = opt$kraken_db_path
seqkit_path = opt$seqkit_path
star_path = opt$star_path
star_genome_dir = opt$star_genome_dir
blast_path = opt$blast_path
blast_db_path = opt$blast_db_path
prism_path = ifelse(str_detect(opt$prism_path, '/$'), opt$prism_path, paste0(opt$prism_path, '/'))
model_org_taxids = ifelse(opt$model_org_taxids == 'NA', paste0(prism_path, opt$model_org_taxids), opt$model_org_taxids)
ranks = 'S'
min_read_per = opt$min_read_per
min_uniq_frac = opt$min_uniq_frac
paired = opt$paired %>% as.logical()
barcode_only = opt$barcode_only %>% as.logical()
max_sample = opt$max_sample
out_path = paste0(ifelse(str_detect(data_path, '/$'), data_path, paste0(data_path, '/')), sample, '_prism/data/')
out_path_final = paste0(ifelse(str_detect(data_path, '/$'), data_path, paste0(data_path, '/')), sample, '_prism/')
fq1_end = opt$fq1_end
fq2_end = opt$fq2_end



cat(paste('\nSample name:', sample, '\n'))

load (paste0(out_path, sample, ".image.line 1088.R"))



############## multimapping and gene/protein mapping 

cat(paste('Started final multimapping and gene/protein mapping', Sys.time(), '\n'))

f = list.files(paste0(out_path), full.names = T) %>% str_subset('-final-blast')
blast2 = list()
for(i in 1:length(f)){
  x = read.csv(f[i], header = F, stringsAsFactors = F)
  colnames(x) = c('id', 'sacc', 'staxids', 'pos', 'ppos', 'bitscore', 'strand', 'qcovs')
  x$staxids = as.character(x$staxids)
  blast2[[i]] = tibble(read = i, x)
}
blast2 = bind_rows(blast2)

blast2 = blast2 %>% 
  tibble() %>% 
  group_by(id) %>% 
  mutate(human = ifelse(any(staxids == 9606), 1, 0)) %>% 
  subset(human == 0) %>% 
  dplyr::select(-human) %>% 
  mutate(staxids = str_remove(staxids, ';.*'))

rank_order = data.frame(rank = c('r','d','k','p','c','o','f','g','s'), order = 1:9)

x = table(blast2$id, blast2$staxids) %>% as.data.frame.matrix() 
x = x %>% mutate(duplicate_group = group_indices(x, across(everything())))
z = distinct(x)


save.image(paste0(out_path, sample, ".image.line 1100.R"))
##future::plan(multisession)
##df = future_map(1:nrow(z), .progress = T, function(i){
##  mpa[which(kr2$V7 %in% colnames(z)[which(z[i,]> 0)]), 'V1'] %>% as.character() %>% 
##    strsplit('\\|') %>% 
##    unlist() %>% 
##    table() %>% 
##    as.data.frame() %>%
##    arrange(-Freq) %>% 
##    filter(Freq == max(Freq)) %>% 
##    mutate(rank = str_remove(., '_.*'), name = str_remove(., '.*__') %>% str_replace_all('[^[:alnum:]]', ' ') %>% str_squish()) %>% 
##    dplyr::select(-.) %>% 
##    dplyr::select(-rank) %>% 
##    left_join(data.frame(name = kr2$V8 %>% str_replace_all('[^[:alnum:]]', ' ') %>% str_squish() , taxid = kr2$V7, rank = kr2$V6), by = c('name')) %>% 
##    mutate(rank = str_to_lower(rank) %>% str_remove_all('[0-9]+')) %>% 
##    left_join(rank_order, by = 'rank') %>% 
##    subset(order == max(order, na.rm = T)) %>% 
##    tail(n=1) %>% 
##    dplyr::select(rank, taxid, name) %>% 
##    mutate(duplicate_group = z$duplicate_group[i], name = str_replace_all(name, '_', ' '))
##}) %>% bind_rows() %>% tibble()

## The block above is patched? 



df3LIST<- list()

for (i in which(kr2$V7 %in% colnames(z))){
df3LIST [[i]] <- mpa[i, 'V1'] %>% as.character() %>%
     strsplit('\\|') %>%
     unlist() %>%
     table() %>%
     as.data.frame() %>%
     arrange(-Freq) %>%
     filter(Freq == max(Freq)) %>%
     mutate(rank = str_remove(., '_.*'), name = str_remove(., '.*__') %>% str_replace_all('[^[:alnum:]]', ' ') %>% str_squish()) %>%
     dplyr::select(-.) %>%
     dplyr::select(-rank) %>%
     left_join(data.frame(name = kr2$V8 %>% str_replace_all('[^[:alnum:]]', ' ') %>% str_squish() , taxid = kr2$V7, rank = kr2$V6), by = c('name')) %>%
     mutate(rank = str_to_lower(rank) %>% str_remove_all('[0-9]+')) %>%
     left_join(rank_order, by = 'rank') %>%
     subset(order == max(order, na.rm = T)) %>%
     tail(n=1) %>%
     dplyr::select(rank, taxid, name) %>%
     mutate(duplicate_group = z$duplicate_group[i], name = str_replace_all(name, '_', ' '))
 }


df <- do.call(rbind, df3LIST) %>% tibble



save.image(paste0(out_path, sample, ".image.line 1171.R"))

print.noquote("end of patch at line 1171")

multidf = x %>% rownames_to_column('id') %>% left_join(df, by = 'duplicate_group') %>% dplyr::select(id,rank,taxid,name)

write.table(multidf, file = paste0(out_path, sample, '-multidf.txt'))

# gene/protein mapping
future::plan(multisession)
uacc = unique(blast2$sacc) %>% as.character()

print.noquote("small future map below line 1160")
ii = future_map(1:length(uacc), .progress = T, function(x){str_which(genbank_dict$accession, uacc[x])}) %>% unlist() %>% unique()
print.noquote("after.small future map line 1160")

genbank = list()
for(i in 1:length(ii)){
  # print(i)
  genbank[[i]] = readRDS(paste0(prism_path, 'genbank/rds/', genbank_dict$file[ii[i]]))
}
genbank = data.table::rbindlist(genbank) %>% tibble()
genbank = mutate_all(genbank, .funs = as.character)
genbank$start = as.numeric(genbank$start)
genbank$end = as.numeric(genbank$end)


x = blast2 %>% data.table::data.table()
x = x[, .(Group_List = .(list(.SD))), by = sacc]
s = x[[1]]; b = x[[2]]; names(b) = x$sacc
dt <- data.table::as.data.table(genbank)
dt_list <- dt[, .(Group_List = .(list(.SD))), by = accession]

bl = list()
for(i in 1:length(b)){
  # print(i)
  new_b = b[[i]] %>% data.frame()
  if(names(b)[i] %in% dt_list$accession == F){bl[[i]] = new_b; next}
  g = dt_list[[2]][which(dt_list$accession == names(b)[i])][[1]] %>% 
    data.frame() %>% distinct() %>% subset(end > start) %>% rownames_to_column('rn') %>% arrange(start)
  if(nrow(g) == 0){bl[[i]] = new_b; next}
  xx = g %>% 
    dplyr::select(-taxon,-locus) %>% 
    pivot_longer(-c(rn,definition,version,organism,gene,product,protein)) %>% 
    subset(!is.na(value))
  re = c(); for(j in 1:(nrow(xx)-1)){re[j] = ifelse(xx$value[j] > xx$value[j+1], 1, 0)}
  if(any(re == 1)){xx = xx[-which(re == 1), ]}
  int = findInterval(new_b$pos, xx$value, rightmost.closed = T)
  if(any(int == 0)){int[int == 0] = NA}
  if(any(names(table(xx$name[int])) == 'end')){int[xx$name[int] == 'end'] = NA}
  bl[[i]] = new_b %>% mutate(sacc = names(b)[i]) %>% cbind(g[xx$rn[int], ] %>% dplyr::select(-c(rn,taxon))) %>% tibble()
}


def = read.delim(paste0(prism_path, 'cog-20.def.tab'), header = F)
fun = read.delim(paste0(prism_path, 'fun-20.tab'), header = F)

bfinal = bind_rows(bl) %>% arrange(id, product, gene) %>% distinct(id, .keep_all = T) %>% 
  mutate(cog = def$V3[match(str_to_lower(gene), str_to_lower(def$V4))], 
         cat = def$V2[match(str_to_lower(gene), str_to_lower(def$V4))]) %>% 
  mutate(cat = fun$V3[match(cat, fun$V1)]) %>% 
  left_join(multidf %>% dplyr::select(id, rank, name) %>% dplyr::rename(tax_name = name), by = 'id') %>% 
  dplyr::select(id, staxids, rank, tax_name, everything()) %>% 
  dplyr::select(-organism, -locus)

print.noquote("get taxonomy for each read")
## get taxonomy for each read
rank_order2 = data.frame(rank = c('R','R1','R2','R3','R4','R5','R6','R7','R8','R9','R10',
                                  'D','D1','D2','D3','D4','D5','D6','D7','D8','D9','D10',
                                  'K','K1','K2','K3','K4','K5','K6','K7','K8','K9','K10',
                                  'P','P1','P2','P3','P4','P5','P6','P7','P8','P9','P10',
                                  'C','C1','C2','C3','C4','C5','C6','C7','C8','C9','C10',
                                  'O','O1','O2','O3','O4','O5','O6','O7','O8','O9','O10',
                                  'F','F1','F2','F3','F4','F5','F6','F7','F8','F9','F10',
                                  'G','G1','G2','G3','G4','G5','G6','G7','G8','G9','G10',
                                  'S','S1','S2','S3','S4','S5','S6','S7','S8','S9','S10'
),
order = 1:99)

tax = unique(bfinal$staxids)

tax_df1 = list()
for(i in 1:length(tax)){
  tax_df1[[i]] = data.frame(name = mpa[which(kr2$V7 == tax[i]),'V1'] %>% str_split('\\|') %>% unlist() %>% str_remove('.*__') %>% str_replace_all('_', ' '), 
                            taxid = mpa[which(kr2$V7 == tax[i]), 'taxid'] %>% str_remove('\\*') %>% strsplit('\\*') %>% unlist() %>% as.character()) %>% 
    left_join(kr2 %>% dplyr::select(V6,V7) %>% dplyr::rename(rank = V6, taxid = V7) %>% mutate(taxid = as.character(taxid)), by = 'taxid') %>% 
    pivot_longer(-rank) %>%
    mutate(name = paste0(rank,'_',name), staxids = tax[i]) 
}

tax_df1 = tax_df1 %>% bind_rows() %>% mutate(name = factor(name, levels = c(paste0(rank_order2$rank, '_taxid'), paste0(rank_order2$rank, '_name')))) %>% 
  distinct(name, value, staxids, .keep_all = T) %>% arrange(name) %>% 
  pivot_wider(id_cols = staxids, names_from = name, values_from = value)

tax_df2 = list()
for(i in 1:length(tax)){
  tax_df2[[i]] = data.frame(name = mpa[which(kr2$V7 == tax[i]),'V1'] %>% str_split('\\|') %>% unlist() %>% str_remove('.*__') %>% str_replace_all('_', ' '), 
                            taxid = mpa[which(kr2$V7 == tax[i]), 'taxid'] %>% str_remove('\\*') %>% strsplit('\\*') %>% unlist() %>% as.character()) %>% 
    left_join(kr2 %>% dplyr::select(V6,V7) %>% dplyr::rename(rank = V6, taxid = V7) %>% mutate(taxid = as.character(taxid)), by = 'taxid') %>% 
    mutate(staxids = tax[i])
}
tax_df2 = bind_rows(tax_df2)

saveRDS(list(tax_df1 = tax_df1, tax_df2 = tax_df2), file = paste0(out_path, sample, '-taxdf.RDS'))

save.image(paste0(out_path, sample, ".image.line 1275.R"))
### counts table and prism predictions 
prod = readRDS(paste0(out_path, sample, '-products.RDS')) %>%   group_by(taxid, name) %>% 
  summarize(n = n(), fprod = n(), fugene = length(unique(gene)), fuprod = length(unique(product)),
            prod_div = sapply(1:10, function(x) sample(product, 100, replace = T) %>% table() %>% vegan::diversity()) %>% mean(),
            gene_div = sapply(1:10, function(x) sample(gene, 100, replace = T) %>% table() %>% vegan::diversity()) %>% mean(),
            .groups = 'drop') %>% 
  mutate(fprod = fprod/sum(fprod), fugene = fugene/sum(fugene), fuprod = fuprod/sum(fuprod))

counts = read.delim(paste0(out_path, sample, '.kraken.report.txt'), header =F) %>% tibble() %>% mutate(V8 = str_squish(V8)) %>% summarize(species = V8, n2=V2, n3=V3, rank = V6, taxid = V7, uniq = V5, fmicro = V2/sum(V2[V7 %in% c(2,4751,10239)]), fcontam = V2/sum(V2[V7 %in% c(2100, 1747, 1282, 75775, 562, 1270, 1423,573,34062,470)])) %>%subset(str_detect(rank, 'S')) %>% dplyr::select(-rank)
kphylo = read.table(paste0(out_path, sample, '-kmerphylo.txt')) %>% pivot_longer(-c(taxid)) %>% group_by(taxid) %>% mutate(value = value/sum(value)) %>% pivot_wider(id_cols = c(taxid), names_from = name, values_from = value) %>% ungroup()              
mc = read.table(paste0(out_path, sample, '-misclass.txt')) %>% pivot_longer(-c(taxid, rank)) %>% group_by(taxid, name) %>% mutate(value = value/value[rank == 'k']) %>%subset(rank != 'k') %>% pivot_wider(id_cols = c(taxid), names_from = c(rank, name), values_from = value) %>% ungroup()
mb = read.table(paste0(out_path, sample, '-multiblast.txt'))  

test = mb %>% 
  left_join(counts, by = c('taxid')) %>% 
  left_join(kphylo, by = c('taxid')) %>% 
  left_join(prod, by = c('taxid')) %>% 
  left_join(mc, by = c('taxid')) %>% 
  dplyr::select(species, taxid, everything())

write.table(test, file = paste0(out_path, sample, '-xgmat.txt'))

xg = readRDS(file = paste0(prism_path, 'prismxg.RDS'))

pred = predict(xg, test %>% dplyr::select(xg$feature_names) %>% as.matrix())

counts = read.delim(paste0(out_path, sample, '.kraken.report.txt'), header =F)

ct = bfinal %>% group_by(staxids) %>% summarize(n = n()) %>% 
  left_join(tax_df2, by = 'staxids') %>% group_by(name, taxid, rank) %>% summarize(n = sum(n), .groups = 'drop') %>% 
  mutate(rank = factor(rank, levels = rank_order2$rank)) %>% arrange(rank)

ct = rbind(data.frame(name = 'Total', taxid = '00', rank = 'T', n = sum(counts$V2[counts$V7 %in% c(0,1)])),
           data.frame(name = 'Homo sapiens', taxid = '9606', rank = 'S', n = sum(counts$V2[counts$V7 %in% c(9606)])),
           ct
)

ct = left_join(ct, data.frame(taxid = test$taxid %>% as.character(), pred = pred), by = 'taxid')

write.csv(ct, file = paste0(out_path_final, sample, '-counts.csv'), row.names = F, quote = F)

# add contamination score to results file
bfinal = left_join(bfinal %>% mutate(staxids = as.character(staxids)), data.frame(staxids = test$taxid %>% as.character(), pred = pred), by = 'staxids')
data.table::fwrite(bfinal, file = paste0(out_path_final, sample, '-results.csv'))

# make PRISM microbiome fasta file
headers <- headers %>% str_remove('kraken.*') %>% str_trim()
fa1_new = subset(fa1_new, header_ids %in% bfinal$id) # subset for reads that were blasted
ShortRead::writeFasta(fa1_new, file = paste0(out_path, sample, '-new_1.fa'))
headers = subset(headers, header_ids %in% bfinal$id)
header_ids = subset(header_ids, header_ids %in% bfinal$id)
df = data.frame(header = headers, id = header_ids) %>% left_join(bfinal %>% select(id, staxids, sacc, pos), by = 'id')
new_headers = paste0(">", df$header, ' | PRISM | staxids:', df$staxids, ' sacc:', df$sacc, ' pos:', df$pos)
write(new_headers, file = paste0(out_path_final, 'new_headers.txt'))
str = paste0("awk 'NR==FNR { h[++i] = $0; next } /^>/ { print h[++j]; next } { print }' ", out_path_final, 'new_headers.txt ', 
             out_path, sample, '-new_1.fa > ', out_path_final, sample, '_1.fa')
system(str)
str = paste0('rm ', out_path, sample, '_1.fa ', out_path, sample, '-new_1.fa ', out_path_final, 'new_headers.txt')
system(str)

if(paired == T){
  headers2 <- headers2 %>% str_remove('kraken.*') %>% str_trim()
  fa2_new = subset(fa2_new, header_ids2 %in% bfinal$id) # subset for reads that were blasted
  ShortRead::writeFasta(fa2_new, file = paste0(out_path, sample, '-new_2.fa'))
  headers2 = subset(headers2, header_ids2 %in% bfinal$id)
  header_ids2 = subset(header_ids2, header_ids2 %in% bfinal$id)
  df = data.frame(header = headers2, id = header_ids2) %>% left_join(bfinal %>% select(id, staxids, sacc, pos), by = 'id')
  new_headers = paste0(">", df$header, ' | PRISM | staxids:', df$staxids, ' sacc:', df$sacc, ' pos:', df$pos)
  write(new_headers, file = paste0(out_path_final, 'new_headers.txt'))
  str = paste0("awk 'NR==FNR { h[++i] = $0; next } /^>/ { print h[++j]; next } { print }' ", out_path_final, 'new_headers.txt ', 
               out_path, sample, '-new_2.fa > ', out_path_final, sample, '_2.fa')
  system(str)
  str = paste0('rm ', out_path, sample, '_2.fa ', out_path, sample, '-new_2.fa ', out_path_final, 'new_headers.txt')
  system(str)
}


cat(paste('Finished at', Sys.time(), '\n\n'))





