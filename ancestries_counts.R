# Create a GDS of my own

# Load packages
library(tidyverse)
library(gdsfmt)
library(SNPRelate)

library(sqldf) # For conditional join

# Counts of each ancestry in REDS-III
## AFR=0; EUR=1; NAM=2
# Reading data with LAI (Gnomix output)
lai_inferences <-
  read_delim("/home/chpassos/Analises/REDS-III/ancestries/lai/chr11_output/query_results.msp",
             skip = 1)

# Saving info
info <- lai_inferences |>
  select(c(1:6))

# Removing the first six columns
lai <- lai_inferences |>
  select(-c(1:6))

# Columns that has .0 in its name
dot0 <- data.frame(names(lai)) |>
  mutate(col_n = row_number()) |>
  filter(grepl("\\.0", names.lai.))
# Columns that has .1 in its name
dot1 <- data.frame(names(lai)) |>
  mutate(col_n = row_number()) |>
  filter(grepl("\\.1", names.lai.))

# Renaming columns
df_dot0 <- lai |>
  select(any_of(dot0$names.lai.)) |>
  rename_with(~ str_remove(., "\\.0"), everything()) 
id_order <- colnames(df_dot0) # Separating column names, to order the DF
df_dot0 <- df_dot0 |>
  select(order(id_order)) # Rearranging the column order of dot0
# Renaming columns
df_dot1 <- lai |>
  select(any_of(dot1$names.lai.)) |>
  rename_with(~ str_remove(., "\\.1"), everything())
df_dot1 <- df_dot1 |>
  select(order(id_order)) # Rearranging the column order of dot1

## African Ancestry (AFR=0)
### Checking if all rows, from all columns, are equal to 0 (in this case, AFR ancestry)
afr_dot0 <- df_dot0 |>
  map(~ . == 0) |>
  map(~ as.double(.))
afr_dot1 <- df_dot1 |>
  map(~ . == 0) |>
  map(~ as.double(.))

## European Ancestry (EUR=1)
eur_dot0 <- df_dot0 |>
  map(~ . == 1) |>
  map(~ as.double(.))
eur_dot1 <- df_dot1 |>
  map(~ . == 1) |>
  map(~ as.double(.))

## Native American Ancestry (NAM=2)
nam_dot0 <- df_dot0 |>
  map(~ . == 2) |>
  map(~ as.double(.))
nam_dot1 <- df_dot1 |>
  map(~ . == 2) |>
  map(~ as.double(.))

# Summing AFR counts
sum_afr <- map2(afr_dot0, afr_dot1, 
                \(afr_dot0, afr_dot1) afr_dot0 + afr_dot1)
# Summing EUR counts
sum_eur <- map2(eur_dot0, eur_dot1, 
                \(eur_dot0, eur_dot1) eur_dot0 + eur_dot1)
# Summing NAM counts
sum_nam <- map2(nam_dot0, nam_dot1, 
                \(nam_dot0, nam_dot1) nam_dot0 + nam_dot1)

# Putting info columns back in the data
## African Ancestry
tbl_afr <- as_tibble(sum_afr)
afr_count <- info |>
  cbind(tbl_afr) |>
  as_tibble()

## European Ancestry
tbl_eur <- as_tibble(sum_eur)
eur_count <- info |>
  cbind(tbl_eur) |>
  as_tibble()

## Native American Ancestry
tbl_nam <- as_tibble(sum_nam)
nam_count <- info |>
  cbind(tbl_nam) |>
  as_tibble()


# What happens now is that we have the ancestries count for EACH ANCESTRY
# But only in the INTERVAL estimated by Gnomix
# I need to get the ancestry count per position
## My choice of doing first by intervals was good because otherwise the data
## would have been very big to work with.
## Now I need to do a conditional join of each ancestry count tbl
## to a file that has the positions in the vcf data :)
### I tried to fit all the positions that are in VCF, but 
### it was taking wayyy too long
### I will need to repeat the same task 3 times (each ancestry), 
### for each chromosome (22).
# So because of this, I will introduce some things below


## Get positions in VCF file used in Gnomix
# bcftools query -f '%POS\n' REDS_chr11_gnomix.vcf.gz > chr11_pos.txt

# Read this file
pos <- 
  read.table("/home/chpassos/Analises/REDS-III/ancestries/lai/chr11_output/chr11_pos.txt") |>
  mutate(chr = "11") |>
  rename("position" = "V1")

# Files of interest
head(pos)
head(afr_count)
head(eur_count)
head(nam_count)

# Conditional join
# I need to divide my dataset - It has so many columns!
afr_count_d1 <- afr_count |> select(1:1376)
afr_part1_joined <- sqldf('SELECT * FROM afr_count_d1 LEFT JOIN pos ON position BETWEEN spos and epos')

afr_count_d2 <- afr_count |> select(c(1:6, 1377:ncol(afr_count)))
afr_part2_joined <- sqldf('SELECT * FROM afr_count_d2 LEFT JOIN pos ON position BETWEEN spos and epos')

  


### European Ancestry
tidy_eur <- eur_count |>
  mutate(row_id = row_number()) |>
  pivot_longer(-c(`#chm`, "spos", "epos", "sgpos", "egpos", `n snps`, "row_id"),
               names_to = "sample_id",
               values_to = "count") |>
  rename("chromosome" = `#chm`,
         "n_snps" = `n snps`)

joined_eur <- sqldf('SELECT chromosome, spos, epos, sample_id, count 
      FROM tidy_eur 
      LEFT JOIN new_pos ON V1 BETWEEN spos and epos')

### Native American Ancestry
tidy_nam <- nam_count |>
  mutate(row_id = row_number()) |>
  pivot_longer(-c(`#chm`, "spos", "epos", "sgpos", "egpos", `n snps`, "row_id"),
               names_to = "sample_id",
               values_to = "count") |>
  rename("chromosome" = `#chm`,
         "n_snps" = `n snps`)

joined_nam <- sqldf('SELECT chromosome, spos, epos, sample_id, count 
      FROM tidy_nam 
      LEFT JOIN new_pos ON V1 BETWEEN spos and epos')


# With the conditional table, now return to wide format
joined_afr |>
  pivot_wider(names_from = "sample_id",
              values_from = "count")


joined_eur
joined_nam




snpgdsCreateGeno(
  # Numeric Matrix for SNP genotypes
)
