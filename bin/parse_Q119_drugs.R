# parse Q119 into drug names and therapeutic classes
# devtools::install_github("nt-williams/rxnorm")
# devtools::install_github("mpancia/RxNormR")
library(tidyverse)
library(stringr)
library(english)
library(tm)
library(topicmodels)
library(lattice)
library(tidytext)
library(rxnorm)
library(RxNormR)
data(stop_words)

# words to filter out
vetmed_stopwords = c("vet","veterinarian","veterinarians",
                     "0","1","2","3","4","5","6","7","8","9","10",
                     "na","N/A","NA","none","no",
                     "dosage","dose","doses",
                     "mg","MG","ng","NG","ml","ML",
                     "1x","2x","3x","4x","5x",
                     "day","days","daily",
                     "week","weeks","weekly","biweekly","weekends",
                     "month","months","monthl","monthly",
                     "year","years","yearly",
                     "medication","medications","meds","medicine",
                     "drug","drugs",
                     "treatment","treatments",
                     "preventive","preventatives",
                     "preventative","preventatives",
                     "prevantative","prevantatives",
                     "prevention","preventions",
                     "repellent","spray","topical","cream","powder",
                     "pill","pills",
                     "tablet","tablets",
                     "syringe","syringes",
                     "drop","drops",
                     "approx","approximately",
                     "scheduled","regime",
                     "cycles","combination","cooperative",
                     "recurrent","chronic",
                     "allergy","allergies","springtime","seasonal",
                     "releif","relief","helped",
                     "holistic","homeopathic",
                     "nex","oster","chloe","plaqu","gard","guard",
                     "heartworm","hookworm","worm",
                     "tick","ticks",
                     "flea","fleas",
                     "fish","oil",
                     "supplement","supplements",
                     "vitamin","vitamins",
                     "probiotic","probiotics",
                     "heart","joint","joints","eye","eyes","ear","ears",
                     "dog","dogs","human","humans","people",
                     "OTC","over-the-counter","counter",
                     "green","greenies",
                     "CBD","cannabis",
                     "k9","K9","collar",
                     "comfort","laser","hurts","pain","foster","question","finate","crazy","daycare","greenie","isn't","adopted","advocate")

med = answers %>%
  filter(question %in% c(119)) %>%
  select(dog,question,answer) %>%
  mutate(answer = gsub("([a-z])([A-Z])","\\1 \\2",answer)) %>%
  unnest_tokens(word, answer) %>%
  count(word, sort = TRUE) %>%
  filter(!word %in% stop_words$word) %>%
  filter(!word %in% vetmed_stopwords) %>%
  filter(nchar(word) > 4) %>%
  filter(!str_detect(word, "[0-9]")) %>%
  filter(n > 1) %>%
  arrange(-n)

medications$search_rxcui = unlist(lapply(lapply(X = medications$word, FUN = function (x) {RxNormR::rx_approximateTerm(x)$approximateGroup$candidate[[1]]$rxcui}), function(x) if (length(x) == 0) 0 else x))

medications = medications %>%
  filter(search_rxcui != 0) %>%
  mutate(generic_rxcui = purrr::map_chr(.x = search_rxcui,
                                        .f = function(x) paste(rx_related_tty(x, "IN")[[1]]$conceptGroup[[1]]$conceptProperties[[1]]$rxcui, collapse = "")),
         brand_rxcui = purrr::map_chr(.x = search_rxcui,
                                      .f = function(x) paste(rx_related_tty(x, "BN")[[1]]$conceptGroup[[1]]$conceptProperties[[1]]$rxcui, collapse = "")))

medications = medications %>%
  rowwise() %>%
  mutate(therapeutic = paste(rxnorm::get_atc(generic_rxcui, "second"), collapse = "|")) %>%
  mutate(generic_name = purrr::map_chr(.x = generic_rxcui,
                                       .f = function(x) get_rx(x)))

medications = medications %>%
  filter(!is.na(therapeutic) & therapeutic != "NA" & !is.na(generic_name)) %>%
  select(word,n,rxcui=generic_rxcui,rxname=generic_name,rxatc=therapeutic)

medications = medications %>%
  separate(rxatc, into = paste("tx",1:10, sep = ""), sep = "\\|") %>%
  pivot_longer(cols = paste("tx",1:10, sep = ""),
               names_to = "txout",
               values_to = "rxatc",
               values_drop_na = T) %>%
  filter(rxatc != "NA") %>%
  select(word,n,rxcui,rxname,rxatc)


med_key = answers %>%
  filter(question %in% c(119)) %>%
  select(dog,question,answer) %>%
  mutate(answer = gsub("([a-z])([A-Z])","\1 \2",answer)) %>%
  unnest_tokens(word, answer) %>%
  filter(!word %in% stop_words$word) %>%
  filter(!word %in% med_stopwords) %>%
  filter(!word %like% "mg") %>%
  filter(nchar(word) > 4) %>%
  filter(!str_detect(word, "[0-9]")) %>%
  filter(word %in% (medications %>%
                      filter(rxatc %in% tx_include | rxname %in% px_include) %>%
                      pull(word)))

med_key = medications %>% merge(med_key, by = "word")

dog_meds = answers %>%
  filter(question == 119) %>%
  select(dog,answer,notes) %>%
  merge((med_key %>%
           select(dog,word,rxcui,rxname,rxatc)), by = "dog", all = T) %>%
  arrange(dog,rxname) %>%
  filter(!is.na(rxname)) %>%
  unique() %>%
  as.data.table()
