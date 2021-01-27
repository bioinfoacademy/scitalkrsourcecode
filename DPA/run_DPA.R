
# Differential proportion analysis of cell populations between conditions

source("diffprop_functions.R");

#############################################
### Run DPA on population counts from TIP ###
#############################################

## Read in file of counts of cells in each population across conditions
obs.counts = as.matrix(read.csv("TIP_population_counts.csv", row.names = 1))
print(obs.counts)

## Run an example using error (p) of 0.1 and with 100,000 iterations
tip.exp <- generateNull(obs.counts, n=100000, p=0.1);     # Generate the null distribution based on sampling

obs.counts/apply(obs.counts, 1, sum)

### P-value tests for Sham vs MI-Day 3
two.class.test(obs.counts, tip.exp, cond.control="Sham", cond.treatment="MI-day 3",to.plot=T);

## Sham vs MI-Day 7
two.class.test(obs.counts, tip.exp, cond.control="Sham", cond.treatment="MI-day 7",to.plot=T)  

## MI-Day3 vs MI-Day 7
two.class.test(obs.counts, tip.exp, cond.control="MI-day 3", cond.treatment="MI-day 7",to.plot=T)  


########################################################
## Get a table of P-values for a range of 'p' values ###
########################################################

res.table.shamVsMID3 = c()
res.table.shamVsMID7 = c()
res.table.MID3VsMID7 = c()

## Go through a series of error probabilities
for (err_prob in c(0.5, 0.4, 0.3, 0.25, 0.2, 0.15, 0.1, 0.05)) {
  tip.exp <- generateNull(obs.counts, n=100000, p=err_prob);
  ## Sham vs MI-day 3
  res.1 = two.class.test(obs.counts, tip.exp, cond.control="Sham", cond.treatment="MI-day 3",to.plot=F)
  res.table.shamVsMID3 = rbind(res.table.shamVsMID3, res.1)
  ## Sham vs MI-day 7
  res.2 = two.class.test(obs.counts, tip.exp, cond.control="Sham", cond.treatment="MI-day 7",to.plot=F)
  res.table.shamVsMID7 = rbind(res.table.shamVsMID7, res.2)
  ## MI-day 3 vs MI-day 7
  res.3 = two.class.test(obs.counts, tip.exp, cond.control="MI-day 3", cond.treatment="MI-day 7",to.plot=F)
  res.table.MID3VsMID7 = rbind(res.table.MID3VsMID7, res.3)
}

rownames(res.table.shamVsMID3) = as.character(c(0.5, 0.4, 0.3, 0.25, 0.2, 0.15, 0.1, 0.05))
rownames(res.table.shamVsMID7) = as.character(c(0.5, 0.4, 0.3, 0.25, 0.2, 0.15, 0.1, 0.05))
rownames(res.table.MID3VsMID7) = as.character(c(0.5, 0.4, 0.3, 0.25, 0.2, 0.15, 0.1, 0.05))

#######################
### Run on GFP data ###
#######################

obs.counts = as.matrix(read.csv("GFP_population_counts.csv", row.names = 1))
print(obs.counts)

all.exp1<-generateNull(obs.counts, n=100000, p=0.1);  # Generate the null distribution based on sampling

obs.counts/apply(obs.counts,1,sum)

### P-value tests comparing the two sham conditions
two.class.test(obs.counts, all.exp1, cond.control="Sham-day 3", cond.treatment="Sham-day 7",to.plot=T); 

## Get a table of P-values for a range of 'p' values
res.table = c()
for (err_prob in c(0.5, 0.4, 0.3, 0.25, 0.2, 0.15, 0.1, 0.05, 0.01, 0.001)) {
  all.exp1<-generateNull(obs.counts, n=10000, p=err_prob);
  res = two.class.test(obs.counts, all.exp1, cond.control="Sham-day 3", cond.treatment="Sham-day 7",to.plot=F)
  res.table = rbind(res.table, res)
}
rownames(res.table) = as.character(c(0.5, 0.4, 0.3, 0.25, 0.2, 0.15, 0.1, 0.05, 0.01, 0.001))

################################################################################
### Now merge the two sham conditions together and compare to the MI samples ###
################################################################################

obs.counts.merged = obs.counts[1, ] + obs.counts[2, ]
obs.counts.merged = rbind(obs.counts.merged, obs.counts[3:4, ])
rownames(obs.counts.merged) = c("Sham", "MI-day 3", "MI-day 7")
print(obs.counts.merged)

# Generate the null distribution based on sampling
gfp.exp <- generateNull(obs.counts.merged, n=100000, p=0.1);  

obs.counts.merged/apply(obs.counts.merged, 1, sum)

## P-value tests comparing sham to MI-day 3
two.class.test(obs.counts.merged, gfp.exp, cond.control="Sham", cond.treatment="MI-day 3", to.plot=T); 

## P-value tests comparing sham to MI-day 7
two.class.test(obs.counts.merged, gfp.exp, cond.control="Sham", cond.treatment="MI-day 7", to.plot=T);

## P-value tests comparing MI-day 3 to MI-day 7
two.class.test(obs.counts.merged, gfp.exp, cond.control="MI-day 3", cond.treatment="MI-day 7", to.plot=T);

gfp.res.table.shamVsMID3 = c()
gfp.res.table.shamVsMID7 = c()
gfp.res.table.MID3VsMID7 = c()

## Go through a series of error probabilities
for (err_prob in c(0.5, 0.4, 0.3, 0.25, 0.2, 0.15, 0.1, 0.05)) {
  gfp.exp<-generateNull(obs.counts.merged, n=100000, p=err_prob);
  ## Sham vs MI-day 3
  res.1 = two.class.test(obs.counts.merged, gfp.exp, cond.control="Sham", cond.treatment="MI-day 3",to.plot=F)
  gfp.res.table.shamVsMID3 = rbind(gfp.res.table.shamVsMID3, res.1)
  ## Sham vs MI-day 7
  res.2 = two.class.test(obs.counts.merged, gfp.exp, cond.control="Sham", cond.treatment="MI-day 7",to.plot=F)
  gfp.res.table.shamVsMID7 = rbind(gfp.res.table.shamVsMID7, res.2)
  ## MI-day 3 vs MI-day 7
  res.3 = two.class.test(obs.counts.merged, gfp.exp, cond.control="MI-day 3", cond.treatment="MI-day 7",to.plot=F)
  gfp.res.table.MID3VsMID7 = rbind(gfp.res.table.MID3VsMID7, res.3)
}

rownames(gfp.res.table.shamVsMID3) = as.character(c(0.5, 0.4, 0.3, 0.25, 0.2, 0.15, 0.1, 0.05))
rownames(gfp.res.table.shamVsMID7) = as.character(c(0.5, 0.4, 0.3, 0.25, 0.2, 0.15, 0.1, 0.05))
rownames(gfp.res.table.MID3VsMID7) = as.character(c(0.5, 0.4, 0.3, 0.25, 0.2, 0.15, 0.1, 0.05))

