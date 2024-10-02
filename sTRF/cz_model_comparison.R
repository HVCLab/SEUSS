# obtain p values by participant for each model comparing against permuted null distribution

library(tidyverse)


# channels of interest
chanind = R.matlab::readMat('../util/channel_ind_mtr.mat') 
ind = chanind[["chanindmtr"]][[1]]
mag = chanind[["chanindmtr"]][[2]]

 roi = ind[1:10]
 # roi = 1:124


# permuted distribution, Careful, you should only keep Xmtr elec, same as you do
aux1 = R.matlab::readMat('./mtrPermDistr_model1_091224.mat') %>% as.data.frame() %>%
  filter(complete.cases(.)) %>% select(all_of(roi))
aux2= R.matlab::readMat('./mtrPermDistr_model2_091224.mat') %>% as.data.frame() %>%
  filter(complete.cases(.)) %>% filter(complete.cases(.)) %>% select(all_of(roi))
aux3 = R.matlab::readMat('./mtrPermDistr_model3_091224.mat') %>% as.data.frame() %>%
  filter(complete.cases(.)) %>%  filter(complete.cases(.)) %>% select(all_of(roi))
aux4 = R.matlab::readMat('./mtrPermDistr_model4_091224.mat') %>% as.data.frame() %>%
  filter(complete.cases(.)) %>%  filter(complete.cases(.)) %>% select(all_of(roi))
aux5 = R.matlab::readMat('./mtrPermDistr_model5_091224.mat') %>% as.data.frame() %>%
  filter(complete.cases(.)) %>%  filter(complete.cases(.)) %>% select(all_of(roi))
 aux6 = R.matlab::readMat('./mtrPermDistr_model6_091224.mat') %>% as.data.frame() %>%
   filter(complete.cases(.)) %>%  filter(complete.cases(.)) %>% select(all_of(roi))

aux = list(aux1,aux2,aux2,aux4,aux5,aux6)
(thre = aux %>% map_dbl(\(aux) quantile(aux, 1, probs = 0.95)))

## load data
datos =  R.matlab::readMat('./modelcomparison_fitAllchan_091224_all.mat')

# average across ALL channels 
modelSbj = as.data.frame(apply(datos$fitsAllchan[, roi, ], c(1, 3), mean))
modelGen = as.data.frame(apply(datos$fitsAllchanGen[, roi, ], c(1, 3), mean))


modelGen$modelo = 'gen'
modelSbj$modelo = 'sbj'
modelGen$subjectNmb = 1:26
modelSbj$subjectNmb = 1:26

nombres = c('m1:SO','m2:SO+PR','m3:SO+S',
            'm4:SO+R','m5:SO+R+S','m6:SO+R*S',
            'modelo','sbjNmb')


names(modelGen) = nombres
names(modelSbj) = nombres
model = rbind(modelSbj,modelGen)


model_long = model %>% 
  pivot_longer(c(-modelo,-sbjNmb),
               names_to = 'modelNmb',
               values_to = 'mtr')
model_long$R2 = model_long$mtr*model_long$mtr

# add p95 threshold for each model
model_long = model_long %>% 
  mutate(thre = recode(modelNmb, !!!setNames(thre, nombres[1:6])))

# inspect
# 
model_long %>% 
  group_by(modelNmb, modelo) %>% 
  ggplot(aes(modelNmb,mtr*mtr, col = modelNmb)) +
  geom_jitter(size = 4,alpha = 0.3,width = 0.1) + 
  geom_point(stat = 'summary', fun.data =   mean_se, size = 6, col = 'grey30') +
  geom_errorbar(stat = 'summary',fun.data =   mean_se,size = 1, col = 'grey30', width = 0.1) +
  facet_grid(~modelo) +
  geom_point(aes(y = thre*thre), size = 4, col = 'red', shape = 18) +
  # ylim(0,0.011) +
  theme_minimal() +
  xlab('model number') + ylab('R2') + ggtitle('R squared')


model_long %>%
  group_by(modelNmb, modelo) %>%
  ggplot(aes(modelNmb,mtr, col = modelNmb)) +
  geom_jitter(size = 4,alpha = 0.3,width = 0.1) +
  geom_point(stat = 'summary', fun.data =   mean_se, size = 6, col = 'grey30') +
  geom_errorbar(stat = 'summary',fun.data =   mean_se,size = 1, col = 'grey30', width = 0.1) +
  facet_grid(~modelo) +
  geom_point(aes(y = thre), size = 4, col = 'red', shape = 18) +
   # ylim(0,0.11) +
  theme_minimal() +
  xlab('model number') + ylab('r') + ggtitle('Pearsons r') 
# 
# # 
# # # Calculate the mean values for each sbjNmb
# # means <- model_long %>%
# #   group_by(sbjNmb) %>%
# #   summarize(mean_value = mean(mtr))
# # 
# # # Reorder sbjNmb based on the mean values
# # model_long <- model_long %>%
# #   mutate(sbjNmb = factor(sbjNmb, levels = means$sbjNmb[order(means$mean_value)]))
# # 
# # 
# # model_long %>% 
# #   filter(modelNmb != 'model1') %>% 
# #   group_by(modelNmb,modelo,sbjNmb) %>% 
# #   ggplot(aes(sbjNmb,R2, col = modelNmb)) +
# #   geom_point(size = 4,alpha = 0.6,width = 0.1) + 
# #   facet_grid(~modelo) +
# #   # ylim(0,0.011) +
# #   theme_minimal() +
# #   xlab('sbj nmb') +  ylab('R2') + ggtitle('R squared')
# # 
# # model_long %>% 
# #   filter(modelNmb != 'model1') %>% 
# #   group_by(modelNmb,modelo,sbjNmb) %>% 
# #   ggplot(aes(sbjNmb,mtr, col = modelNmb)) +
# #   geom_point(size = 4,alpha = 0.6,width = 0.1) + 
# #   facet_grid(~modelo) +
# #   # ylim(0,0.11) +
# #   theme_minimal() +
# #   xlab('sbj nmb') +  ylab('r') + ggtitle('Pearson r')
# # 
# # # obtain means by model type and model number
# # model_long %>% 
# #   group_by(modelo,modelNmb) %>% 
# #   summarise(across(c(mtr,R2),list(media = mean,desvio = sd)))
# 
# # pairwise model comparison, friedman
# # model_long %>%
# #   filter(modelo == 'sbj',modelNmb != 'model1') %>% 
# #   rstatix::friedman_test(R2 ~ modelNmb | sbjNmb)

model_long %>%
  filter(modelo == 'gen',modelNmb != 'm1:SO') %>% 
  rstatix::friedman_test(R2 ~ modelNmb | sbjNmb)



# pairwise model comparison, wilcox
model_long %>%
  filter(modelo == 'sbj',modelNmb != 'm1:SO') %>%
  rstatix::wilcox_test(R2 ~ modelNmb,paired = T,p.adjust = 'fdr')

# all against all
model_long %>%
  filter(modelo == 'gen',modelNmb != 'm1:SO') %>%
  rstatix::wilcox_test(R2 ~ modelNmb,paired = T,p.adjust = 'fdr')

# against model 2
model_long %>%
  filter(modelo == 'gen',modelNmb %in% c('m3:SO+S','m4:SO+R','m2:SO+PR','m5:SO+R+S')) %>%
  rstatix::wilcox_test(R2 ~ modelNmb,paired = T,p.adjust = 'fdr',
                       ref.group = 'm2:SO+PR')

# 3 against 4
model_long %>%
  filter(modelo == 'gen',modelNmb %in% c('m3:SO+S','m4:SO+R')) %>%
  rstatix::wilcox_test(R2 ~ modelNmb,paired = T,p.adjust = 'fdr')

# 3 and 4 against 5
# 
model_long %>%
  filter(modelo == 'gen',modelNmb %in% c('m3:SO+S','m4:SO+R','m5:SO+R+S')) %>%
  rstatix::wilcox_test(R2 ~ modelNmb,paired = T,p.adjust = 'fdr',
                       ref.group = 'm5:SO+R+S')



stat.test = model_long %>%
  filter(modelo == 'gen',modelNmb != 'm1:SO') %>% 
  rstatix::wilcox_test(R2 ~ modelNmb,paired = T,
                       p.adjust = 'fdr')

library(rstatix)
library(ggpubr)
p <-  model_long %>% 
  group_by(sbjNmb) %>% 
  filter(modelo == 'sbj',modelNmb != 'm1:SO') %>% 
  ggpaired(x = 'modelNmb', y = 'R2', id = 'sbjNmb',color = 'modelNmb', palette = "jco", 
  line.color = "gray", line.size = 0.1)

p + stat_pvalue_manual(stat.test, label = "p.adj", y.position = c(0.0042,0.0044, 0.0046,0.0048))

# Yulias model diagnostics, by participant
# everything against model2
# 
# 
model_longm2 = model %>% 
  select(-`m1:SO`) %>% 
  filter(modelo == 'gen') %>% 
  pivot_longer(-c(`m2:SO+PR`,modelo,sbjNmb),
               names_to = 'modelNmb',
               values_to = 'mtr')
model_longm2 %>% 
  ggplot(aes(`m2:SO+PR`^2, mtr^2, col = as.factor(sbjNmb)), ) +
  geom_abline(slope = 1, intercept = 0, col = 'gray70',size = 1) +
  geom_point(size = 6, alpha = 0.7) +
  theme_minimal() +
  facet_wrap(~ modelNmb,nrow = 2) 

### CHANNEL TOPOGRAPHIC ANALYSIS ###

# pairwise model comparison, wilcox
# model_long %>%
#   filter(modelo == 'sbj',modelNmb != 'm1:SO') %>% 
#   rstatix::wilcox_test(mtr ~ modelNmb,paired = T,p.adjust = 'fdr')


# 
# 
# wiclox comparison of channel across subjects
# average across ALL subjects
datosCan =  R.matlab::readMat('./modelcomparison_fitAllchan_091224_all.mat')



# take 3d matrix and turn it into a long df
modelGenChan_long <- lapply(1:6, function(i) {
  model_df = as.data.frame(datosCan$fitsAllchanGen[,,i])
  names(model_df) = 1:124
  model_df = model_df %>%
    pivot_longer(cols = everything(),
                 names_to = 'chan',
                 values_to = 'mtr') %>%
    mutate(model = paste0("model", i)) # Add dynamic model name
  return(model_df)
}) %>%
  bind_rows() 


# pairwise model comparison, by channel

salida24 = list()  # Initialize salida as an empty list
n = length(levels(as.factor(modelGenChan_long$chan)))  # Define 'n' as the number of unique channels

for (i in 1:n) {
  result = modelGenChan_long %>%
    mutate(R2 = mtr^2) %>%
    group_by(chan) %>%
    filter(model %in% c('model2', 'model4')) %>%
    filter(chan %in% levels(as.factor(modelGenChan_long$chan))[i]) %>%
    rstatix::wilcox_test(R2 ~ model, paired = TRUE, p.adjust.method = 'fdr') %>%
    mutate(chan = levels(as.factor(modelGenChan_long$chan))[i])  # Add channel info
  
  salida24[[i]] = result  # Store the result for each iteration in the list
}


salida23_df = bind_rows(salida) %>%
  mutate(chan =  as.numeric(chan)) %>%
  arrange(chan) 
salida23_df$padj = p.adjust(salida23_df$p,
                            method = 'fdr', n = length(salida23_df$p))

ggplot(salida23_df,aes(1:124,padj)) +
  geom_point(size = 4) +
  ylim(0,0.05) +
  geom_text(aes(label = chan), vjust = -1, size = 3) 

