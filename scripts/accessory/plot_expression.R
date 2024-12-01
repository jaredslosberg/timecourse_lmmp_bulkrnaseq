ntd_long_sub<- ntd_long_sub %>% group_by(external_gene_name,age) %>% mutate(mn = mean(counts), std_err = sd(counts)/sqrt(n()))

ggplot(ntd_long_sub,aes(x = age_numeric, y = counts)) + geom_point(aes(color = sex)) +
  facet_wrap(~external_gene_name, scales = "free_y") +
  stat_summary(fun = mean, geom = "line",position=position_dodge2(width = .3)) + 
  geom_errorbar(aes(ymin = mn-std_err, ymax = mn+std_err), width =10)
  