
N = seq(2,6)
M=length(Chains$alpha)
indices = seq(1,M)



get_risk_by_household_size <- function(k, Chains){

    a <- Chains$alpha[k]
  A = data.frame(size = N, risk = 1/(N-1)^a)
  return(A)
  
}
input.risk = data.frame(size = N, risk = 1/(N-1)^0.5)
indices = seq(1,M)

g7 = indices %>%
  map(get_risk_by_household_size, Chains= Chains) %>%
  bind_rows() %>%
  group_by(size)%>%
  summarise_at(.vars = "risk",
               .funs = c(mean="mean",quantile025 = "quantile025", quantile975="quantile975")) %>%
  ggplot()+
  geom_line(aes(x= size, y= 100*mean), linetype = 'dashed', color= "red")+
  geom_ribbon(aes(x = size, ymin = 100*quantile025, ymax = 100*quantile975), alpha=0.2, fill= "red") +
  geom_line(data = input.risk, aes(x=size,y = 100*risk),color='black')+ylim(c(0,NA))
print(g7)


# 
# +
#   geom_point(aes(x=size, y = 100*mean.seroprevalence),   color= "black",size=1.6) +
#   scale_color_discrete(name="Cohort")+
#   ylab('Seroprevalence (%)')+
#   ylim(c(0,100))+
#   theme_bw() +
#   xlab('Sampling year')+
#   scale_x_continuous(breaks=seq(1995,2012))+
#   theme(axis.text.x = element_text(size=16,angle = 45, vjust=0.5),
#         axis.text.y = element_text(size=16),
#         text=element_text(size=16))+
# print(g7)




## statistics data

xdata$all.HH
xdata$size_households
xdata$n_infected
 

color.dark.blue = "#246575"
n.simulations = 1
#name = 'French_HH_07122024_'

color.within.host = "#20B4B6"
color.fit = "#FF8A5B"
color.data.exp= "#5C9B58"
color.data.constant = "#EA526F"
#color.5 = "#FCEADE"



D = data.frame(hh.size=xdata$size_households, X =xdata$n_infected, epid.model = "exp")  %>%
  mutate(data = ifelse(epid.model =="exp", "Dose-response", "Constant"))  

b = seq(0,5)

stacked_histogram_hhsize =D %>% group_by(data,hh.size) %>%# %>% group_by(epid.model, hh.size,transmission.model) %>%
  count(X) %>%
  mutate(Infections =X)%>%
  # mutate(Number = factor(X, levels = rev(seq(0,5)))) %>%
  ggplot(aes(fill = Infections,y= n, x =hh.size ))+
  geom_bar(position = "fill", stat = "identity", show.legend = FALSE) +
  theme_bw()+
  theme(axis.text.x = element_text(size=14,angle = 0, vjust=0),
        axis.text.y = element_text(size=14),
        text=element_text(size=14),
        strip.background = element_rect(fill=FALSE))+
  xlab("Household size") +
  ylab('Proportion')+ #
  #+ scale_fill_manual(values = c(color.data.constant,color.data.exp))+
  scale_fill_gradientn(limits = c(0,6),
                       #  colours=c("navyblue", "darkmagenta", "darkorange1"),
                       colours=c(color.dark.blue, color.within.host, "darkorange1"),
                       breaks=b, labels=format(b))
print(stacked_histogram_hhsize)



