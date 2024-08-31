library("ggplot2") 

### Dataset from Henao et al. 2019 PNAS
## Case Study 1
rates_time<-read.csv("~/datasets_cases12/summary_tree_results.csv")
head(rates_time)



## Speciation rates with error
speciation.plot<-ggplot(rates_time, aes(x=tree.max.age, y=mean.clade.lambda)) + geom_point() +geom_linerange(aes(ymin =lambda.lower.se , ymax =lambda.upper.se),colour="gray",alpha=0.5)
speciation.plot<-speciation.plot+theme_classic()+xlab("Tree age")+ylab("Average speciation rate for the tree")

## Extinction rates with error
extinction.plot<-ggplot(rates_time, aes(x=tree.max.age, y=mean.clade.mu)) + geom_point() +geom_linerange(aes(ymin =mu.lower.se , ymax =mu.upper.se),colour="gray",alpha=0.5)
extinction.plot<-extinction.plot+theme_classic()+xlab("Tree age")+ylab("Average speciation rate for the tree")

### Number of tips against uncertainty
p<-ggplot(rates_time, aes(x=mean.clade.lambda, y=ntips))+geom_point()+geom_linerange(aes(xmin =lambda.lower.se , xmax =lambda.upper.se))

#Creating categories by tree size
color_cat<-rep(0,length(rates_time$ntips))
size_cat<-rep(0,length(rates_time$ntips))
lessthan100<-which(rates_time$ntips<100)
bet100300<-which((rates_time$ntips>99.9&rates_time$ntips<300))
bet300500<-which(rates_time$ntips>299.9&rates_time$ntips<500)
morethan500<-which(rates_time$ntips>499.9)
color_cat[lessthan100]="#FFBE0B"
size_cat[lessthan100]=3
color_cat[bet100300]="#3f88c5"
size_cat[bet100300]=10
color_cat[bet300500]="#a5be00"
size_cat[bet300500]=17
color_cat[morethan500]="#d00000"
size_cat[morethan500]=24

# Figure 2- Reproducing Henao-Diaz et al 2019 with uncertainty
rates_time<-cbind(rates_time,size_cat,color_cat)
p2.v2<-ggplot(rates_time,aes(x=tree.max.age, y=mean.clade.lambda)) + geom_point(alpha=0.95,color=color_cat)+geom_linerange(aes(ymin =lambda.lower.se , ymax =lambda.upper.se),alpha=0.5,colour="gray") +theme_classic()+xlab("Tree age")+ylab("Average speciation rate")

### Inset speciation and extinction- Figure 2
specextinct<- ggplot(rates_time, aes(x=mean.clade.lambda, y=mean.clade.mu))+geom_linerange(aes(xmin=lambda.lower.se , xmax =lambda.upper.se), alpha=0.3, colour="gray")+theme_classic()+xlab("Speciation average rate")+ylab("Extinction average rate")
specextinct<-specextinct+geom_linerange(aes(ymin=mu.lower.se, ymax=mu.upper.se),alpha=0.3, colour="gray")+geom_point(alpha=0.5, color=color_cat, size=4)
specextinct+geom_abline(slope=1)


# "Scaling" of standard error by mean= How big are errors?
# This is the summary statistic for relative standard error for estimates of speciation
y=(rates_time$lambda.upper.se- rates_time$lambda.lower.se)/rates_time$mean.clade.lambda
rates_time<-cbind(rates_time,y)

rangevtree<-ggplot(rates_time, aes(x=tree.max.age, y=y))+
  geom_point(color=color_cat)+theme_bw()+xlab("Tree Age")+
  geom_smooth(aes(color=color_cat), method="lm", se=F) +
  ylab("Relative range of standard error")+
  geom_smooth(aes(tree.max.age,y,colour=color_cat,)) 


# A simpler model investigating if clade age and the different sizes of trees are correlatied with relative standard error

linear.model<-lm(y~tree.max.age+size_cat, data=rates_time)
summary(linear.model)


# This is the correct linear model presented in the main manuscript
linear.model2<-lm(y~tree.max.age+ size_cat+tree.max.age*size_cat, data=rates_time)
summary(linear.model2)

# Effect can be see also with number of tips and not only with the categories
linear.model3<-lm(y~tree.max.age+ ntips, data=rates_time)
summary(linear.model3)


#What happens with sampling fraction? Is it informative of uncertainty 
prop.sampled<-rates_time$ntips/rates_time$n.clade

rates_time=cbind(rates_time,prop.sampled)

# This is the response to reviews- Sampling fraction is not informative of relative error-uncertainty.
linear.model4<-lm(y~tree.max.age+ size_cat+tree.max.age*size_cat+prop.sampled, data=rates_time)
summary(linear.model4)


# This is the response to reviews- What happens if we only choose well sampled trees
high_sampled=rates_time%>%filter(prop.sampled>0.75)

# Answer the correlation breaks mostly because we have remove a lot of trees that are older or
# with many tips.
linear.model7<-lm(y~tree.max.age, data=high_sampled)
summary(linear.model7)

#This can be seen in this figure below
ggplot(data=high_sampled, aes(x=tree.max.age,y=y))+
  geom_point(color=cols2)+
  theme_bw()+
  xlab("Tree age")+
  ylab("y=Relative range of standard error")+
  geom_smooth(method="lm",se=T)


# For those highly sampled, there is an exponential decay of uncertainty when we increase the
# number of tips
linear.model9<-lm(log(y)~ntips, data=high_sampled)
summary(linear.model9)

log.model.df <- data.frame(x = high_sampled$ntips,y = exp(fitted(linear.model9)))
rangevtre3<-ggplot(high_sampled, aes(x=ntips, y=y))+
  geom_point(color=cols2)+
  theme_bw()+
  geom_line(data = log.model.df, aes(x, y, color = "Log Model"))  + 
  xlab("Number of tips in the phylogeny")+
  ylab("Relative range of standard error")


# Figure 3- Sampling fraction not correlated with uncertainty 
rangevtre5<-ggplot(rates_time, aes(x=prop.sampled, y=y))+
  geom_point(color=color_cat)+
  geom_smooth(method="lm", se=T) +
  theme_bw()+
  xlab("Sampling fraction")+
  ylab("Relative range of standard error")



##########Distribution by paeleontological period

library("cowplot")

### Paleocene trees
paleocene<-which((rates_time$tree.max.age< 60)&(rates_time$tree.max.age>23))
length(paleocene)

paleocene.dataframe<-data.frame(cbind(rates_time[paleocene,],Era=rep("Paleocene",length(paleocene))))

ggplot(paleocene.dataframe, aes(x=mean.clade.lambda))+geom_density()

## Neogene trees
neogene<-which(rates_time$tree.max.age< 23)
length(neogene)

neogene.dataframe<-data.frame(cbind(rates_time[neogene,],Era=rep("Neogene",length(neogene))))

### Densities by era, what is slow or fast?
era.dataframe<-rbind(paleocene.dataframe,neogene.dataframe)
colors_two<-c("#D90429","#2b2d42")
p4<-ggplot(era.dataframe, aes(x=mean.clade.lambda,color=Era,fill=Era))+geom_density(alpha=0.5)+theme_classic() + scale_fill_manual(values=colors_two)+scale_color_manual(values=colors_two) +xlab("Average speciation rate")+ylab("Density from estimated rates")
ggplot(era.dataframe, aes(x=mean.clade.lambda,color=era))+geom_boxplot()

p4<-ggplot(era.dataframe, aes(x=mean.clade.lambda,color=era,fill=era))+geom_histogram(alpha=0.5)+theme_classic()+scale_colour_manual(values=c("red","blue")) + xlab("Average speciation rate")+ylab("Density")

###### Investigating bradytely or tachytely with point estimaes
##### PALEOCENE ###########
quantile(paleocene.dataframe$mean.clade.lambda, probs=c(0.025,0.975))

#### Paleocene  Slow or fast speciation cut offs
#2.5%    97.5% 
#0.062714 0.723193

# Examples of fast speciation Paleocene
which(paleocene.dataframe$mean.clade.lambda>0.723193)
##
paleocene.dataframe$study.ref[c(17,40)]

## "Neupane, S., Lewis, P. O., Dessein, S., Shanks, H., Paudyal, S. and Lens, F. (2017), Evolution of woody life form on tropical mountains in the tribe Spermacoceae (Rubiaceae). American Journal of Botany, 104: 419\xd0438. doi:10.3732/ajb.1600248"
## "Uribe-Convers S, Tank DC (2015) Shifts in diversification rates linked to biogeographic movement into new areas: an example of a recent radiation in the Andes. American Journal of Botany 102(11): 1854-1869. https://doi.org/10.3732/ajb.1500229" 


# Examples of slow speciation Paleocene
which(paleocene.dataframe$mean.clade.lambda<0.062714)
##
paleocene.dataframe$study.ref[c( 3,39)]
# "Thijmen Breeschoten, Camiel Doorenweerd, Sergei Tarasov, Alfried P. Vogler, Phylogenetics and biogeography of the dung beetle genus Onthophagus inferred from mitochondrial genomes, In Molecular Phylogenetics and Evolution, Volume 105, 2016, Pages 86-95, ISSN 1055-7903, https://doi.org/10.1016/j.ympev.2016.08.016."
# "Weston Testo, Michael Sundue, A 4000-species dataset provides new insight into the evolution of ferns, Molecular Phylogenetics and Evolution, Volume 105, December 2016, Pages 200-211, ISSN 1055-7903, https://doi.org/10.1016/j.ympev.2016.09.003."   


### NEOGENE ###
quantile(neogene.dataframe$mean.clade.lambda, probs=c(0.025,0.975))
### Neogene slow or fast speciation cut offs 
##2.5%    97.5% 
##  0.127091 1.365704  

# Examples of fast speciation neogene
which(neogene.dataframe$mean.clade.lambda>1.365704)
##
neogene.dataframe$study.ref[5]
# "Lagomarsino, L. P., Forrestel, E. J., Muchhala, N. and Davis, C. C. (2017), Repeated evolution of vertebrate pollination syndromes in a recently diverged Andean plant clade. Evolution, 71: 1970-1985. doi:10.1111/evo.13297"

# Examples of slow speciation neogene
which(neogene.dataframe$mean.clade.lambda<0.127091)
##
neogene.dataframe$study.ref[13]
#Felicien Tosso, Olivier J. Hardy, Jean-Louis Doucet, Kasso Danou, Esra Kaymak, Jeremy Migliore, Evolution in the Amphi-Atlantic tropical genus Guibourtia (Fabaceae, Detarioideae), combining NGS phylogeny and morphology, Molecular Phylogenetics and Evolution, Volume 120,2018,Pages 83-93,https://doi.org/10.1016/j.ympev.2017.11.026."


########### Extinction #############

##### PALEOCENE ###########
quantile(paleocene.dataframe$mean.clade.mu, probs=c(0.025,0.975))

#### Paleocene  Slow or fast extinction
#2.5%     97.5% 
# 0.0057635 0.5960834

# Examples of fast extinction Paleocene
which(paleocene.dataframe$mean.clade.mu>0.5960834)
##
paleocene.dataframe$study.ref[c(27, 40)]

# "Soares, Andr_ E. R. , Ben J. Novak, James Haile, Tim H. Heupink, Jon Fjelds, M. Thomas P. Gilbert, Hendrik Poinar, George M. Church, Beth Shapiro. 2016. Complete mitochondrial genomes of living and extinct pigeons revise the timing of the columbiform radiation. BMC Evolutionary Biology 16: 230"
# "Uribe-Convers S, Tank DC (2015) Shifts in diversification rates linked to biogeographic movement into new areas: an example of a recent radiation in the Andes. American Journal of Botany 102(11): 1854-1869. https://doi.org/10.3732/ajb.1500229"   

# Examples of slow extinction Paleocene
which(paleocene.dataframe$mean.clade.lambda<0.0057635)
##

##None- this could be an issue with extinction estimates

### NEOGENE ###
quantile(neogene.dataframe$mean.clade.mu, probs=c(0.025,0.975))
### Neogene slow or fast extinction
##  2.5%     97.5% 
#0.0523765 0.6853915   

# Examples of fast extinction in the neogene
which(neogene.dataframe$mean.clade.mu>0.6853915)
##
neogene.dataframe$study.ref[10]
# "Ribas, C. C., A. Aleixo, A. C. R. Nogueira, C. Y. Miyaki, J. Cracraft, 2011. A palaeobiogeographic model for biotic diversification within Amazonia over the past three million years. Proceedings of the Royal Society B: Biological Sciences 279 (1729): 681-689."

# Examples of slow extinction in the neogene
which(neogene.dataframe$mean.clade.lambda<0.0523765)
## none- this could be an issue with extinction estimates

### CASE STUDY 2 
### Dataset from Helmstetter et al. 2022 (SSE review)
sse_review<-read.csv("~/datasets_cases12/helmstetter_review.csv")
head(sse_review)

#Fastest transition rate
which.max(sse_review$transition_rate)
sse_review[c(13:14),]
# maximum rate of state transition 
#nacker_et_al_2011_Origins_and_consequences_of_serpentine_endemism_in 3 endemic to non-endemic transition BiSSE

plot(sse_review$age, sse_review$transition_rate) #same patterns! This is crazy!

boxplot(sse_review$div_rate, na.rm=TRUE, ylim=c(-3,3))
# which are the extremes of net-diversifications
which.max(sse_review$div_rate) #411
sse_review[c(410:412),]  #Frenzke_et_al_2016_Evolution_of_epiphytism_and_fruit_traits_act_uneve.txt BiSSE Piperaceae     transition rate 27.67  ? age sampling 0.07
which.min(sse_review$div_rate) # 16
sse_review[c(15:17),] # Anacker_et_al_2011_Origins_and_consequences_of_serpentine_endemism_in.txt -10 GeoSSE
plot(sse_review$age, sse_review$div_rate, ylim=c(-4,6)) # same patterns though more constancy

#### Creating the same categories as Case 1
tips_sampled<-sse_review$tips*sse_review$perc_sampling
lessthan100<-which(tips_sampled<100)
bet100300<-which((tips_sampled>99.9&tips_sampled<300))
bet300500<-which(tips_sampled>299.9&tips_sampled<500)
morethan500<-which(tips_sampled>499.9)
color_cat=rep(0,length(tips_sampled))
size_cat<-rep(0,length(tips_sampled))
color_cat[lessthan100]="#FFBE0B"
size_cat[lessthan100]<-3
color_cat[bet100300]="#3f88c5"
size_cat[bet100300]<-10
color_cat[bet300500]="#a5be00"
size_cat[bet300500]<-17
color_cat[morethan500]="#d00000"
size_cat[morethan500]<-24
sse_review<-cbind(sse_review,tips_sampled)
sse_review<-cbind(sse_review,size_cat)

# Figures 5 and 6 main manuscript
p5<-ggplot(sse_review,aes(x=age, y=transition_rate,size=size_cat)) +
  geom_point(color=color_cat,alpha=0.5)+ scale_size(name="Number of tips",breaks=c(100, 300, 500, 1000))+
  theme_classic()+xlab("Tree age")+ylab("Transition rates between states")+
  coord_cartesian(xlim=c(0,150))+
  geom_vline(xintercept =c(23,60),lty="dashed") 

p6<-ggplot(sse_review,aes(x=age, y=div_rate,size=size_cat)) +
  geom_point(color=color_cat,alpha=0.5)+ 
  scale_size(name="Number of tips",breaks=c(100, 300, 500, 1000))+
  theme_classic()+xlab("Tree age")+
  ylab("Net diversification rates per state")+coord_cartesian(xlim=c(0,150),ylim=c(-2,3.3))+
  geom_hline(yintercept=0)+geom_vline(xintercept =c(23,60),lty="dashed") 

# Inset of speciation v. extinction
specextinct_Helms<- ggplot(sse_review, aes(x=mean.clade.lambda, y=mean.clade.mu))+geom_point(alpha=0.5, color=color_cat)+geom_linerange(aes(xmin=lambda.lower.se , xmax =lambda.upper.se), alpha=0.3, colour="gray")+theme_classic()+xlab("Speciation average rate")+ylab("Extinction average rate")
specextinct_Helms<-specextinct+geom_linerange(aes(ymin=mu.lower.se, ymax=mu.upper.se),alpha=0.3, colour="gray")

# Checking the amount of variation per sampling fraction, no patterns
plot(sse_review$perc_sampling,sse_review$div_rate) # plenty of variation independently from sampling, a litte broader when sampling <0.5
plot(sse_review$perc_sampling,sse_review$transition_rate)  #plenty of variation independently from sampling,at all levels

### Paleocene trees
paleocene<-which((sse_review$age< 60)&(sse_review$age>23))
length(paleocene) #341

paleocene.dataframe<-data.frame(cbind(sse_review[paleocene,],era=rep("Paleocene",length(paleocene))))

ggplot(paleocene.dataframe, aes(x=div_rate))+geom_density()

## Neogene trees- Checking differences in this example of neogenea and paleocene. 
neogene<-which(sse_review$age< 23)
length(neogene)

neogene.dataframe<-data.frame(cbind(sse_review[neogene,],era=rep("Neogene",length(neogene))))
ggplot(neogene.dataframe, aes(x=div_rate))+geom_density()

### Densities by era, what is slow or fast?
era.dataframe<-rbind(paleocene.dataframe,neogene.dataframe)
ggplot(era.dataframe, aes(x=div_rate,color=era))+geom_density()
ggplot(era.dataframe, aes(x=div_rate,color=era))+geom_boxplot()

p7<-ggplot(era.dataframe, aes(x=div_rate,color=era))+geom_density()+theme_classic()+scale_colour_manual(values=c("red","blue"))+ xlab("Net Diversification")+ylab("Density")

### Densities by era, what is slow or fast in state change?
ggplot(era.dataframe, aes(x=transition_rate,color=era))+geom_density()
ggplot(era.dataframe, aes(x=transition_rate,color=era))+geom_boxplot()


