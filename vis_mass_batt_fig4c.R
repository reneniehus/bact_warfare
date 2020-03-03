library(tidyverse)
library(R.matlab)
#library(gridExtra)
library(ggthemes)

which_sens <- "THREE" # [Nu TB TA QS] or Z

jj_max <- 5; 
endtime <- 50; 
#
all_n <- 60
top_n <-  4
mut_n <- 46

# loop to load and wrangle the data
#
df_runs <- list()
for (jj in 1:jj_max) {
                load_name <- paste0( 
                                "/Users/reneniehus/Documents/AAUni/1-ASystemsBiology/15.2-ZoologyProject/DPhil/ZoologyProject/ToxinRegulation2015/Results/allevolvePNAS",which_sens, jj , 't',endtime,'_onlyWARR.mat' )
                data <- readMat(load_name)
                
                dataM <- data$WARR.saveS
                dim(dataM)[3]
                df_list <- list()
                for (i in 1:dim(dataM)[3]) {
                                df_temp <- as_tibble( dataM[,,i])
                                #[fnull | fN  UN | fTB  UTB | fTA  UTA | fQS QS]
                                names(df_temp) <- c('fnull','fN','UN','fTB','UTB','fTA','UTA','fQS','QS',
                                                    'tag_pool','tag_top','tox','tox2','fit')
                                df_temp <- df_temp %>% mutate( rank="migrant",
                                                               rank=ifelse( (1:n()<=(top_n+mut_n)) , "mutant" , rank  ),
                                                               rank=ifelse( (1:n()<=top_n) , "top" , rank  ) )
                                
                                df_list[[i]] <-   df_temp    
                }
                df <- bind_rows(df_list , .id="time" )
                #
                #
                df <- df %>% mutate( sensor="none",
                                     sensor=ifelse(fN!=0, "Nu",sensor),
                                     sensor=ifelse(fTB!=0, "TB", sensor),
                                     sensor=ifelse(fTA!=0, "TA" , sensor),
                                     sensor=ifelse(fQS!=0, "QS" , sensor))
                df <- df %>% mutate( time=as.numeric(time))
                
                # put into
                df_runs[[jj]] <- df
                
}
#
df <- bind_rows(df_runs , .id="run") %>% mutate( run=as.numeric(run))


#
df_mean <- df %>% filter(  ) %>%  # rank!="migrant"
                group_by( run,time ) %>% summarise( f_Nu=sum(sensor=="Nu")/n(),
                                                    f_TB=sum(sensor=="TB")/n(),
                                                    f_QS=sum(sensor=="QS")/n()) %>% 
                gather(key,value,-time,-run) %>% ungroup() %>% 
                group_by(time, key ) %>% summarise( mean_value=mean(value) ) %>% 
                spread( key=key, value=mean_value ) %>% mutate( run=1)
# 
(p1 <- df %>% filter(  ) %>%  # rank!="migrant"
                                group_by( run,time ) %>% summarise( f_Nu=sum(sensor=="Nu")/n(),
                                                                    f_TB=sum(sensor=="TB")/n(),
                                                                    f_QS=sum(sensor=="QS")/n()) %>% 
                                #gather(key,value,-time,-run) %>% 
                                ggplot( aes(x=time,group=run) ) + 
                                geom_line(aes(y=f_Nu),color="blue",alpha=0.3) +
                                geom_line(aes(y=f_TB),color="red",alpha=0.3) +
                                geom_line(aes(y=f_QS),color="orange", alpha=0.3) ) 
(p2 <- p1 + geom_line( data=df_mean , aes(x=time,y=f_Nu),color="blue" , size=1.5 )  +
                                geom_line( data=df_mean , aes(x=time,y=f_TB),color="red" , size=1.5 )  +
                                geom_line( data=df_mean , aes(x=time,y=f_QS),color="orange" , size=1.5 )+
                                geom_vline(xintercept = 21))

p2 + theme_solid()
