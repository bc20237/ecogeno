library(devtools)
library(tidyverse)
library(ggplot2)
require(plyr)
ggplot(p1, aes(x=V1, y=V3)) + geom_point() + geom_density_2d()+
  labs(title=".01 initial expression",x="Time", y = "genes expressed, p")+theme_classic() 
ggplot(p5, aes(x=V1, y=V3)) + geom_point() + geom_density_2d()+
  labs(title=".05 initial expression",x="Time", y = "genes expressed, p")+theme_classic() 
ggplot(p45, aes(x=V1, y=V3)) + geom_point() + geom_density_2d()+
  labs(title=".45 initial expression",x="Time", y = "genes expressed, p")+theme_classic() 
ggplot(p55, aes(x=V1, y=V3)) + geom_point() + geom_density_2d()+
  labs(title=".55 initial expression",x="Time", y = "genes expressed, p")+theme_classic() 
ggplot(p99, aes(x=V1, y=V3)) + geom_point() + geom_density_2d()+
  labs(title=".99 initial expression",x="Time", y = "genes expressed, p")+theme_classic() 
ggplot(p100, aes(x=V1, y=V3)) + geom_point() + geom_density_2d()+
  labs(title="1 initial expression",x="Time", y = "genes expressed, p")+theme_classic() 
