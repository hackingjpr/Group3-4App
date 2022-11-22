y <- readRDS(file = "./y.vals.rds")
ecdf(y) -> model
model(y) -> y2


plot( y2, y, pch=0, type="l", col="black", yaxt="n",  xlab  = "frequency",ylab="G3/G4 Score", axes=F)
rect(xleft=0, xright=0.68, ybottom= 0, ytop=1, col="green", border = NA )
rect(xleft=0.68, xright=1, ybottom= 0, ytop=1, col="yellow", border = NA)
lines( y2, y, lwd = 1.5 )
axis(side=2, at=c(-0.1,0,0.2,0.4,0.6,0.8,1,1.1))
axis(side=1, at=c(-0.1,0,0.2,0.4,0.6,0.8,1,1.1))
abline(h=c(0,0.5,1), lty = 2)
abline(v=0.68, lty = 1)


generate_figure_highlight_ecrt <- function(new.sample.meta.score, indexRow) {
  if(is.null(indexRow)){indexRow=1}
  temp.df <- readRDS(file = "./y.vals.rds")
  ecdf(y) -> model
  model(y) -> y2
  plot( y2, y, pch=0, type="l", col="black", yaxt="n",  xlab  = "frequency",ylab="G3/G4 Score", axes=F)
  rect(xleft=0, xright=0.68, ybottom= 0, ytop=1, col="green", border = NA )
  rect(xleft=0.68, xright=1, ybottom= 0, ytop=1, col="yellow", border = NA)
  lines( y2, y, lwd = 1.5 )
  axis(side=2, at=c(-0.1,0,0.2,0.4,0.6,0.8,1,1.1))
  axis(side=1, at=c(-0.1,0,0.2,0.4,0.6,0.8,1,1.1))
  abline(h=c(0,0.5,1), lty = 2)
  abline(v=0.68, lty = 1)


  
  df.lines.hor <-
    foreach(i = 1:length(new.sample.meta.score),
            .combine = rbind) %do% {
              data.frame(
                x = 0,
                xend = max(which(
                  y2 < new.sample.meta.score[i]
                )),
                y = new.sample.meta.score[i],
                yend = new.sample.meta.score[i]
              )
            }
  df.lines.hor$labels <- names(new.sample.meta.score)
  
  df.lines.ver <-
    foreach(i = 1:length(new.sample.meta.score),
            .combine = rbind) %do% {
              data.frame(
                x = max(which(
                  y2 < new.sample.meta.score[i]
                )),
                xend = max(which(
                  y2 < new.sample.meta.score[i]
                )),
                y = new.sample.meta.score[i],
                yend = min(y2)
              )
            }
  df.lines.ver$perc <-
    paste0(round(df.lines.ver$xend / length(y2) * 100), "th")
  
  df.lines.ver$colour <- factor(ifelse(1:nrow(df.lines.ver)==indexRow,"highlight","no.highlight"), levels = c("highlight","no.highlight"))
  df.lines.hor$colour <- factor(ifelse(1:nrow(df.lines.hor)==indexRow,"highlight","no.highlight"), levels = c("highlight","no.highlight"))
  
  b <- b +
    geom_segment(
      aes(
        x = x,
        y = y,
        xend = xend,
        yend = yend,
        colour = as.character(colour)
      ),
      #colour = "red",
      linetype = "dashed",
      data = df.lines.hor
    ) +
    geom_segment(
      aes(
        x = x,
        y = y,
        xend = xend,
        yend = yend,
        colour = colour
      ),
      #colour = "red",
      linetype = "dashed",
      data = df.lines.ver
    ) +
    #geom_text_repel(aes(x = x+10, y = y+0.1, label = labels), direction = "y", data = df.lines.hor)
    geom_text(aes(
      x = x + 10,
      y = y + 0.1,
      label = labels,
      colour = colour
    ), data = df.lines.hor)
  #scale_color_discrete("red", "lightgrey")

  df <- data.frame()
  c <- ggplot() + theme_void()


  #ggarrange(a,a2,ggarrange(
  # c,b, ncol = 2, nrow = 1, widths = c(0.015,1)),ncol=1,nrow=3)
  ggarrange(c,
            b,
            ncol = 2,
            nrow = 1,
            widths = c(0.015, 1))
  #return(data.frame(perc = df.lines.ver[,5], row.names=rownames(df.lines.ver)))
}

generate_figure_highlight_ecrt(0.5,0.3)
d


library(ggplot2)

temp.df <- readRDS(file = "./y.vals.rds")
ecdf(y) -> model
model(y) -> y2
b <- plot( y2, y, pch=0, type="l", col="black", yaxt="n",  xlab  = "frequency",ylab="G3/G4 Score", axes=F, log = "")
rect(xleft=0, xright=0.68, ybottom= 0, ytop=1, col="green", border = NA )
rect(xleft=0.68, xright=1, ybottom= 0, ytop=1, col="yellow", border = NA)
lines( y2, y, lwd = 1.5 )
axis(side=2, at=c(-0.1,0,0.2,0.4,0.6,0.8,1,1.1))
axis(side=1, at=c(-0.1,0,0.2,0.4,0.6,0.8,1,1.1))
abline(h=c(0,0.5,1), lty = 2)
abline(v=0.68, lty = 1)
abline(v=0.9, lty=10)


df <- as.data.frame(cbind(y2,y))

library(plotly)


b <- 
   ggplot()+
  geom_rect(xmin=0, xmax=0.68, ymin= 0, ymax=1, col="green", fill = "green") + 
  geom_rect(xmin=0.68, xmax=1, ymin= 0, ymax=1, col="Yellow", fill = "Yellow")+
  geom_line(data = df, aes(x=y2, y=y)) 

fig <- ggplotly(b)


b <- qplot( y2, y, pch=0, type="l", col="black", yaxt="n",  xlab  = "frequency",ylab="G3/G4 Score", axes=F)
            

df.lines.hor <-
  foreach(i = 1:length(0.5),
          .combine = rbind) %do% {
            data.frame(
              x = 0,
              xend = max(which(
                y2 < 0.5[i]
              )),
              y = 0.5[i],
              yend = 0.5[i]
            )
          }
df.lines.hor$labels <- names(0.5)

df.lines.ver <-
  foreach(i = 1:length(0.5),
          .combine = rbind) %do% {
            data.frame(
              x = max(which(
                y2 < 0.5[i]
              )),
              xend = max(which(
                y2 < 0.5[i]
              )),
              y = 0.5[i],
              yend = min(y2)
            )
          }
df.lines.ver$perc <-
  paste0(round(df.lines.ver$xend / length(y2) * 100), "th")

df.lines.ver$colour <- factor(ifelse(1:nrow(df.lines.ver)==0.5,"highlight","no.highlight"), levels = c("highlight","no.highlight"))
df.lines.hor$colour <- factor(ifelse(1:nrow(df.lines.hor)==0.5,"highlight","no.highlight"), levels = c("highlight","no.highlight"))

b1 <- b +
  geom_segment(
    aes(
      x = x,
      y = y,
      xend = xend,
      yend = yend,
      colour = as.character(colour)
    ),
    #colour = "red",
    linetype = "dashed",
    data = df.lines.hor
  ) +
  geom_segment(
    aes(
      x = x,
      y = y,
      xend = xend,
      yend = yend,
      colour = colour
    ),
    #colour = "red",
    linetype = "dashed",
    data = df.lines.ver
  ) +
  #geom_text_repel(aes(x = x+10, y = y+0.1, label = labels), direction = "y", data = df.lines.hor)
  geom_text(aes(
    x = x + 10,
    y = y + 0.1,
    label = labels,
    colour = colour
  ), data = df.lines.hor)
#scale_color_discrete("red", "lightgrey")
c <- ggplot() + theme_void()

b
#ggarrange(a,a2,ggarrange(
# c,b, ncol = 2, nrow = 1, widths = c(0.015,1)),ncol=1,nrow=3)
ggarrange(c,
          b,
          ncol = 2,
          nrow = 1,
          widths = c(0.015, 1))







library(plotly)

library(dplyr)

df <- data.frame(name = c("Nixon", "Ford", "Carter", "Reagan", "Bush", "Clinton", "Bush", "Obama"),
                 start = as.Date(c("1969-01-20", "1974-08-09", "1977-01-20", "1981-01-20",
                                   "1989-01-20", "1993-01-20", "2001-01-20", "2009-01-20")),
                 end = as.Date(c("1974-08-09", "1977-01-20", "1981-01-20", "1989-01-20", 
                                 "1993-01-20", "2001-01-20", "2009-01-20", "2017-01-20")),
                 party = c("R", "R", "D", "R", "R", "D", "R", "D"),
                 stringsAsFactors = FALSE) %>%
  mutate(median_x = start + floor((end-start)/2))

p <- ggplot(economics, aes(x=date,y=unemploy)) +
  geom_rect(data=df, aes(NULL,NULL,xmin=start,xmax=end,fill=party),
            ymin=0,ymax=16000, colour="white", size=0.5, alpha=0.2) +
  scale_fill_manual(values=c("R" = "red", "D" = "blue")) +
  geom_line() +
  geom_text(data=df,aes(x=median_x,y=3000,label=name), size=3) +
  labs(title = "Unemmployment numbers since 1967",
       y = "No. unemployed (x 1000)")
fig <- ggplotly(p)

fig


df1 <- data.frame(green = c(0,0.68),
                  yellow = c(0.68,1),
                  colour = c("G", "Y"))

g <- ggplot(df, aes(x=y2, y=y)) +
  geom_rect(data = df1, aes(NULL,NULL, xmin=green, xmax=yellow, fill=colour),
            ymin=0,ymax=1, colour="white", size=0.5, alpha=0.2) +
  scale_fill_manual(values=c("G"= "green", "Y" = "yellow")) +
  geom_hline(yintercept=0, linetype="dashed") +
  geom_hline(yintercept=0.5, linetype="dashed") +
  geom_hline(yintercept=1, linetype="dashed") +
  geom_vline(xintercept=0.68, linetype="dashed") +
  geom_line() 

figure <- ggplotly(g)

figure

generate_figure_highlight_g3g4 <- function(new.sample.meta.score, indexRow) {
  if(is.null(indexRow)){indexRow=1}
  temp.df <- readRDS(file = "./y.vals.rds")
  ecdf(y) -> model
  model(y) -> y2
  ggplot(df, aes(x=y2, y=y)) +
    geom_rect(data = df1, aes(NULL,NULL, xmin=green, xmax=yellow, fill=colour),
              ymin=0,ymax=1, colour="white", size=0.5, alpha=0.2) +
    scale_fill_manual(values=c("G"= "green", "Y" = "yellow")) +
    geom_hline(yintercept=0, linetype="dashed") +
    geom_hline(yintercept=0.5, linetype="dashed") +
    geom_hline(yintercept=1, linetype="dashed") +
    geom_vline(xintercept=0.68, linetype="dashed") +
    geom_line() 
 
  
  
  
  df.lines.hor <-
    foreach(i = 1:length(new.sample.meta.score),
            .combine = rbind) %do% {
              data.frame(
                x = 0,
                xend = max(which(
                  y2 < new.sample.meta.score[i]
                )),
                y = new.sample.meta.score[i],
                yend = new.sample.meta.score[i]
              )
            }
  df.lines.hor$labels <- names(new.sample.meta.score)
  
  df.lines.ver <-
    foreach(i = 1:length(new.sample.meta.score),
            .combine = rbind) %do% {
              data.frame(
                x = max(which(
                  y2 < new.sample.meta.score[i]
                )),
                xend = max(which(
                  y2 < new.sample.meta.score[i]
                )),
                y = new.sample.meta.score[i],
                yend = min(y2)
              )
            }
  df.lines.ver$perc <-
    paste0(round(df.lines.ver$xend / length(y2) * 100), "th")
  
  df.lines.ver$colour <- factor(ifelse(1:nrow(df.lines.ver)==indexRow,"highlight","no.highlight"), levels = c("highlight","no.highlight"))
  df.lines.hor$colour <- factor(ifelse(1:nrow(df.lines.hor)==indexRow,"highlight","no.highlight"), levels = c("highlight","no.highlight"))
  
  b <- b +
    geom_segment(
      aes(
        x = x,
        y = y,
        xend = xend,
        yend = yend,
        colour = as.character(colour)
      ),
      #colour = "red",
      linetype = "dashed",
      data = df.lines.hor
    ) +
    geom_segment(
      aes(
        x = x,
        y = y,
        xend = xend,
        yend = yend,
        colour = colour
      ),
      #colour = "red",
      linetype = "dashed",
      data = df.lines.ver
    ) +
    #geom_text_repel(aes(x = x+10, y = y+0.1, label = labels), direction = "y", data = df.lines.hor)
    geom_text(aes(
      x = x + 10,
      y = y + 0.1,
      label = labels,
      colour = colour
    ), data = df.lines.hor)
  #scale_color_discrete("red", "lightgrey")
  
  df <- data.frame()
  c <- ggplot() + theme_void()
  
  
  #ggarrange(a,a2,ggarrange(
  # c,b, ncol = 2, nrow = 1, widths = c(0.015,1)),ncol=1,nrow=3)
  ggarrange(c,
            b,
            ncol = 2,
            nrow = 1,
            widths = c(0.015, 1))
  #return(data.frame(perc = df.lines.ver[,5], row.names=rownames(df.lines.ver)))
}

generate_figure_highlight_g3g4(0.5,0.3)
d

new.sample.meta.score = 1
indexRow = 1
generate_figure_highlight_g3g4 <- function(new.sample.meta.score, indexRow){
if(is.null(indexRow)){indexRow=1}
  temp.df <- readRDS(file = "./y.vals.rds")
  ecdf(y) -> model
  model(y) -> y2
  b <- ggplot(df, aes(x=y2, y=y)) +
    geom_rect(data = df1, aes(NULL,NULL, xmin=green, xmax=yellow, fill=colour),
              ymin=0,ymax=1, colour="white", size=0.5, alpha=0.2) +
    scale_fill_manual(values=c("G"= "green", "Y" = "yellow")) +
    geom_hline(yintercept=0, linetype="dashed") +
    geom_hline(yintercept=0.5, linetype="dashed") +
    geom_hline(yintercept=1, linetype="dashed") +
    geom_vline(xintercept=0.68, linetype="dashed") +
    geom_line() 
  
  
  
  
  df.lines.hor <-
    foreach(i = 1:length(new.sample.meta.score),
            .combine = rbind) %do% {
              data.frame(
                x = 0,
                # xend = max(which(
                #   y2 < new.sample.meta.score[i]
                # )),
                xend = 1,
                y = new.sample.meta.score[i],
                yend = new.sample.meta.score[i]
              )
            }
  df.lines.hor$labels <- names(new.sample.meta.score)
  
  df.lines.ver <-
    foreach(i = 1:length(new.sample.meta.score),
            .combine = rbind) %do% {
              data.frame(
                x = 1,
                # x = max(which(
                #   y2 < new.sample.meta.score[i]
                # )),
                # xend = max(which(
                #   y2 < new.sample.meta.score[i]
                # )),
                xend = 1,
                y = new.sample.meta.score[i],
                yend = min(y2)
              )
            }
  df.lines.ver$perc <-
    paste0(round(df.lines.ver$xend / length(y2) * 100), "th")
  
  df.lines.ver$colour <- factor(ifelse(1:nrow(df.lines.ver)==indexRow,"highlight","no.highlight"), levels = c("highlight","no.highlight"))
  df.lines.hor$colour <- factor(ifelse(1:nrow(df.lines.hor)==indexRow,"highlight","no.highlight"), levels = c("highlight","no.highlight"))
  
  
  b <- b +
    geom_segment(
      aes(
        x = x,
        y = y,
        xend = xend,
        yend = yend,
        colour = "green"
      ),
      #colour = "red",
      linetype = "dashed",
      data = df.lines.hor
    ) +
    geom_segment(
      aes(
        x = x,
        y = y,
        xend = xend,
        yend = yend,
        colour = "red"
      ),
      #colour = "red",
      linetype = "dashed",
      data = df.lines.ver
    ) 
    #geom_text_repel(aes(x = x+10, y = y+0.1, label = labels), direction = "y", data = df.lines.hor)
    geom_text(aes(
      x = x + 10,
      y = y + 0.1,
      label = labels,
      colour = "green"
    ), data = df.lines.hor)
  #scale_color_discrete("red", "lightgrey")
  
  df <- data.frame()
  c <- ggplot() + theme_void()
  
  
  #ggarrange(a,a2,ggarrange(
  # c,b, ncol = 2, nrow = 1, widths = c(0.015,1)),ncol=1,nrow=3)
 d<- ggarrange(c,
            b,
            ncol = 2,
            nrow = 1,
            widths = c(0.015, 1))
 
 d
  #return(data.frame(perc = df.lines.ver[,5], row.names=rownames(df.lines.ver)))
}
  figuretest <- ggplotly(d)
figuretest
generate_figure_highlight_g3g4(0.5,0.3)

g

################################################################### survival

time.comb <- readRDS(file = "./time.comb.rds")
age.comb <- readRDS(file = "./age.comb.rds")
status.comb <- readRDS(file = "./status.comb.rds")
comb.cont <-readRDS(file = "./comb.cont.rds")

library(survival)

coxph(Surv(time.comb, status.comb) ~ comb.cont) -> train.fit
summary(survfit(train.fit, data.frame(g3g4.values=comb.cont)), time = 5) -> x

df2 <- data.frame(pred = comb.cont,
                  surv = as.numeric(x$surv),
                  up = as.numeric(x$upper),
                  lo = as.numeric(x$lower)
)

ggplot(df2, aes(x=pred, y=surv)) +
  geom_line() +
  geom_point(alpha = 1/20, size = 2) +
  geom_line(aes(x=pred, y=lo),linetype="dotted") +
  geom_line(aes(x=pred, y=up),linetype="dotted") +
  theme_classic() + xlab("Prediction Metagene") + ylab("Survival") +
  labs(title = "New plot title", subtitle = "A subtitle") +
  ylim(0,1)


#############################

coxph(Surv(time.comb, c(status.comb)) ~ comb.cont + age.comb) -> train.fit.age

summary(survfit(train.fit.age, data.frame(comb.cont=comb.cont,
                                          age.comb = age.comb)), time = 5) -> x

summary(survfit(train.fit.age, data.frame(comb.cont=c(seq(0,1,0.1),seq(0,1,0.1)),
                                          age.comb=c(rep("TRUE",11),rep("FALSE",11)))), time = 5) -> y

df2 <- data.frame(pred = comb.cont[row.names(x$table)],
                  surv = as.numeric(x$surv),
                  up = as.numeric(x$upper),
                  lo = as.numeric(x$lower),
                  age = age.comb[row.names(x$table)]
)


df2.y <- data.frame(pred = c(seq(0,1,0.1),seq(0,1,0.1)),
                    surv = as.numeric(y$surv),
                    up = as.numeric(y$upper),
                    lo = as.numeric(y$lower),
                    age = c(rep("TRUE",11),rep("FALSE",11))
)

df2$point <- rep("yes",nrow(df2))
df2.y$point <- rep("no",nrow(df2.y))

ggplot(df2, aes(x=pred, y=surv, group=age, color = age)) +
  #geom_line() +
  geom_point(alpha = 1/10) +
  geom_line(data = df2.y, aes(x=pred, y=surv, group=age, color = age)) +
  geom_line(data = df2.y, aes(x=pred, y=lo, group=age),linetype="dotted") +
  geom_line(data = df2.y, aes(x=pred, y=up, group=age),linetype="dotted") +
  theme_classic() + xlab("Prediction Metagene") + ylab("Survival") +
  scale_color_manual(values=c('red','dodgerblue')) +
  labs(title = "New plot title", subtitle = "A subtitle") +
  ylim(0,1)

################



coxph(Surv(time.comb, c(status.comb)) ~ comb.cont + age.comb) -> train.fit.age

summary(survfit(train.fit.age, data.frame(new.sample.meta.score=comb.cont,
                                          age.comb = age.comb)), time = 5) -> x

summary(survfit(train.fit.age, data.frame(comb.cont=c(seq(0,1,0.1),seq(0,1,0.1)),
                                          age.comb=c(rep("TRUE",11),rep("FALSE",11)))), time = 5) -> y
