args = commandArgs(trailingOnly=TRUE)

# Load in required packages
suppressMessages(library("plyr"))
suppressMessages(library("MCMCpack"))
suppressMessages(library("R.utils"))
suppressMessages(library("ggplot2"))
suppressMessages(library("MBESS"))
suppressMessages(library("tidyverse"))
suppressMessages(library("ReplicationBF"))
suppressMessages(library("ggthemes"))
suppressMessages(library("BayesFactor"))
suppressMessages(library("gridExtra"))

##########################################################################
#------------------------------------------------------------------------#
##########################################################################


                  ###  BFrep Simulation 1 ### 

temp <- paste0(args[1], "\\BFrep_Sim1.rds")
results <- readRDS(file=temp)

# Prepare plot
results$n.orig = factor(results$n.orig,
                        labels = paste0("n[plain(orig)] == ",
                                        unique(results$n.orig)))
results$n.rep = factor(results$n.rep,
                       labels = paste0("n[plain(rep)] == ",
                                       unique(results$n.rep)))
results$f2.orig.val = as.numeric(results$f2.orig)
results$f2.orig = factor(results$f2.orig,
                         labels = paste0("f^2 == ", unique(results$f2.orig)))

# Plot
fig1 <- paste0(args[1], "\\Group_9_Figure_1.pdf")
pdf(file=fig1)

ggplot(results, aes(x = f2.rep, color = as.factor(k)  , y = BFrep)) +
  geom_hline(aes(yintercept = 1)) +
  geom_point() +
  geom_line() +
  scale_y_log10(name = "Replication Bayes factor (log_10-scaled)") +
  scale_x_continuous(name ="Replication Effect Size") +
  theme(legend.position = "bottom", legend.direction = "horizontal") +
  facet_grid(f2.orig ~ n.orig + n.rep, scales = "free_y", labeller = label_parsed) +
  labs(color = "# of Groups")

dev.off()
                ###  BFrep Simulation 1 (END) ### 





##########################################################################
#------------------------------------------------------------------------#
##########################################################################



                     ###  BFrep Simulation 3  ###


# Load sim 3 results
temp2 <- paste0(args[1], "\\BFrep_Sim3.rds")
data <- readRDS(file=temp2)#temp

# Generate plot
fig2 <- paste0(args[1], "\\Group_9_Figure_2.pdf")
pdf(file=fig2)

ggplot(data=data, 
       aes(x=BFrep.t.ISe, y=BFrep.F.ISe, color=as.factor(d.rep), shape = as.factor(d.orig))) +
  geom_abline(slope=1, intercept=0) +
  geom_point() +
  scale_x_log10() +
  scale_y_log10() +
  scale_shape_discrete(name =expression('d'["orig"])) +
  scale_color_discrete(name =expression('d'["rep"])) +
  facet_grid(n.rep ~ n.orig, labeller = label_both, drop = T) +
  labs(x = "Bayes Replication factor for (t-Test)", y = "Bayes Replication factor for (F-Test)")

dev.off()
# Correlation
#cor(data$BFrep.t.ISe, data$BFrep.F.ISe)


               ###  BFrep Simulation 3 (END)  ###





##########################################################################
#------------------------------------------------------------------------#
##########################################################################


                ###  BFrep Example 3 ###


# Load in data from example 3
temp3 <- paste0(args[1], "\\BFrep_Ex3.rds")
Ex3results <- readRDS(file=temp3)#temp

# Plot data 
fig3 <- paste0(args[1], "\\Group_9_Figure_3.pdf")
pdf(file=fig3)

suppressMessages(ggplot(data = Ex3results, aes(x = x)) +
  geom_line(aes(y = mean, group = study)) +
  geom_point(aes(y = mean, group = study)) +
  facet_grid(. ~ study) +
  scale_y_continuous(name = "Mean") +
  scale_x_discrete(name = "Group")) 

dev.off()

# Load in table data from example 3
temp4 <- paste0(args[1], "\\BFrep_Ex3data.rds")
Ex3data <- readRDS(file=temp4)#temp

# Plot table
fig4 <- paste0(args[1], "\\Group_9_Table_1.pdf")
pdf(file=fig4)

p<-tableGrob(Ex3data)
grid.arrange(p)

dev.off()

               ###  BFrep Example 3 (END) ###



##########################################################################
#------------------------------------------------------------------------#
##########################################################################

