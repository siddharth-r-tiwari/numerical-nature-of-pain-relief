library(ggplot2)
library(patchwork)
library(ggdist)
library(wesanderson)
library(pbmcapply)
library(ggh4x)

gen_additive <- function(mu, tau, sigma, desired_n, n_points, delta) {
  n_sub <- desired_n*3*n_points
  
  alpha <- rnorm(n_sub, mu, tau)
  X <- rnorm(n_sub * n_points, rep(alpha, each=n_points), sigma) 
  Y <- rnorm(n_sub * n_points, rep(alpha, each=n_points) + delta, sigma) 
  
  df <- data.frame(sub = rep(1:n_sub, each = n_points),
                   pre = X,
                   post = Y)
  
  to_remove <- unique(df[which(df[,2] < 0 | df[,3] < 0 | df[,2] > 100 | df[,3] > 100),"sub"])
  df <- df[-to_remove,]
  
  aggregate(cbind(pre,post) ~ sub, df[1:desired_n,], mean)
}

gen_mult <- function(mu, tau, sigma, desired_n, n_points, delta) {
  n_sub <- desired_n*3*n_points
  
  alpha <- rlnorm(n_sub, 
                  log(mu^2/sqrt(mu^2 + tau^2)), 
                  sqrt(log(1+(tau^2/mu^2))))
  
  X <- rlnorm(n_sub * n_points, 
              rep(log(alpha), each=n_points), 
              sigma/mu)
  
  
  alpha2 <- alpha * (1 + delta/mu)
  Y <- rlnorm(n_sub * n_points, 
              rep(log(alpha2), each=n_points), 
              sigma/mu)
  
  df <- data.frame(sub = rep(1:n_sub, each = n_points),
                   pre = X,
                   post = Y)
  
  
  to_remove <- unique(df[which(df[,2] < 0 | df[,3] < 0 | df[,2] > 100 | df[,3] > 100),"sub"])
  df <- df[-to_remove,]
  
  aggregate(cbind(pre,post) ~ sub, df[1:desired_n,], function(x) exp(mean(log(x))))
}

##specify params for both models
mu <- 70
tau <- 15
sigma <- 7.5
n_sub <- 500
n_points <- 1
delta <- -30

#####
# FIGURE 1: Hierarchical explanation
set.seed(1000)
alphas <- rnorm(4,0,tau) + mu
ggplot() +
  stat_dist_slab(aes(y = 0,
                     dist = "norm",
                     arg1 = mu,
                     arg2 = tau),
                 color = NA,
                 fill = wes_palette("BottleRocket2", 1)) +
  geom_segment(aes(x = mu,
                   xend = mu,
                   y = 0,
                   yend = 0.90 )) +
  stat_dist_slab(aes(y = 0,
                     dist = "norm",
                     arg1 = alphas,
                     arg2 = sigma),
                 height = .35,
                 color = NA,
                 alpha = 0.5,
                 fill = wes_palette("BottleRocket2", 2)[2]) +
  geom_segment(aes(x = alphas,
                   xend = alphas,
                   y = 0,
                   yend = .35 * 0.90 ),
               color = "white") +
  geom_segment(aes(x = 78,
                   xend = 63,
                   y = .75,
                   yend = 0.6),
               size = 1,
               arrow = arrow(length = unit(12, "pt")),
               lineend = "round",
               linejoin = "bevel") +
  geom_text(aes(x = 80,
                y = 0.8,
                label = "Distribution of patients' mean pain scores")) +
  geom_segment(aes(x = 15,
                   xend = 20,
                   y = .5,
                   yend = 0.2),
               size = 1,
               arrow = arrow(length = unit(12, "pt")),
               lineend = "round",
               linejoin = "bevel") +
  geom_text(aes(x = 15,
                y = .55,
                label = "Single patient distribution")) +
  
  coord_cartesian(xlim=c(0,100),
                  expand = F) +
  xlab("Pain rating (VAS)") +
  theme_minimal() +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        panel.grid = element_blank(),
        plot.margin = margin(2,12,2,12))

ggsave("figures/figure1.pdf", width=8, height=3, device = cairo_pdf)


#####
# Figure 2. Properties of additive and multiplicative data

set.seed(1)
df.add <- gen_additive(mu, tau, sigma, n_sub, n_points, delta)

set.seed(1)
df.mult <- gen_mult(mu, tau, sigma, n_sub, n_points, delta)

dfs <- rbind(cbind(type="Additive", df.add),
             cbind(type="Multiplicative", df.mult))

B1 <-
  ggplot() +
  geom_point(data = dfs, aes(x = pre, y = post),
             color = wes_palette("FantasticFox1", 5)[3],
             alpha = 0.25,
             size = 0.5) +
  facet_wrap(. ~ type, nrow = 1) +
  ylim(0,100) +
  xlim(0,100) +
  xlab("Pre-intervention pain (VAS)") +
  ylab("Post-intervention pain (VAS)")

C1 <-
  ggplot() +
    geom_ribbon(data = data.frame(x = 0:100,
                                  ymax = -1 * 0:100,
                                  ymin = -100),
              aes(ymin = ymin,
                  ymax = ymax,
                  x = x),
              fill = "grey") +
    geom_ribbon(data = data.frame(x = 0:100,
                                  ymin = 100 - 0:100,
                                  ymax = 100),
                aes(ymin = ymin,
                    ymax = ymax,
                    x = x),
                fill = "grey") +
    geom_point(data = dfs, 
               aes(x = pre, y = post-pre),
               color = wes_palette("FantasticFox1", 5)[3],
               alpha = 0.25,
               size = 0.5) +
    facet_wrap(. ~ type, nrow = 1) +
    ylim(-100,100) +
    xlim(0,100) +
    xlab("Pre-intervention pain (VAS)") +
    ylab("Change in pain (ΔVAS)")

B1 + C1 + plot_layout(nrow = 2) + plot_annotation(tag_levels = "A") & theme_minimal() & coord_cartesian(expand = F) & 
  theme(panel.background = element_rect(fill = "#f9f9f9",
                                        color = "#ffffff"),
        panel.spacing.x = unit(20, "pt"),
        aspect.ratio = 1)

ggsave("figures/figure2.pdf", width = 4, height = 5, device = cairo_pdf)

#####
# Figure 2. Slopes, # of points, and their relationship with ICC (RTM)

# Different ICC's
set.seed(1)
iccs <- seq(0.05,0.95,by=0.05)
sigmas <- sqrt(-(iccs-1)*tau^2)/sqrt(iccs)

df.mult <- df.add <- c()
for(i in 1:length(iccs)) {
  sigma <- sigmas[i]
  icc <- iccs[i]
  print(icc)
  tmp <- gen_additive(mu, tau, sigma, n_sub*1000, n_points, delta)
  df.add <- rbind(df.add,
                  data.frame(icc = icc, slope = coef(lm(tmp$post - tmp$pre ~ tmp$pre))[2]))
  
  tmp <- gen_mult(mu, tau, sigma, n_sub*1000, n_points, delta)
  df.mult <- rbind(df.mult,
                   data.frame(icc = icc, slope = coef(lm(tmp$post - tmp$pre ~ tmp$pre))[2]))
}

dfs <- rbind(cbind(type = "Additive", df.add),
             cbind(type = "Multiplicative", df.mult))

write.csv(dfs, "data/icc_slopes.csv", row.names = F, quote = F)

dfs <- read.csv("data/icc_slopes.csv")

ggplot(dfs,
       aes(x = icc,
           y = slope)) +
  geom_point() +
  ylim(c(-1,0)) +
  xlim(c(0,1)) +
  xlab("Intraclass correlation coefficient") +
  ylab("Slope of change vs. pre") +
  facet_wrap(. ~ type) +
  theme_minimal() +
  theme(panel.background = element_rect(fill = "#f9f9f9",
                                        color = "#ffffff"))

ggsave("figures/figure3.pdf", width = 8, height = 2.5, device = cairo_pdf)



#####
# ancovas

set.seed(1)
df.add <- gen_additive(mu, tau, sigma, n_sub, n_points, delta)

set.seed(1)
df.mult <- gen_mult(mu, tau, sigma, n_sub, n_points, delta)

plot(df.add)
plot(df.mult)

df.ancova <- rbind(
  data.frame(type = "Additive",
             model = "Raw",
             resids = resid(lm(df.add$post ~ df.add$pre)),
             fitted = fitted(lm(df.add$post ~ df.add$pre))),
  data.frame(type = "Additive",
             model = "Log",
             resids = resid(lm(log(df.add$post) ~ log(df.add$pre))),
             fitted = exp( fitted(lm(log(df.add$post) ~ log(df.add$pre)))) ),
  data.frame(type = "Multiplicative",
             model = "Raw",
             resids = resid(lm(df.mult$post ~ df.mult$pre)),
             fitted = fitted(lm(df.mult$post ~ df.mult$pre))),
  data.frame(type = "Multiplicative",
             model = "Log",
             resids = resid(lm(log(df.mult$post) ~ log(df.mult$pre))),
             fitted = exp( fitted(lm(log(df.mult$post) ~ log(df.mult$pre)))) )
)

df.ancova.trimmed <- df.ancova[-which(abs(df.ancova$resids) > 1 & df.ancova$type == "Additive" & df.ancova$model == "Log"),]
df.ancova$model <- factor(df.ancova$model, c("Raw","Log"))
df.ancova.trimmed$model <- factor(df.ancova.trimmed$model, c("Raw","Log"))

blanks <- rbind(data.frame(type = "Additive",
                           model = "Log",
                           x = c(40,50),
                           y = c(0,1)),
                data.frame(type = "Multiplicative",
                           model = "Raw",
                           x = c(40,50),
                           y = c(0,30)))

blanks$model <- factor(blanks$model, levels = c("Raw","Log"))

ggplot() +
  geom_point(data = df.ancova.trimmed,
             aes(x = fitted,
                 y = abs(resids)),
             alpha = 0.25,
             color = wes_palette("FantasticFox1", 5)[3]) +
  geom_smooth(data = df.ancova,
              aes( x = fitted,
                   y = abs(resids)),
              method = "lm",
              color = wes_palette("FantasticFox1", 5)[4],
              fill = wes_palette("FantasticFox1", 5)[4]) +
  geom_blank(data = blanks, aes(x,y)) +
  facet_grid(model ~ type, scales = "free") +
  coord_cartesian(ylim = c(0,NA), expand = F) +
  xlab("Fitted post-intervention pain") +
  ylab("Absolute value of residuals") +
  scale_x_continuous(breaks = seq(0,100,by=10)) +
  theme_minimal() +
  theme(panel.background = element_rect(fill = "#f9f9f9",
                                        color = "#ffffff"),
        panel.spacing = unit(20, "pt"))

ggsave("figures/figure4.pdf", width = 8, height = 4, device = cairo_pdf)

#####
# placebo 2
ratings_pl2 <- read.csv("data/placeboII_ratings.csv", strip.white = TRUE, stringsAsFactors = FALSE)
groups_pl2 <- readxl::read_excel("data/Data_paulo.xlsx")

get_pl2 <- function(npoints) {
  pain_pl2 <- sapply(unique(ratings_pl2$subject), function(sub) {
    tmp = subset(ratings_pl2, subject == sub)
    sub_group = groups_pl2[groups_pl2$SID == sub,]$group
    tmp <- na.omit(tmp)
    
    if(sum(as.numeric(tmp$phase == 1)) > 0 & length(sub_group) == 1){
      data.frame(pre  = mean(head(tmp[tmp$phase == 1, "pain"], npoints))*10,
                 post = mean(tail(tmp$pain, npoints))*10,
                 group = sub_group) 
    }
  })
  pain_pl2 <- do.call(rbind, pain_pl2)
  pain_pl2 <- na.omit(pain_pl2)
  pain_pl2 <- as.data.frame(pain_pl2)
  colnames(pain_pl2) <- c("pre", "post", "group")
  
  pain_pl2 <- cbind(pain_pl2,
                    logpre = log(pain_pl2$pre + 1),
                    logpost = log(pain_pl2$post + 1))
  
  pain_pl2$group <- factor(pain_pl2$group, levels = c("notx","placebo","naproxen"))
  
  pain_pl2
}


#####
# placebo 1
ratings_pl1 <- read.csv("data/placeboI_ratings.csv", strip.white = TRUE, stringsAsFactors = FALSE)
groups_pl1 <- readxl::read_excel("data/placeboI_groups.xlsx")[,c(1,78)]
groups_pl1 <- as.data.frame(groups_pl1)

get_pl1 <- function(npoints) {
  
  pain_pl1 <- (lapply(unique(ratings_pl1$subject), function(sub) {
    tmp = subset(ratings_pl1, subject == sub)
    if(sub %in% groups_pl1[,1]){
      sub_group = groups_pl1[[which(groups_pl1[,1] == sub), 2]]
      tmp <- na.omit(tmp)
      
      if(sum(as.numeric(tmp$phase == 1)) > 0 & !is.null(sub_group)){
        cbind(pre = mean(head(tmp[tmp$phase == 1, "pain"], npoints))*10,
              post = mean(tail(tmp$pain, npoints))*10, group=sub_group) 
      }
    }
    
  }))
  pain_pl1 <- do.call(rbind, pain_pl1)
  pain_pl1 <- na.omit(pain_pl1)
  pain_pl1 <- as.data.frame(pain_pl1)
  colnames(pain_pl1) <- c("pre", "post", "group")
  
  pain_pl1 <- cbind(pain_pl1,
                    logpre = log(pain_pl1$pre + 1),
                    logpost = log(pain_pl1$post + 1))
  
  pain_pl1 <- pain_pl1[-c(6),]
  pain_pl1$group <- factor(pain_pl1$group, levels = c(3,1), labels = c("notx","placebo"))
  
  pain_pl1
}

#####
# sbp
ratings_sbp <- read.csv("data/sbp_ratings.csv", strip.white = TRUE, stringsAsFactors = FALSE)
ratings_sbp$pain_vas  <- as.numeric(ratings_sbp$pain_vas)
ratings_sbp <- na.omit(ratings_sbp)

rownames(ratings_sbp) <- 1:nrow(ratings_sbp)

for(i in 1:nrow(ratings_sbp)) {
  if(ratings_sbp[i,3] <= 1){
    ratings_sbp[i,3] <- ratings_sbp[i,3]*100
  }
}

sel <- subset(ratings_sbp, session=="interview")$sub
ratings_sbp <- subset(ratings_sbp, sub %in% sel)

pain_sbp <- (lapply(unique(ratings_sbp$sub), function(SId) {
  tmp = subset(ratings_sbp, sub == SId)
  tmp <- na.omit(tmp)
  
  if(nrow(tmp) > 3){
    cbind(pre = tmp$pain_vas[1],
          post = tmp$pain_vas[5]) 
  }
}))
pain_sbp <- do.call(rbind, pain_sbp)
pain_sbp <- na.omit(pain_sbp)
pain_sbp <- as.data.frame(pain_sbp)
colnames(pain_sbp) <- c("pre", "post")

pain_sbp <- cbind(pain_sbp,
                  group = as.factor(0),
                  logpre = log(pain_sbp$pre + 1),
                  logpost = log(pain_sbp$post + 1))



#####
# ldopa
ratings_ldopa <- read.csv("data/ldopa_ratings.csv", strip.white = TRUE, stringsAsFactors = FALSE)
ratings_ldopa <- na.omit(ratings_ldopa)

groups_ldopa <- read.csv("data/ldopa_groups.csv")[,c(1,4)]
groups_ldopa$group <- factor(groups_ldopa$treat, levels = 1:3, labels = c("l_dopa", "placebo", "N/A"))
groups_ldopa$id <- paste0("SAT",sprintf("%03d", groups_ldopa$id))

get_ldopa <- function(npoints) {
  
  pain_ldopa <- (lapply(unique(ratings_ldopa$id), function(sub) {
    tmp = subset(ratings_ldopa, sub == id)
    if(sub %in% groups_ldopa[,1]){
      sub_group = groups_ldopa[which(groups_ldopa[,1] == sub), 3]
      
      if(nrow(tmp) > 30 & !is.null(sub_group)){
        tmp <- na.omit(tmp)
        cbind(pre = mean(head(tmp$pain, npoints))*10,
              post = mean(tail(tmp$pain, npoints))*10,
              group = sub_group) 
      } 
    }
  }))
  
  pain_ldopa <- do.call(rbind, pain_ldopa)
  pain_ldopa <- na.omit(pain_ldopa)
  pain_ldopa <- as.data.frame(pain_ldopa)
  colnames(pain_ldopa) <- c("pre", "post", "group")
  
  pain_ldopa <- cbind(pain_ldopa,
        logpre = log(pain_ldopa$pre + 1),
        logpost = log(pain_ldopa$post + 1))
  
  pain_ldopa$group <- factor(pain_ldopa$group, levels = 1:2, labels = c("l_dopa", "placebo"))
  
  pain_ldopa
}




pain_merged <- rbind(cbind(get_pl1(5), study = "Placebo 1"),
                     cbind(get_pl2(5), study = "Placebo 2"),
                     cbind(pain_sbp, study = "Longitudinal cohort"),
                     cbind(get_ldopa(5), study = "Levodopa"))

pain_merged$study <- factor(pain_merged$study, levels = unique(pain_merged$study))

######
# SCATTERS

pain_tmp <- rbind(cbind(pain_merged, type="Delta"),
                  cbind(pain_merged, type="Post"))
pain_tmp$type <- factor(pain_tmp$type, c("Delta","Post"))

ggplot() +
  geom_jitter(data = subset(pain_tmp, type == "Post"),
              aes(x = pre,
                  y = post),
              width = 0,
              height = 0,
              color = wes_palette("FantasticFox1", 5)[3],
              alpha = 0.5) +
  geom_smooth(data = subset(pain_tmp, type == "Post"),
              aes(x = pre,
                  y = post),
              method = "lm",
              color = wes_palette("FantasticFox1", 5)[4],
              fill = wes_palette("FantasticFox1", 5)[4]) +
  
  geom_jitter(data = subset(pain_tmp, type == "Delta"),
              aes(x = pre,
                  y = post-pre),
              width = 0,
              height = 0,
              color = wes_palette("FantasticFox1", 5)[3],
              alpha = 0.5) +
  geom_smooth(data = subset(pain_tmp, type == "Delta"),
              aes(x = pre,
                  y = post-pre),
              method = "lm",
              color = wes_palette("FantasticFox1", 5)[4],
              fill = wes_palette("FantasticFox1", 5)[4]) +
  
  facet_grid(type ~ study, scales = "free_y") +
  ylab("Pain") +
  xlab("Pre-intervention pain") +
  scale_x_continuous(breaks = seq(0,100,by=20)) +
  theme_minimal() +
  theme(panel.background = element_rect(fill = "#f9f9f9",
                                        color = "#ffffff"),
        panel.spacing = unit(20, "pt"))

ggsave("figures/figure5.pdf", width = 8, height = 4, device = cairo_pdf)




######
# ANCOVAS

pain_tmp <- rbind(cbind(pain_merged, type = "Raw"),
                  cbind(pain_merged, type = "Log"))
pain_tmp$type <- factor(pain_tmp$type, c("Raw","Log"))
pain_tmp[which(pain_tmp$type == "Log"), c("pre","post")] <- pain_tmp[which(pain_tmp$type == "Log"), c("logpre","logpost")]

ggplot() +
  geom_jitter(data = subset(pain_tmp, type == "Raw"),
              aes(x = (fitted(lm(post ~ pre + group, subset(pain_tmp, type == "Raw")))),
                  y = (abs(resid(lm(post ~ pre + group, subset(pain_tmp, type == "Raw")))))),
              color = wes_palette("FantasticFox1", 5)[3],
              width = 0,
              height = 0,
              alpha = 0.5) +
  geom_smooth(data = subset(pain_tmp, type == "Raw"),
    aes(x = fitted(lm(post ~ pre + group, subset(pain_tmp, type == "Raw"))),
        y = abs(resid(lm(post ~ pre + group, subset(pain_tmp, type == "Raw"))))),
    method = "lm",
    color = wes_palette("FantasticFox1", 5)[4],
    fill = wes_palette("FantasticFox1", 5)[4]) +
  
  geom_jitter(data = subset(pain_tmp, type == "Log"),
              aes(x = exp(fitted(lm(logpost ~ logpre + group, subset(pain_tmp, type == "Log")))),
                  y = abs(resid(lm(logpost ~ logpre + group, subset(pain_tmp, type == "Log"))))),
              color = wes_palette("FantasticFox1", 5)[3],
              width = 0,
              height = 0,
              alpha = 0.5) +
  geom_smooth(data = subset(pain_tmp, type == "Log"),
              aes(x = exp(fitted(lm(logpost ~ logpre + group, subset(pain_tmp, type == "Log")))),
                  y = abs(resid(lm(logpost ~ logpre + group, subset(pain_tmp, type == "Log"))))),
              method = "lm",
              color = wes_palette("FantasticFox1", 5)[4],
              fill = wes_palette("FantasticFox1", 5)[4]) +
  
  facet_grid(type ~ study, scales = "free_y") +
  ylab("Absolute value of residuals") +
  xlab("Fitted post-intervention pain") +
  scale_x_continuous(breaks = seq(0,100,by=20)) +
  theme_minimal() +
  theme(panel.background = element_rect(fill = "#f9f9f9",
                                        color = "#ffffff"),
        panel.spacing = unit(20, "pt"))

ggsave("figures/figure7.pdf", width = 8, height = 4, device = cairo_pdf)




pain_output <- rbind(cbind(get_pl1(5), study = "Placebo 1"),
                     cbind(get_pl2(5), study = "Placebo 2"),
                     cbind(get_ldopa(5), study = "Ldopa"))

for(i in unique(pain_output$study)) {
  print(i)
  tmp <- subset(pain_output, study == i)
  print(cbind(coef(lm(post ~ pre + group, tmp)),confint(lm(post ~ pre + group, tmp))))
  print(exp(cbind(coef(lm(logpost ~ logpre + group, tmp)),confint(lm(logpost ~ logpre + group, tmp)))))
}


slopes <- c()
for(i in 1:10) {
  pain_output <- rbind(cbind(get_pl1(i), study = "Placebo 1"),
                       cbind(get_pl2(i), study = "Placebo 2"),
                       cbind(get_ldopa(i), study = "Levodopa"))
  
  slope <- lapply(unique(pain_output$study), function(j) {
    data.frame(study = j,
               points = i,
               slope = coef(lm(post - pre ~ pre, subset(pain_output, study == j)))[2])
  })
  slope <- do.call(rbind.data.frame, slope)
  
  slopes <- rbind(slopes, slope)
  
}

slopes$study <- factor(slopes$study, levels = unique(slopes$study))

ggplot(slopes,
       aes(x = points,
           y = slope)) +
  geom_point() +
  facet_grid(. ~ study) + 
  ylim(-1,0) +
  xlab("Number of points averaged") +
  scale_x_continuous(breaks = seq(1,10,by=1)) +
  ylab("Slope of change vs. pre") +
  theme_minimal() +
  theme(panel.background = element_rect(fill = "#f9f9f9",
                                        color = "#ffffff"),
        panel.spacing = unit(20, "pt"),
        panel.grid.minor = element_blank())


ggsave("figures/figure6.pdf", width = 8, height = 2.5, device = cairo_pdf)



