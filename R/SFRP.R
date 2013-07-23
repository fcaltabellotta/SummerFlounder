# SFRP.R
# Steve Martell
# July 20, 2013

# GROWTH PARS
# Current assessment

# 38,173

# 12,14
# male female both
# linf
# 73.5, 80.9, 87.2
# k
# 0.14, 0.18, 0.14


# |----------------------------------------------------------------------------------|
# | LIBRARIES
# |----------------------------------------------------------------------------------|
# |
setwd("/Users/stevenmartell1/Documents/CONSULTING/SummerFlounder/R/")
require(ggplot2)
require(reshape2)
source("Selex.R")

# |----------------------------------------------------------------------------------|
# | DATA & VARIABLES
# |----------------------------------------------------------------------------------|
# |
A    <- 10
G    <- 31
S    <- 2
dl   <- 2.5
age  <- 0:A
dim  <- c(length(age),G,S)
pg   <- dnorm(seq(-1.96,1.96,length=G),0,1)
pg   <- pg/sum(pg)

Theta<- list(
RA   = "Two Sex",
age  = age ,
pg   = pg,
dim  = dim,
S    = S,
linf = c(78.48,65.25),   
vbk  = c(0.22,0.23),     
to   = c(-1.12,-1.50),
cvla = rep(0.1,S),
a    = rep(1e-6,S),
b    = rep(3,S),
a50  = rep(1.5,S),
g50  = rep(0.25,S),
m    = c(0.2,0.3),
bo   = 100,
h    = 0.80
)

S <- 2
dim  <- c(length(age),G,S)
T1<- list(
RA   = "Single Sex",
age  = age ,
pg   = pg,
S    = S,
dim  = dim,
linf = rep(mean(c(78.48,65.25)),2),   
vbk  = rep(mean(c(0.22,0.23)),2),     
to   = rep(mean(c(-1.12,-1.50)),2),
cvla = rep(0.1,S),
a    = rep(1e-6,S),
b    = rep(3,S),
a50  = rep(1.5,S),
g50  = rep(0.25,S),
m    = c(0.25,0.25),
bo   = 100,
h    = 0.80
)
# |----------------------------------------------------------------------------------|
# | MANAGEMENT VARIABLES
# |----------------------------------------------------------------------------------|
# |
fe   <- seq(0,1.00,by=0.01)
dm   <- 0.10
slim <- 33
ulim <- 1000
cvlm <- 0.1

# |----------------------------------------------------------------------------------|
# | Calculate survivorship, growth and fecundity (unfished)
# |----------------------------------------------------------------------------------|
# |
.calcLifeTable <- function(Theta)
{
	with(Theta, {
		lx     <- array(1, dim)
		la	   <- array(0, dim)
		sd_la  <- array(0, dim)
		wa	   <- array(0, dim)
		fa	   <- array(0, dim)
		M      <- array(0, dim)

		for( i in 1:S )
		{
			# maturity at age
			ma      <- plogis(age, a50[i], g50[i])

			# length-at-age
			mu      <- linf[i]*(1.-exp(-vbk[i]*(age-to[i])))
			sigma   <- cvla[i]*mu
			dev     <- seq(-1.96,1.96,length=G)
			if(G==1) dev <- 0
			la[,,i]    <- sapply(dev,fn<-function(dev){la=mu+dev*sigma})
			sd_la[,,i] <- sqrt(1/G*(cvla[i]*mu)^2)
			wa[,,i]    <- a[i]*la[,,i]^b[i]
			fa[,,i]    <- ma*wa[,,i]	

			# Sex-specific natural mortality rate
			# Should also look at size-dependent M (cm != 0 )
			cm      <- -0.1
			l_r     <- mean(la[,,i])
			delta   <- (la[,,i]/l_r)^cm / mean((la[,,i]/l_r)^cm)
			M[,,i]	<- m[i] * delta

			# Survivorship
			AA <- length(age)
			for( j in 2:AA )
			{
				lx[j,,i] <- lx[j-1,,i]*exp(-M[j-1,,i])
			}
			lx[AA,,i] <- lx[AA,,i]/(1-exp(-M[AA,,i]))
		}
		Theta$lx    = lx
		Theta$la    = la
		Theta$sd_la = sd_la
		Theta$wa    = wa
		Theta$fa    = fa 
		Theta$M     = M  
		return( Theta )
		})
}


# |----------------------------------------------------------------------------------|
# | Size-based selectivities and joint capture probability
# |----------------------------------------------------------------------------------|
# |
.calcSelectivities <- function(Theta)
{
	with (Theta,{
		# Length-interval midpoints for integration
		xl  <- seq(5,100,by=dl)
		
		
		# Length-based selectivity (length-based -> age-based)
		sc	<- array(0, dim)
		sr	<- array(0, dim)
		sd	<- array(0, dim)
		va	<- array(0, dim)
		vd	<- array(0, dim)  #discard fishery
		std	<- cvlm*slim+1.e-30	

		for( i in 1:S )
		{
			# probability of capturing a fish of length l
			pl      <- plogis(xl,30,6.5)  

			# probability of capturing a fish of age a given l
			sc[,,i] <- .calcPage(la[,,i],sd_la[,,i],pl,xl)

			# retention probability
			pr      <- plogis(xl,slim,0.1) - plogis(xl,ulim,0.1)
			pr      <- pr / max(pr)

			# probabilty of retaining a fish of age a
			sr[,,i] <- .calcPage(la[,,i],sd_la[,,i],pr,xl)

			# probabilty of discarding a fish of age a
			sd[,,i]  <- 1-sr[,,i]

			# vulnerability to death by fishing
			va[,,i]  <- sc[,,i]*(sr[,,i]+sd[,,i]*dm)

			# discard fishery selecitvity
			pd       <- plogis(xl, 5,  0.1) - plogis(xl, slim, 0.1)
			vd[,,i]  <- .calcPage(la[,,i],sd_la[,,i], pd, xl)
		}
		Theta$sc <- sc	# Length-based commercial selectivity.
		Theta$sr <- sr	# Age-specific retention probability.
		Theta$sd <- sd	# Age-specific discard probability.
		Theta$va <- va	# Joint capture probability.
		Theta$vd <- vd	# Discard probability in trawl fishery.
		return(Theta)
		})
}
# |---------------------------------------------------------------------------|
# | Calculate stock recruitment relationship.                       
# |---------------------------------------------------------------------------|
# |
.calcSRR <- function(Theta)
{
	with(Theta, {
		# Unfished SPR  (phi.E)
		phi.E	<- sum(t(t(lx[,,1]*fa[,,1])*pg))
		
		# Unfished recruitment (ro)
		ro		<- bo/phi.E
		
		# Beverton-Holt model
		kap <- 4*h/(1-h)

		# Ricker Model
		# kap <- (5*h)^(5/4)
		
		Theta$phi.E <- phi.E
		Theta$ro    <- ro
		Theta$kap   <- kap
		
		return(Theta)
	})
}

# |---------------------------------------------------------------------------|
# | Age-structure equilibrium model asem                            
# |---------------------------------------------------------------------------|
# | fe is the equilibrium fishing mortality rate.
.asem <- function(fe=0, Theta, ct=0)
{
	# | Psuedocode:
	# | 1. Calculate age-specific total mortality, retention, and discard rates
	# | 2. Calculate survivorship with fe>0.
	# | 3. Calculate equilibrium recruitment (re) and biomass (be)
	# | 4. Calculate yield per recruit,  spawning biomass per recruit,  yield, discards.
	# | 5. Calculate average weight-at-age.
	with(Theta, {
		# 1. Age-specific total mortality, survival, retention, and discard rate.
		za	<- array(0, dim)
		sa	<- array(0, dim)
		qa	<- array(0, dim)
		da	<- array(0, dim)
		
		bycatch <- ct
		bapprox <- bo * 0.15/(0.15+fe)
		fd      <- bycatch/bapprox
		#if(ct>0)
		#cat("fe = ", fe, " fd = ", fd, "\n")
		
		for(i in 1:S)
		{
			za[,,i]  <- M[,,i] + fe*va[,,i] + fd*vd[,,i]
			sa[,,i]  <- exp(-za[,,i])
			qa[,,i]  <- (sc[,,i]*sr[,,i]) * (1-sa[,,i])/za[,,i]
			da[,,i]  <- (sc[,,i]*sd[,,i]) * (1-sa[,,i])/za[,,i]
		}
		
		# 2. Survivorship under fished conditions lz(A, G, S)
		lz	<- array(1, dim)
		AA <- length(age)
		for(i in 1:S)
		{
			for(j in 2:AA)
			{
				lz[j,,i] <- lz[j-1,,i]*exp(-za[j-1,,i])
			}
			lz[AA,,i] <- lz[AA,,i]/(1-exp(-za[AA,,i]))
		}
		
		# 3. Calculate equilibrium recruitment and biomass
		phi.E	<- sum(t(t(lx[,,1]*fa[,,1])*pg))
		phi.e	<- sum(t(t(lz[,,1]*fa[,,1])*pg))
		
		# Beverton-Holt model
		t1      <- phi.E/phi.e
		t2      <- (kap-t1)
		re      <- max(0, ro*t2/(kap-1))
		
		# Ricker model
		#t1		<- log(phi.E/(kap*phi.e))
		#t2		<- (log(kap)*phi.e)
		#re		<- max(0, -(t1*ro*phi.E)/t2)
		
		be		<- re * phi.e
		
		# 4. Calculate yield per recruit,  spawning biomass per recruit,  yield, discards.
		ye		<- 0
		ye.sex  <- rep(0,S)
		de		<- 0
		ypr		<- 0
		dpr     <- 0
		bpr     <- 0
		spr		<- phi.e/phi.E
		for(i in 1:S)
		{
			ye  <- ye + sum( re * fe * t(lz[,,i]*wa[,,i]*qa[,,i])*pg )
			ye.sex[i] <- sum( re * fe * t(lz[,,i]*wa[,,i]*qa[,,i])*pg )
			de	<- de + sum( re * fe * dm * t(lz[,,i]*wa[,,i]*da[,,i])*pg )
			bpr <- bpr + sum( t(lz[,,i]*wa[,,i]*qa[,,i])*pg )
			ypr <- ypr + sum( fe * t(lz[,,i]*wa[,,i]*qa[,,i])*pg )
			dpr <- dpr + sum( fe * dm * t(lz[,,i]*wa[,,i]*da[,,i])*pg )
		}
		
		# 5. Calculate average weight-at-age
		wbar <- matrix(0, nrow=S, ncol=length(age))
		for(i in 1:S)
		{
			tmp      <- lz[,,i]/rowSums(lz[,,i])
			wbar[i,] <- rowSums(wa[,,i]*tmp)
		}
		
		Theta$lz   <- lz
		Theta$re   <- re
		Theta$be   <- be
		Theta$ye   <- ye
		Theta$ye.sex <- ye.sex
		Theta$de   <- de
		Theta$bpr  <- bpr
		Theta$ypr  <- ypr
		Theta$spr  <- spr
		Theta$dpr  <- dpr
		Theta$wbar <- wbar
		
		return(Theta)
	})
}



# |---------------------------------------------------------------------------|
# | Calculate equilibrium values for a given fe vector              
# |---------------------------------------------------------------------------|
# | This function calls asem several times and constructs a data 
# | frame with columns: fe ye be re de ypr spr
.calcEquilibrium <- function(Theta, Scenario=NULL, bycatch=0)
{
	with(Theta, {
		# Proto-type function to get equilibrium vector.
		fn <- function(fe)
		{
			tmp <- .asem(fe, Theta, bycatch)
			out <- c(fe=fe, ye=tmp$ye, be=tmp$be, de=tmp$de, 
				re=tmp$re, spr=tmp$spr, ypr=tmp$ypr, 
				bpr=tmp$bpr, dpr=tmp$dpr,ye.sex=tmp$ye.sex)
				
			# average weight arrays
			wbar_f <- c(wbar=tmp$wbar[1,])
			if(dim(tmp$wbar)[1] > 2) 
				wbar_m <- c(wbar=tmp$wbar[2,])
			else
				wbar_m <- c(wbar=tmp$wbar[1,])

			out <- c(out, wbar_f=wbar_f, wbar_m=wbar_m, Scenario=Scenario)
			return(out)
		}
		
		xx <- sapply(fe, fn)
		Theta$equil <- as.data.frame(t(xx))
		
		return(Theta)
	})
}

# |---------------------------------------------------------------------------|
# | Construct a results data frame for use with ggplot.  
# |---------------------------------------------------------------------------|
# |
.makeDataFrame <- function(Theta)
{
	cat(".makeDataFrame\n")
	n  <- length(Theta)
	df <- data.frame()
	print(n)
	for(i in 1:n)
	{
		# | Find Fspr = 0.3
		xx      <- as.double(Theta$equil$spr)
		yy      <- as.double(Theta$equil$fe)
		fspr.30 <- approx(xx,yy,0.30,yleft=0,yright=max(yy))$y
		fspr.35 <- approx(xx,yy,0.35,yleft=0,yright=max(yy))$y
		
		# | Find F0.1
		xx      <- as.double(Theta$equil$bpr)
		yy      <- as.double(Theta$equil$fe)
		bpr0    <- 0.1*xx[1]
		f0.1    <- approx(xx, yy, bpr0, yright=max(yy))$y
		
		# | Find Fmsy
		yy      <- as.double(Theta$equil$ye)
		ii      <- which.max(yy)
		fmsy    <- as.double(Theta$equil$fe)[ii]
		msy     <- as.double(Theta$equil$ye)[ii]
		
		u1 <- round(1-exp(-f0.1), 3)
		u2 <- round(1-exp(-fspr.30), 3)
		u3 <- round(1-exp(-fmsy), 3)
		if(i==1)
		cat( "\t", "Fspr30", "\t", "Umsy\n")
		cat( "\t", u2, "\t", u3, "\n")
		
		# print(Theta)
		tmp <- data.frame(Stock  = Theta$RA, 
						Fspr.30 = fspr.30,
						Fspr.35 = fspr.35,
						F0.1    = f0.1,  
						Fmsy    = fmsy, 
						msy     = msy, 
						Theta$equil)
		df  <- rbind(df, tmp)
	}
	return(df)	
}

# |----------------------------------------------------------------------------------|
# | MAIN
# |----------------------------------------------------------------------------------|
# |
# Theta = T1
Theta <- .calcLifeTable(Theta)
Theta <- .calcSelectivities(Theta)
Theta <- .calcSRR(Theta)
Theta <- .calcEquilibrium(Theta)
TwoSex<- .makeDataFrame(Theta)

T1 <- .calcLifeTable(T1)
T1 <- .calcSelectivities(T1)
T1 <- .calcSRR(T1)
T1 <- .calcEquilibrium(T1)
OneSex<- .makeDataFrame(T1)

# RBind two data frames

DF <- rbind(TwoSex,OneSex)


# |----------------------------------------------------------------------------------|
# | GRAPHICS
# |----------------------------------------------------------------------------------|
# |

df.ye <- subset(DF,select=c(fe,ye,msy,Stock))
mdf   <- melt(df.ye,id.vars=c("fe","Stock"))
p     <- ggplot(df.ye,aes(fe,ye/msy,color=Stock)) + geom_line()
p.ye  <- p + labs(x="Fishing mortality",y="Equilibrium yeild relative to MSY",color="Model")
print(p.ye)

df.spr<- subset(DF,select=c(fe,spr,Stock))
mdf   <- melt(df.spr,id.vars=c("fe","Stock"))
p     <- ggplot(mdf,aes(fe,value,color=Stock)) + geom_line() + ylim(c(0,1))
p     <- p + geom_segment(aes(x=0,y=0.35,xend=0.310,yend=0.35))
p     <- p + geom_segment(aes(x=0.310,y=0,xend=0.310,yend=0.35))
p.spr <- p + labs(x="Fishing mortality",y="Female spawning potential ratio",color="Model")
print(p.spr)

df.ypr<- subset(DF,select=c(fe,ypr,Stock))
p     <- ggplot(df.ypr,aes(fe,ypr,color=Stock)) + geom_line()
p     <- p + labs(x="Fishing mortality",y="Yield per recruit",color="Model")
print(p)
