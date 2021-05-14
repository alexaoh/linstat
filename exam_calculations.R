# Multiple choice 7.
x <- matrix(c(1,1,1,1,0,0,1,1), ncol = 2, byrow = F)
h <- x%*%solve(t(x)%*%x)%*%t(x)
y <- c(-1,0,3,5)
sse <- t(y)%*%(diag(4)-h)%*%y
sse

df <- data.frame(cbind(y, c(0,0,1,1)))
colnames(df) <- c("resp", "covar")

fit <- lm(resp~as.factor(covar), data = df)
sse2 <- (4-2)*summary(fit)$sigma^2
sse2

# Problem 2b).
options(contrasts=c("contr.treatment", "contr.poly"))
bloodcelldata<-read.table("https://www.math.ntnu.no/emner/TMA4267/2021v/blodcells.csv")
fit1 <- lm(bloodcells~.,data = bloodcelldata)
model.matrix(fit1)
sse <- (105-8)*summary(fit1)$sigma^2
fitred<-lm(bloodcells~.-sport,data=bloodcelldata)
anova(fitred, fit1)

# Manually
r<-4 # rank of c
C2<-cbind(rep(0,r), rep(0,r), rep(0,r), rep(0,r), c(1,0, 0, 0), c(0,1,0,0), c(0,0,1,0), c(0,0,0,1))
d <- matrix(rep(0,r), ncol=1)

betahat3 <- matrix(fit1$coefficients, ncol=1)
sigma2hat3 <- summary(fit1)$sigma^2
X3 <- model.matrix(fit1)
F3 <- (t(C2%*%betahat3-d)%*%solve(C2%*%solve(t(X3)%*%X3)%*%t(C2))%*%(C2%*%betahat3-d))/(r*sigma2hat3)
F3

sse0 <- (105-4)*summary(fitred)$sigma^2
sse0 # f manually
((sse0-sse)/r)/(sse/(97))
