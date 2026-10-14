# Run after all inference tasks finish ----
dir.create("summary",showWarnings=FALSE)
dir.create("plots",showWarnings=FALSE)
m <- read.csv("manifest.csv")
groups <- c("50-99","100-249","250-499","500-999","1000-2000")
scenarios <- c("Slow","Medium","Fast")
rates <- c(Slow=.01,Medium=.1,Fast=.3)
parameters <- c("s_recent","d_recent","x_recent")
starts <- data.frame()
boundary <- data.frame()
for(i in m$tree.index) {
  task <- m$task[i]
  position <- (i-1) %% 10
  for(start in 1:3) {
    j <- position*3+start
    path <- sprintf("results/task-%03d/fit-%03d.rds",task,j)
    f <- if(file.exists(path)) readRDS(path) else list(error="missing")
    ok <- is.null(f$error) && length(f$lnLik)==1 && is.finite(f$lnLik) && f$lnLik > -1e49 &&
      length(f$par)==3 && all(parameters %in% names(f$par)) && all(is.finite(f$par))
    row <- data.frame(tree.index=i,start=start,valid=ok,timeout=isTRUE(f$timeout),
      error=if(is.null(f$error)) "" else f$error,lnLik=if(ok) f$lnLik else NA_real_,
      convergence=if(ok) f$convergence else NA_integer_,elapsed=if(ok) f$elapsed else NA_real_)
    for(p in parameters) row[[p]] <- if(ok) unname(f$par[p]) else NA_real_
    starts <- rbind(starts,row)
    if(!is.null(f$boundary)) boundary <- rbind(boundary,data.frame(tree.index=i,start=start,
      convergence=f$boundary$convergence,lnLik=f$boundary$lnLik,gain=f$boundary$lnLik-f$lnLik))
  }
}
write.csv(starts,"summary/all-starts.csv",row.names=FALSE)
write.csv(boundary,"summary/boundary-diagnostics.csv",row.names=FALSE)
best <- data.frame()
for(i in m$tree.index) {
  a <- starts[starts$tree.index==i & starts$valid,]
  if(nrow(a)==0) next
  b <- a[which.max(a$lnLik),]
  b$n.valid <- nrow(a)
  b$ll.spread <- diff(range(a$lnLik))
  b$strict <- nrow(a)==3 && all(a$convergence==0) && b$ll.spread<1e-4
  best <- rbind(best,b)
}
stopifnot(nrow(best)>0)
best <- merge(best,m,by="tree.index")
write.csv(best,"summary/best-fits.csv",row.names=FALSE)
write.csv(best[!best$strict,],"summary/flagged-fits.csv",row.names=FALSE)
metrics <- data.frame()
for(scenario in scenarios) for(group in groups) for(p in parameters) {
  a <- best[best$scenario==scenario & best$tip.group==group,]
  truth <- unname(rates[scenario]) * if(p=="x_recent") .5 else 1
  error <- abs(a[[p]]/truth-1)
  q <- quantile(a[[p]],c(.025,.5,.975),names=FALSE)
  metrics <- rbind(metrics,data.frame(scenario=scenario,tip.group=group,parameter=p,
    expected=100,n=nrow(a),truth=truth,mean=mean(a[[p]]),median=q[2],lower=q[1],upper=q[3],
    median.relative.error=median(error),within25.of.valid=mean(error<=.25),
    within25.of.expected=sum(error<=.25)/100,near.zero=sum(a[[p]]<1e-5)))
}
write.csv(metrics,"summary/recovery.csv",row.names=FALSE)

# Full-range and focused boxplots; all observations enter box calculations ----
for(view in c("full","focused")) {
  png(paste0("plots/recovery-",view,".png"),width=2700,height=2400,res=250,type="cairo")
  par(mfrow=c(3,3),mar=c(5,4,3,1),oma=c(0,0,2,0))
  for(scenario in scenarios) for(p in parameters) {
    a <- best[best$scenario==scenario,]
    truth <- unname(rates[scenario]) * if(p=="x_recent") .5 else 1
    values <- split(a[[p]],factor(a$tip.group,levels=groups))
    upper <- if(view=="focused") truth*3 else max(a[[p]],truth)*1.1
    boxplot(values,col="#A6CEE3",ylim=c(0,upper),pch=16,cex=.35,las=2,cex.axis=.8,
      main=paste(scenario,sub("_recent","",p)),ylab="Estimated rate")
    abline(h=truth,col="#B2182B",lty=2,lwd=2)
    if(view=="focused") for(g in 1:5) {
      number <- sum(values[[g]]>upper)
      if(number>0) text(g,upper*.94,paste0(number," above"),cex=.65)
    }
  }
  mtext("Constant rates in two linked epochs | red dashed line = truth",outer=TRUE)
  dev.off()
}
png("plots/recovery-vs-tips.png",width=2100,height=800,res=190,type="cairo")
par(mfrow=c(1,3),mar=c(5,4,3,1))
colors <- c("#0072B2","#E69F00","#009E73")
for(p in parameters) {
  a <- metrics[metrics$parameter==p,]
  plot(NA,xlim=c(1,5),ylim=c(0,max(30,100*a$median.relative.error,na.rm=TRUE)*1.1),
    xaxt="n",xlab="",ylab="Median relative error (%)",main=sub("_recent","",p))
  axis(1,1:5,labels=groups,las=2,cex.axis=.8)
  for(k in 1:3) {
    b <- a[a$scenario==scenarios[k],]
    lines(match(b$tip.group,groups),100*b$median.relative.error,type="b",col=colors[k],pch=16,lwd=2)
  }
  legend("topright",legend=scenarios,col=colors,pch=16,lty=1,bty="n",cex=.7)
}
dev.off()

# Short report, always stating completeness and the meaning of intervals ----
lines <- c("# Constant-rate recovery: publication benchmark", "",
  sprintf("Valid fits: %d/4500; optimizer code 0: %d; timeouts: %d.",sum(starts$valid),
    sum(starts$convergence==0,na.rm=TRUE),sum(starts$timeout)),
  sprintf("Trees with a valid fit: %d/1500; strict three-start likelihood agreement: %d/%d.",
    nrow(best),sum(best$strict),nrow(best)), "",
  "Median absolute relative error (%) in the largest size group:", "",
  "| Scenario | s | d | x |", "|---|---:|---:|---:|")
for(scenario in scenarios) {
  a <- metrics[metrics$scenario==scenario & metrics$tip.group=="1000-2000",]
  e <- a$median.relative.error[match(parameters,a$parameter)]*100
  lines <- c(lines,sprintf("| %s | %.1f | %.1f | %.1f |",scenario,e[1],e[2],e[3]))
}
lines <- c(lines,"", "The highest finite primary likelihood is used per tree; boundary diagnostics are separate. Strict agreement requires three converged starts and a likelihood spread below 0.0001. Convergence is not proof of accurate recovery or a global optimum.","",
  "Both epochs share s, d and x, with a known boundary at half the prescribed stem duration. This tests the episodic implementation under constant rates, not independent epoch-rate recovery. Slow/Medium/Fast are proportional time-scaled processes. Results are conditional on size selection and the living-lineage cap; no crossing or tip-state filter was imposed.","",
  "The lower/upper columns are central 95% replicate ranges, NOT confidence intervals or coverage. Focused plots retain all data in calculations and count estimates above the displayed range. See README.md for methods and summary/ for all numerical results.")
writeLines(lines,"report.md")
