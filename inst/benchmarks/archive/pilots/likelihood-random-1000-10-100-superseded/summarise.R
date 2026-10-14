# Numerical agreement and compact publication plots ----
checks <- read.csv("results/likelihood-checks.csv")
summary <- data.frame()
for (model in unique(checks$comparison)) {
  for (survival in c(FALSE, TRUE)) {
    z <- checks[checks$comparison == model & checks$survival == survival, ]
    finite <- is.finite(z$difference)
    summary <- rbind(summary, data.frame(comparison = model, survival = survival,
      n = nrow(z), finite = sum(finite), passed = sum(z$pass),
      max.absolute.error = max(abs(z$difference[finite])),
      median.absolute.error = median(abs(z$difference[finite]))))
  }
}
write.csv(summary, "summary.csv", row.names = FALSE)
pdf("plots/likelihood-agreement.pdf", width = 10, height = 6)
par(mfrow = c(2,3), mar = c(5.5,5.5,2,1), mgp = c(4,1,0), las = 1, bty = "l")
for (model in unique(checks$comparison)) {
  z <- checks[checks$comparison == model & checks$survival, ]
  plot(z$reference, z$episodic, pch = 16, cex = .4, col = "#0072B250",
    xlab = "diversitree log likelihood", ylab = "DivExtra log likelihood", main = model)
  abline(0,1, col = "grey40", lty = 2)
}
for (model in unique(checks$comparison)) {
  z <- checks[checks$comparison == model & checks$survival, ]
  plot(z$n.tips, pmax(abs(z$difference), 1e-16), log = "y", pch = 16,
    cex = .4, col = "#0072B250", xlab = "Number of tips",
    ylab = "Absolute difference", main = model)
  abline(h = 1e-5, col = "#D55E00", lty = 2)
}
dev.off()
report <- c("# Fixed-parameter likelihood validation", "",
  sprintf("%d of %d paired checks passed the prespecified absolute tolerance of 1e-5.",
    sum(checks$pass), nrow(checks)), "",
  sprintf("Maximum finite absolute difference: %.8g.", max(abs(checks$difference), na.rm = TRUE)),
  "", "Results by model and survival conditioning:", "", "```",
  capture.output(print(summary, row.names = FALSE)), "```", "",
  "These checks establish numerical equivalence only over the tested cases. They do not validate genuinely changing epoch rates, parameter recovery, or stochastic mapping.",
  "Both epochs use identical parameters. Root treatment is ROOT.OBS, with complete sampling. The main plots show survival-conditioned results; the tables retain both settings.")
writeLines(report, "report.md")
print(summary)
