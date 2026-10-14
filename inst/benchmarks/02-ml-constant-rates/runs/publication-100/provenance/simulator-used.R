simulate <-
function (pars.tb, k, max.taxa = Inf, max.t1 = Inf, max.t2 = Inf, 
    x0, single.lineage = TRUE) 
{
    pars1 <- pars.tb[1, ]
    if (x0 < 1 || x0 > k) 
        stop("x0 must be an integer in [1,k]")
    pars1.list <- diversitree:::inflate.pars.classe(pars1, k)
    pars1.list$lambda[which(is.na(pars1.list$lambda))] <- 0
    to1 <- matrix(unlist(lapply(1:k, function(i) (1:k)[-i])), 
        k, k - 1, TRUE)
    pars1 <- cbind(rowSums(pars1.list$lambda), pars1.list$mu, 
        matrix(pars1[-seq_len(k * k * (k + 1)/2 + k)], k, k - 
            1, TRUE))
    r.i.1 <- rowSums(pars1)
    extinct <- FALSE
    split <- FALSE
    parent <- 0
    n.i <- rep(0, k)
    len <- 0
    t <- 0
    hist <- list()
    if (single.lineage) {
        states <- x0
        n.taxa <- lineages <- n.i[x0] <- 1
        start <- 0
    }
    else {
        stop("Nope.")
    }
    while (n.taxa <= max.taxa && n.taxa > 0) {
        r.n <- r.i.1 * n.i
        r.tot <- sum(r.n)
        dt <- rexp(1, r.tot)
        t <- t + dt
        if (t > max.t1) {
            dt <- dt - (t - max.t1)
            len[lineages] <- len[lineages] + dt
            t <- max.t1
            break
        }
        len[lineages] <- len[lineages] + dt
        state <- sample(k, 1, FALSE, r.n/r.tot)
        j <- sample(n.i[state], 1)
        lineage <- lineages[states[lineages] == state][j]
        type <- sample(k + 1, 1, FALSE, pars1[state, ])
        if (type == 1) {
            if (n.taxa == max.taxa) 
                break
            new.i <- length(extinct) + 1:2
            split[lineage] <- TRUE
            extinct[new.i] <- split[new.i] <- FALSE
            lam <- pars1.list$lambda[state, , ]
            j <- sample(k * k, 1, FALSE, lam)
            s.daught <- c((j - 1)%%k + 1, (j - 1)%/%k + 1)
            states[new.i] <- s.daught
            n.i[state] <- n.i[state] - 1
            n.i[s.daught[1]] <- n.i[s.daught[1]] + 1
            n.i[s.daught[2]] <- n.i[s.daught[2]] + 1
            parent[new.i] <- lineage
            start[new.i] <- t
            len[new.i] <- 0
            n.taxa <- n.taxa + 1
            lineages <- which(!split & !extinct)
        }
        else if (type == 2) {
            extinct[lineage] <- TRUE
            lineages <- which(!split & !extinct)
            n.i[state] <- n.i[state] - 1
            n.taxa <- n.taxa - 1
        }
        else {
            states[lineage] <- state.new <- to1[state, type - 
                2]
            n.i[c(state.new, state)] <- n.i[c(state.new, state)] + 
                c(1, -1)
            hist[[length(hist) + 1]] <- c(lineage, t, state, 
                state.new)
        }
    }
    info1 <- data.frame(idx = seq_along(extinct), len = len, 
        parent = parent, start = start, state = states, extinct = extinct, 
        split = split)
    attr(info1, "hist") <- .classe_td_history_table(hist, info1)
    attr(info1, "k") <- k
    tb.e1 <- diversitree:::me.to.ape.bisse(info1[-1, ], info1$state[1])
    if (n.taxa == 0) 
        return(NULL)
    t.epoch1 <- t
    pars2 <- pars.tb[2, ]
    pars2.list <- diversitree:::inflate.pars.classe(pars2, k)
    pars2.list$lambda[which(is.na(pars2.list$lambda))] <- 0
    to2 <- matrix(unlist(lapply(1:k, function(i) (1:k)[-i])), 
        k, k - 1, TRUE)
    pars2 <- cbind(rowSums(pars2.list$lambda), pars2.list$mu, 
        matrix(pars2[-seq_len(k * k * (k + 1)/2 + k)], k, k - 
            1, TRUE))
    r.i.2 <- rowSums(pars2)
    while (n.taxa <= max.taxa && n.taxa > 0) {
        r.n <- r.i.2 * n.i
        r.tot <- sum(r.n)
        dt <- rexp(1, r.tot)
        t <- t + dt
        if (t > max.t2) {
            dt <- dt - (t - max.t2)
            len[lineages] <- len[lineages] + dt
            t <- max.t2
            break
        }
        len[lineages] <- len[lineages] + dt
        state <- sample(k, 1, FALSE, r.n/r.tot)
        j <- sample(n.i[state], 1)
        lineage <- lineages[states[lineages] == state][j]
        type <- sample(k + 1, 1, FALSE, pars2[state, ])
        if (type == 1) {
            if (n.taxa == max.taxa) 
                break
            new.i <- length(extinct) + 1:2
            split[lineage] <- TRUE
            extinct[new.i] <- split[new.i] <- FALSE
            lam <- pars2.list$lambda[state, , ]
            j <- sample(k * k, 1, FALSE, lam)
            s.daught <- c((j - 1)%%k + 1, (j - 1)%/%k + 1)
            states[new.i] <- s.daught
            n.i[state] <- n.i[state] - 1
            n.i[s.daught[1]] <- n.i[s.daught[1]] + 1
            n.i[s.daught[2]] <- n.i[s.daught[2]] + 1
            parent[new.i] <- lineage
            start[new.i] <- t
            len[new.i] <- 0
            n.taxa <- n.taxa + 1
            lineages <- which(!split & !extinct)
        }
        else if (type == 2) {
            extinct[lineage] <- TRUE
            lineages <- which(!split & !extinct)
            n.i[state] <- n.i[state] - 1
            n.taxa <- n.taxa - 1
        }
        else {
            states[lineage] <- state.new <- to2[state, type - 
                2]
            n.i[c(state.new, state)] <- n.i[c(state.new, state)] + 
                c(1, -1)
            hist[[length(hist) + 1]] <- c(lineage, t, state, 
                state.new)
        }
    }
    info <- data.frame(idx = seq_along(extinct), len = len, parent = parent, 
        start = start, state = states, extinct = extinct, split = split)
    attr(info, "info.epoch1") <- info1
    attr(info, "t.epoch1") <- t.epoch1
    attr(info, "t.total") <- t
    attr(info, "sim.pars") <- pars.tb
    attr(info, "hist") <- .classe_td_history_table(hist, info)
    attr(info, "k") <- k
    info
}
