library(testthat)
test_that("Compare output of mams() to m1 in JSS Jaki et al", {
        m1 <- mams(K = 3, J = 1, p = 0.65, p0 = 0.55, r = 1, r0 = 1,
                                        alpha = 0.05, power = 0.9)

                expect_equal(m1$u, 2.062, tolerance = 0.001)
                expect_equal(m1$l, 2.062, tolerance = 0.001)
                expect_equal(m1$N, 316)
                expect_equal(m1$n, 79)
})

test_that("Compare output of mams() to m1d in JSS Jaki et al", {
                m1d <-  mams(K = 3, J = 1, p = NULL, p0 = NULL, 
                                delta = 0.545, delta0 = 0.178, sd = 1, r = 1, 
                                r0 = 1, alpha = 0.05, power = 0.9)

                expect_equal(m1d$u, 2.062, tolerance = 0.001)
                expect_equal(m1d$l, 2.062, tolerance = 0.001)
                expect_equal(m1d$N, 316)
                expect_equal(m1d$n, 79)
})

test_that("Compare output of mams() to m2 in JSS Jaki et al", {
        m2 <-  mams(K = 3, J = 2, p = 0.65, p0 = 0.55, r = 1:2,
                        r0 = c(2, 4),
                        alpha = 0.05, power = 0.9, ushape = "triangular",
                        lshape = "triangular")

                expect_equal(m2$u, c(2.359, 2.225), tolerance = 0.001)
                expect_equal(m2$l, c(0.786, 2.225), tolerance = 0.001)
                expect_equal(m2$N, 380)
                expect_equal(m2$n, 76) #! FIXME (was 38)
                expect_equal((m2$rMat * m2$n)[1:2,],
                                matrix(c(76,152, 38,76), nrow = 2, ncol = 2,
                                byrow = TRUE, dimnames = list(c("Control",
                                        "T1"), c("Stage 1", "Stage 2"))))
})

test_that("Compare output of mams() to m3 in JSS Jaki et al", {
                m3 <- mams(K = 3, J = 3, p = 0.65, p0 = 0.55, alpha = 0.05,
                                power = 0.9, r = 1:3, r0 = 1:3,
                                ushape = function(x) return(x:1),
                                lshape = "fixed", lfix = 0)

                expect_equal(m3$u, c(6.125, 4.084, 2.042), tolerance = 0.001)
                expect_equal(m3$l, c(0.000, 0.000, 2.042), tolerance = 0.001)
                expect_equal(m3$N, 324)
                expect_equal(m3$n, 27)
                expect_equal((m3$rMat * m3$n)[1:2,],
                                matrix(c(27,54,81, 27,54,81),
                                nrow = 2, ncol = 3, byrow = TRUE,
                                dimnames = list(c("Control",
                                        "T1"), c("Stage 1", "Stage 2", 
                                        "Stage 3"))))
})

test_that("Compare output of mams.sim() to m2sim in JSS Jaki et al", {
        m2 <-  mams(K = 3, J = 2, p = 0.65, p0 = 0.55, r = 1:2,
                                r0 = c(2, 4), alpha = 0.05, power = 0.9,
                                ushape = "triangular",
                                lshape = "triangular")

        suppressWarnings(m2sim <- mams.sim(nsim = 1e5,
                                                nMat = t(m2$n * m2$rMat),
                                                u = m2$u, l = m2$l,
                                                pv = rep(0.5, 3), ptest = 1:2))

        expect_equal(round(m2sim$sim$H0$main$efficacy[4,2],3), 0.050,
                                                        tolerance = 0.01)
        expect_equal(round(m2sim$sim$H0$main$efficacy[5,2],3), 0.016,
                                                        tolerance = 0.02)
        expect_equal(round(m2sim$sim$H0$main$efficacy[7,2],3), 0.034,
                                                        tolerance = 0.04)
                expect_equal(round(sum(m2sim$sim$H0$main$ess[1]),3), 
                                                        244.907,
                                                        tolerance = 0.002)
        })

test_that("Compare output of Pocock, O’Brien-Fleming, and triangular boundaries
                with  Table 1 in JSS Jaki et al", {

        poc <- mams(K = 3, J = 3, p = 0.65, p0 = 0.55, r = 1:3, r0 = 1:3,
                        alpha = 0.05, power = 0.9, ushape = "pocock",
                        lshape = "pocock")

        obf <- mams(K = 3, J = 3, p = 0.65, p0 = 0.55, r = 1:3, r0 = 1:3,
                        alpha = 0.05, power = 0.9, ushape = "obf",
                        lshape = "obf")

        tri <- mams(K = 3, J = 3, p = 0.65, p0 = 0.55, r = 1:3, r0 = 1:3,
                        alpha = 0.05, power = 0.9, ushape = "triangular",
                        lshape = "triangular")

        suppressWarnings(pocsim <- mams.sim(nsim = 1e5,
                         nMat = t(poc$n * poc$rMat), u = poc$u,
                        l = poc$l, pv = c(0.65, rep(0.55, 2)), ptest = 1))

        suppressWarnings(obfsim <- mams.sim(nsim = 1e5,
                        nMat = t(obf$n * obf$rMat), u = obf$u,
                        l = obf$l, pv = c(0.65, rep(0.55, 2)), ptest = 1))

        suppressWarnings(trisim <- mams.sim(nsim = 1e5,
                         nMat = t(tri$n * tri$rMat), u = tri$u,
                        l = tri$l, pv = c(0.65, rep(0.55, 2)), ptest = 1))

                expect_equal(pocsim$N, 396)
                expect_equal(obfsim$N, 336)
                expect_equal(trisim$N, 408)

        expect_equal(round(
        sum(pocsim$sim[["H1"]]$main$ess[,"ess"]),1), 232.2,
                                                                tolerance = 0.1)
        expect_equal(
        round(sum(obfsim$sim[["H1"]]$main$ess[,"ess"]),1), 259.1,
                                                                tolerance = 0.1)
        expect_equal(
        round(sum(trisim$sim[["H1"]]$main$ess[,"ess"]),1), 217.8,
                                                                tolerance = 0.1)
})

test_that("Compare output of new.bounds() with m2.nb in JSS Jaki et al", {
                m2 <-  mams(K = 3, J = 2, p = 0.65, p0 = 0.55,
                r = 1:2, r0 = c(2, 4),
                alpha = 0.05, power = 0.9,
                ushape = "triangular",
                lshape = "triangular")

                m2.nb <- new.bounds(K = 3, J = 2,
                        nMat = matrix(c(75, 152, 40, 76, 35, 76, 41, 76),
                        nrow = 2, ncol = 4), alpha = 0.05, u = m2$u[1],
                        l = m2$l[1], ushape = "triangular",
                        lshape = "triangular")

                expect_equal(m2$u, c(2.359, 2.224), tolerance = 0.001)
                expect_equal(m2$l, c(0.786, 2.224), tolerance = 0.001)
})

test_that("Compare output of m2.all and Appendix A, JSS Jaki et al", {
                m2 <-  mams(K = 3, J = 2, p = 0.65, p0 = 0.55,
                                r = 1:2, r0 = c(2, 4),
                                alpha = 0.05, power = 0.9,
                                ushape = "triangular",
                                lshape = "triangular")

                m2.all <- stepdown.mams(
                                nMat = matrix(c(76, 152,rep(c(38, 76), 3)),
                                nrow = 2, ncol = 4), lb = m2$l[1],
                                alpha.star = c(0.026, 0.05),
                                selection = "all.promising")

                expect_equal(m2.all$sample.sizes,
                                matrix(c(76,38,38,38, 152,76,76,76),
                                nrow = 2, ncol = 4, byrow = TRUE))

                expect_equal(matrix(unlist(m2.all$u)),
                                matrix(c(1.94, 1.72, 1.94, 1.72, 2.21, 2.06,
                                        1.94, 1.72, 2.21, 2.06, 2.21, 2.06,
                                        2.36, 2.22), nrow = 14,
                                ncol = 1, byrow = TRUE), tolerance = 0.01)

                expect_equal(matrix(unlist(m2.all$l)),
                        matrix(c(0.7864987, 1.72, 0.7864987, 1.72, 0.7864987,
                                2.06, 0.7864987, 1.72, 0.7864987, 2.06,
                                0.7864987, 2.06, 0.7864987, 2.22),
                                nrow = 14, ncol = 1, byrow = TRUE),
                                tolerance = 0.01)

                expect_equal(matrix(unlist(m2.all$`alpha.star`)),
                                matrix(c(0.026, 0.05, 0.026, 0.05, 0.026, 0.05,
                                        0.026, 0.05, 0.026, 0.05, 0.026, 0.05,
                                        0.026, 0.05),
                                nrow = 14, ncol = 1, byrow = TRUE))
})

test_that("Compare output of m2.best and Appendix A, JSS Jaki et al", {
                m2 <-  mams(K = 3, J = 2, p = 0.65, p0 = 0.55,
                                r = 1:2, r0 = c(2, 4),
                                alpha = 0.05, power = 0.9,
                                ushape = "triangular",
                                lshape = "triangular")

                m2.best <- stepdown.mams(
                                nMat = matrix(c(76, 152, rep(c(38, 76), 3)),
                                nrow = 2, ncol = 4),
                                lb = m2$l[1], alpha.star = c(0.026, 0.05),
                                selection = "select.best")

                expect_equal(m2.best$sample.sizes,
                        matrix(c(76, 152, rep(c(38, 76), 3)),
                                nrow = 2, ncol = 4))

                expect_equal(matrix(unlist(m2.best$u)),
                        matrix(c(1.94, 1.71, 1.94, 1.71, 2.21, 2.02,
                                        1.94, 1.71, 2.21, 2.02, 2.21, 2.02,
                                        2.36, 2.17),
                                nrow = 14, ncol = 1, byrow = TRUE))

                expect_equal(matrix(unlist(m2.best$l)),
                        matrix(c(0.7864987, 1.71, 0.7864987, 1.71, 0.7864987,
                        2.02, 0.7864987, 1.71, 0.7864987, 2.02, 0.7864987,
                                2.02, 0.7864987, 2.17),
                                nrow = 14, ncol = 1, byrow = TRUE))

                expect_equal(matrix(unlist(m2.best$`alpha.star`)),
                        matrix(c(0.026, 0.05, 0.026, 0.05, 0.026, 0.05,
                                        0.026, 0.05, 0.026, 0.05, 0.026, 0.05,
                                        0.026, 0.05),
                                nrow = 14, ncol = 1, byrow = TRUE))
})

test_that("Compare output of m2.update and Appendix B, JSS Jaki et al", {

                m2 <-  mams(K = 3, J = 2, p = 0.65, p0 = 0.55,
                                r = 1:2,r0 = c(2, 4),
                                alpha = 0.05, power = 0.9,
                                ushape = "triangular",
                                lshape = "triangular")

                m2.all <- stepdown.mams(
                                nMat = matrix(c(76, 152, rep(c(38, 76), 3)),
                                        nrow = 2, ncol = 4), lb = m2$l[1],
                                        alpha.star = c(0.026, 0.05),
                                        selection = "all.promising")

                m2.update <- stepdown.update(current.mams = m2.all,
                                                nobs = c(75, 40, 35, 41),
                                                zscores = c(1.1, 0.9, 0.9),
                                                selected.trts = c(1, 3),
                                        nfuture = matrix(c(228, 114, 35, 114),
                                                        nrow = 1, ncol = 4))
                expect_equal(m2.update$sample.sizes,
                                matrix(c(75,40,35,41, 228,114,35,114),
                                nrow = 2, ncol = 4, byrow = TRUE))

                expect_equal(matrix(unlist(m2.update$u)),
                                matrix(c(1.94, 1.73, 1.94, 1.71, 2.21, 1.92,
                                        1.94, 1.79, 2.21, 2.14, 2.21, 1.90,
                                        2.36, 2.22),
                                nrow = 14, ncol = 1, byrow = TRUE))

                expect_equal(matrix(unlist(m2.update$l)),
                        matrix(c(0.7864987, 1.73, 0.7864987, 1.71, 0.7864987,
                                        1.92, 0.7864987, 1.79, 0.7864987, 2.14,
                                        0.7864987, 1.90, 0.7864987, 2.22),
                                        nrow = 14, ncol = 1, byrow = TRUE))

                expect_equal(matrix(unlist(m2.update$`alpha.star`)),
                        matrix(c(0.00, 0.08843773, 0.00, 0.06896401, 0.00,
                                0.05613374, 0.00, 0.0576666, 0.00, 0.0510653,
                                        0.00, 0.04326494, 0.00, 0.04109772),
                                        nrow = 14, ncol = 1, byrow = TRUE),
                                        tolerance = 0.001)
})

test_that("Compare output of Non-normal endpoints in JSS Jaki et al", {
                prob <- c(0.075, 0.182, 0.319, 0.243, 0.015, 0.166)
                mord <- ordinal.mams(prob = prob, or = 3.06, or0 = 1.32,
                                        K = 3, J = 2,
                                alpha = 0.05, power = 0.9, r = 1:2, r0 = 1:2,
                                        ushape = "triangular",
                                        lshape = "triangular")

                expect_equal(mord$n, 34, tolerance = 0.001)
                expect_equal(mord$u, c(2.330, 2.197), tolerance = 0.001)
                expect_equal(mord$l, c(0.777, 2.197), tolerance = 0.001)

                mtite <- tite.mams(hr = 1.5, hr0 = 1.1,
                                K = 3, J = 2, alpha = 0.05,
                                power = 0.9, r = 1:2, r0 = 1:2,
                                ushape = "triangular",
                                lshape = "triangular")

                expect_equal(mtite$n, 81, tolerance = 0.001)
                expect_equal(mtite$u, c(2.330, 2.197), tolerance = 0.001)
                expect_equal(mtite$l, c(0.777, 2.197), tolerance = 0.001)
        })
