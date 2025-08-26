library(testthat)


test_that("Check 'method = dtl' works with parameter H0=FALSE", {
    expect_no_warning(mams(method = "dtl", K = c(4,1), H0 = FALSE))
})

test_that("Check summary works with parameter H0=FALSE", {
    expect_no_warning(mams(method = "dtl", K = c(4,1), H0 = FALSE))
})

test_that("dtl with K=3, J=2 returns N=282", {
    expect_equal(
        mams(method = "dtl", K = c(3, 1), delta = 0.545, delta0 = 0.178,
        sd = 1, p = NULL,  p0 = NULL)$N,
        282)
})

test_that("mams with K=4, J=2 returns N=364", {
    expect_equal(
        mams(method = "dtl", K = c(4, 1), delta = 0.545, delta0 = 0.178,
        sd = 1, p = NULL, p0 = NULL)$N,
        364)
})

test_that("mams with K=4, J=3 returns N=330", {
    expect_equal(
        mams(method = "dtl", K = c(4, 2, 1), J = 3, delta = 0.545, 
                delta0 = 0.178,
                sd = 1, p = NULL, p0 = NULL, r = 1:3, r0 = 1:3)$N,
        330)
})

test_that("mams with K=6, J=2 returns N=531", {
    expect_equal(mams(method = "dtl", K = c(6, 1), delta = 0.545,
        delta0 = 0.178, sd = 1, p = NULL, p0 = NULL)$N,
        531, tolerance = 9)
})

    test_that("mams with K=6, J=3 returns N=455", {
    expect_equal(mams(method = "dtl", K = c(6, 3, 1), J = 3, delta = 0.545, 
        delta0 = 0.178,sd = 1, p = NULL, p0 = NULL, r = 1:3, r0 = 1:3)$N,
        455)
})

test_that("mams with K=8, J=3 returns N=585", {
    expect_equal(mams(method = "dtl", K = c(8, 3, 1), J =3, delta = 0.545, 
        delta0 = 0.178, sd = 1, p = NULL, p0 = NULL, r = 1:3, r0 = 1:3)$N,
        585)
})

testthat::test_that("Separate with user-define delta and delta0", {
    testthat::expect_no_warning(
        mams(method = "dtl", K = c(4,1), ushape = "obf",
                    nsim = 1000, J = 2,
                    delta = 0.8, delta0 = 0, sd = 1, 
                    p = NULL, p0 = NULL)
    )
})

testthat::test_that("Separate with user-define p and p0", {
    testthat::expect_no_warning(
        mams(method = "dtl", ushape = "obf", K = c(4,1),
                    nsim = 1000, J = 2,
                    p = 0.73, p0 = 0.5)
    )
})

toto <- mams(method = "dtl", nsim = 1000, K = c(4,1))
toto5 <- mams(method = "dtl", K = c(5, 1))

testthat::test_that("Testing mams switching from prob.scale to cohen.d", {
        
    testthat::expect_no_warning(
        toto5delta <- mams(obj = toto, method = "dtl", K = c(5, 1), 
                                delta = 0.5, delta0 = 0, sd = 1, nsim = 1000)
    )
    
    testthat::expect_no_warning(
        toto5p <- mams(obj = toto5delta, p = 0.63, p0 = 0.5,
                            nsim = 1000)
    )
})

test_that("Simulation with default parameters except K=5", {
    
    expect_no_warning(toto.sim <- mams.sim(obj = toto5, nsim = 1000, 
                                            K = c(5,1)))
})

test_that("DTL simulation with specified `pv` parameters", {
        pv <- c(0.7, rep(0.5, 4))

    expect_silent(toto.sim.pv <- mams.sim(obj = toto5, 
                                                pv = pv, nsim = 1000))
    expect_equal(toto.sim.pv$input$par$pv, pv) 
})

test_that("DTL simulation with specified `deltav` parameters", {
    deltav <- c(0.5, rep(0, 4))
    
    expect_silent(toto.sim.deltav <- mams.sim(obj = toto5, 
                                                deltav = deltav, nsim = 1000))
    expect_equal(toto.sim.deltav$input$par$deltav, deltav)
})

test_that("Simulation with incorrect number of deltav parameters", {
    deltav <- c(0.5, rep(0, 3)) # Intentionally wrong number of reps
    
    expect_error(
        toto.sim.deltav.wrong <- mams.sim(obj = toto5, 
                                                deltav = deltav)
    )
})

test_that("mams error when provided object altered by mams_sim()", {
    toto.sim.deltav <- mams.sim(obj = toto5, 
                                    deltav = c(0.5, rep(0, 4)))
    expect_error(mams(obj = toto.sim.deltav, method = "dtl"))
})

test_that("Example workflow without warnings", {
    expect_no_warning({
        toto.sim <- mams.sim(obj = toto5, nsim = 1000)
        toto.sim.pv <- mams.sim(obj = toto5, 
                                    pv = c(0.7, rep(0.5, 4)), nsim = 1000)
        toto.sim.deltav <- mams.sim(obj = toto5, 
                                        deltav = c(0.5, rep(0, 4)), nsim = 1000)
    })
})

test_that("mams() an error when nsim is less than 1000", {
    expect_error(
        mams(method = "dtl", nsim = 999)
    )
})

test_that("mams handles additional parameters correctly", {
    expect_no_warning({
        toto5delta <- mams(obj = toto, method = "dtl", K = c(5, 1), 
        delta = 0.5, delta0 = 0, sd = 1, nsim = 1000)
        toto5p <- mams(obj = toto5delta, method = "dtl", p = 0.63, p0 = 0.5,
                            nsim = 1000)
    })
    
    expect_equal(toto5delta$input$par$delta, 0.5)
    expect_equal(toto5delta$input$par$delta0, 0)
    expect_equal(toto5delta$input$par$sigma[1], 1)
    expect_equal(toto5p$input$par$p, 0.63)
    expect_equal(toto5p$input$par$p0, 0.5)
})


test_that("Sequential simulations run without errors or warnings", {
    expect_no_warning(toto <- mams(method = "dtl", nsim = 1000, K=c(4,1)))

    expect_no_warning(toto.sim <- mams.sim(obj = toto, nsim = 1000))
    expect_no_warning(mams.sim(obj = toto.sim, nsim = 1000))
    expect_no_warning(toto.sim3 <- mams.sim(obj = toto.sim, 
                                                deltav = c(0.5, rep(0, 3)),
                                                nsim = 1000))
})

test_that("Parameter validation in mams_sim", {
    invalid_deltav <- "invalid_type"
    expect_error(
        mams.sim(obj = toto, deltav = invalid_deltav)
    )
})

test_that("Test summary function with complex pv or deltav", {
    deltav <- c(0.5, 0.4, 0, 0)
    toto.sim.deltav <- mams.sim(obj = toto, deltav = deltav, nsim = 1000)
    expect_no_warning(summary(toto.sim.deltav))

    pv <- c(0.73, 0.65, 0.5, 0.5)
    toto.sim.pv <- mams.sim(obj = toto, pv = pv, nsim = 1000)
    expect_no_warning(summary(toto.sim.pv))
})

test_that("Test automatic allocation for r and r0", {
    expect_error(toto1J  <- mams(toto, method = "dtl", J = 1))
})

test_that("Default call test for automatically allocation of r and r0", {   
    expect_s3_class(toto, "MAMS")
    expect_equal(toto$J, 2)
    expect_equal(toto$input$r, 1:2)
    expect_equal(toto$input$r0, 1:2)
})

test_that("Passing obj as an unnamed argument without specifying J", {
    expect_no_warning(result_unnamed <- mams(toto, method = "dtl",
                                                nsim = 1000, K=c(4,1)))
    
    expect_s3_class(result_unnamed, "MAMS")
    expect_equal(result_unnamed$J, 2)
    expect_equal(result_unnamed$input$r, 1:2)
    expect_equal(result_unnamed$input$r0, 1:2)
})


test_that("Overriding J and providing custom r=1:2 and r0=2:3", {
    expect_no_warning(result_custom_r <- mams(toto, method = "dtl", J = 2, 
                        r = 1:2, r0 = 2:3, nsim = 1000))
    
    expect_equal(result_custom_r$J, 2)
    expect_equal(result_custom_r$input$r, 1:2)
    expect_equal(result_custom_r$input$r0, 2:3)
})

test_that("Switch between methods", {
    expect_warning({
    
    mams.fit1 = mams(toto, method="dtl",ushape="pocock", nsim=1000)
    mams.fit2 = mams(mams.fit1,method = "simultaneous", nsim = 1000)
    
    })
})

test_that("mams function calls with method='dtl' 
                                                without predefined delta0", {
    expect_no_error({
    mams.fit1a <- mams(method = "sep", ushape = "pocock", nsim = 1000)
    mams.fit1b <- mams(mams.fit1a, ushape = "obf", nsim = 1000)
    suppressWarnings({
    mams.fit1e <- mams(mams.fit1b, J = 3, K = c(4,2,1), 
                                delta = c(0.25), sd = 0.5, nsim = 1000, 
                                method = "dtl")
    })
    })
})

test_that("Overriding J=2 without specifying r and r0", {
    expect_no_error({
    suppressWarnings({
mams.fit6a = mams(method = "dtl", J=3, r=1:3,r0=1:3, K=c(4,2,1),
                        nsim=1000)
mams.fit6b = mams(mams.fit6a, delta0=0.075)
mams.fit6b$sim$H1$main$ess[1,4] == 27
mams.fit6b2 = mams(mams.fit6a, p0=0.521)
mams.fit6b2$sim$H1$main$ess[1,4] == 27
    })
    })
})

test_that("mams checks for r0 and r vectors", {
    expect_error(
    mams(K = c(3,2,1), J = 3, alpha = 0.05, power = 0.9, r = c(1, 2), 
    r0 = c(1, 2), method = "dtl"),
    "Length of allocation ratios does not match number of stages."
    )

    expect_error(
    mams(K = c(3, 1), J = 2, alpha = 0.05, power = 0.9, r = c(1, 2), 
    r0 = c(3, 2), method = "dtl"),
    "`r0` must be a monotonically increasing vector."
    )

    expect_error(
    mams(K = c(3, 1), J = 2, alpha = 0.05, power = 0.9, r = c(3, 2), 
    r0 = c(1, 2), method = "dtl"),
    "`r` must be a monotonically increasing vector."
    )

    expect_error(
        mams(K = c(3, 1), J = 2, alpha = 0.05, power = 0.9, 
        r = c(0, 2), 
        r0 = c(1, 2), method = "dtl")
    )

    expect_error(
        mams(K = c(3, 2), J = 2, alpha = 0.05, power = 0.9, 
        r = c(1, 2), 
        r0 = c(0, 2), method = "dtl")
    )

    expect_error(
        mams(K = c(3, 1), J = 2, alpha = 0.05, power = 0.9, 
        r = c(1, 2), 
        r0 = c(1.5, 2), method = "dtl"),
        "First element of `r0` must be integers."
    )
})

test_that("Peter Greenstreet's example of allocation ratios", {
toto_1 <- mams(
    K = c(3, 1), J = 2,
    power = 0.8, r = c(3, 7), r0 = c(3, 7), ushape = "obf",
    lshape = "fixed", lfix = -9999, nstart = 1, nsim = 1000, method = "dtl"
)
expect_equal(toto_1$n, 9)
expect_equal(unname(toto_1$rMat[1,1]), 1)
expect_equal(unname(toto_1$rMat[1,2]), 2.333333, tolerance = 0.01)

})

### Extra tests for input parameters checks

test_that("sep method: allocation ratio vector length is checked", {
    expect_error(
    mams(K = c(4,1), J = 2, alpha = 0.05, power = 0.9, r = c(1), r0 = c(1), 
            method = "dtl"),
    "Length of allocation ratios does not match number of stages."
    )
    expect_error(
    mams(K = c(4,1), J = 2, alpha = 0.05, power = 0.9, r = c(1, 2, 3), 
    r0 = c(1,2,3), 
    method = "dtl"),
    "Length of allocation ratios does not match number of stages."
    )
})

test_that("sep method: allocation ratio vectors must be monotonically increasing
    and valid", {
# r0 not monotonically increasing
    expect_error(
    mams(K = c(4,1), J = 2, alpha = 0.05, power = 0.9, r = c(1, 2), 
                r0 = c(2, 1),
                method = "dtl"),
    "`r0` must be a monotonically increasing vector."
    )
# r values must be >= 1
    expect_error(
    mams(K = c(4,1), J = 2, alpha = 0.05, power = 0.9, r = c(0, 2), 
                r0 = c(1, 2),
                method = "dtl"),
    "`r` values must be >= 1."
    )
# r0[1] must be >= 1
    expect_error(
    mams(K = c(4,1), J = 2, alpha = 0.05, power = 0.9, r = c(1, 2), 
                r0 = c(0, 2),
                method = "dtl"),
    "`r0\\[1\\]` must be >= 1."
    )
# All elements of r0 must be integers (via first element check message)
    expect_error(
    mams(K = c(4,1), J = 2, alpha = 0.05, power = 0.9, r = c(1,2), 
                r0 = c(1.5, 2), method = "dtl"),
    "First element of `r0` must be integers."
    )
})

test_that("sep method: effect size parameter check", {
# Both p and delta, sd provided => error
    expect_warning(expect_error(
    mams(K = c(4,1), J = 2, alpha = 0.05, power = 0.9,
                r = 1:2, r0 = 1:2, 
                p = 0.75, delta = 0.5, sd = 1, method = "dtl")
    ))

# p out-of-bound error: p < 0 or > 1
    expect_error(
    mams(K = c(4,1), J = 2, alpha = 0.05, power = 0.9,
            r = 1:2, r0 = 1:2, 
            p = -0.1, method = "dtl"),
    "Treatment effect parameter not within 0 and 1."
)

# sd must be positive when using delta scale.
    expect_warning(expect_error(
    mams(K = c(4,1), J = 2, alpha = 0.05, power = 0.9,
            r = 1:2, r0 = 1:2, 
            delta = 0.5, sd = -1, method = "dtl"),
    "Standard deviation must be positive."
    ))
})

test_that("sep method: alpha must be between 0 and 1", {
    expect_error(
        mams(K = c(4,1), J = 2, alpha = -0.05, power = 0.9, r = 1:2, 
                    r0 = 1:2, 
                    method = "dtl"),
        "Error rate or power not between 0 and 1."
    )
    expect_error(
    mams(K = c(4,1), J = 2, alpha = 1.05, power = 0.9, r = 1:2, r0 = 1:2, 
    method = "dtl"),
    "Error rate or power not between 0 and 1."
    )
})

test_that("K=3, J=2 yields correct power, type I error, and ESS", {
  design1a <- mams(method = "dtl", K = c(3,1), J = 2, alpha = 0.05, 
                    power = 0.8,  delta = 2, delta0 = 0, sd = 10)

  # The power should be ≈ 0.8 
  eff_H1 <- as.numeric(design1a$sim$H1$main$efficacy[5, 2])
  expect_equal(eff_H1, 0.8, tolerance = 0.02)

  # The p-value should be ≈ 0.05
  eff_H0 <- as.numeric(design1a$sim$H0$main$efficacy[4, 2])
  expect_equal(eff_H0, 0.05, tolerance = 0.01)

  # The ESS should be 420
  ess_H1 <- design1a$sim$H1$main$ess[1, 1]
  expect_equal(ess_H1, 420)
})

test_that("DTL method with K=c(6,1) works without cumulative sample size issues", {
# Issue: DTL method showing cumulative sample size at stage 1 (n=798) 
# larger than final sample size at stage 2 (n=570)

    skip_on_cran()  # This test may take time

    expect_no_error({
    fit.mams <- mams(K=c(6,1),J=2,delta=0.545,delta0=0.178,sd=1,alpha=0.05,
                        power=0.9, method="dtl",r=c(1,2),r0=2:3)
})

    captured.output <- capture.output(summary(fit.mams))
    expect_true(grep(456, captured.output[29]) == 1)
})