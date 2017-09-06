#run mcmc to estimate peak locations from a given ADAM competition

## to test:
# source("../mcselpeak/mcselpeak.R")
# all.pars = mcmc.adam.peak.v1(test.data, 5000)

#Comment to check that git works

# pull in prerequisites
#library("mvtnorm")
library(DierckxSpline)


#some universal constants
#GENOME.LENGTH = 4639675
GENOME.LENGTH = 4000
#GENOME.LENGTH = 1000
NUM.MOVE.TYPES = 13

## parameters controlling move sizes -- tune for optimal acceptance
peak.move.sigma = 10
#peak.move.p = 0.05
peak.size.change.sigma = 10
peak.height.change.sigma = 1
mu0.change.sigma=0.05
sigma.base.change.sigma = 0.05
sigma.divergence.change.sigma=0.01
nu.change.sigma=0.5

peaksplit.delta.sd = 100
peaksplit.epsilon.sd = .1
peaksplit.gamma.sd = 0.1

## parameters controlling the prior distributions
#min.peak.sigma = 30
#max.peak.sigma = 250
#peak.sigma.prior.int = log(max.peak.sigma) - log(min.peak.sigma)


## Unneeded, given that we will be instead tracking binding energy (alimf)
height.prior.shape = 5
height.prior.scale=1

##Priors for binding energy 
energy.prior.shape = 1.5
energy.prior.scale = 18
mismatch.prior.shape = 3.125
mismatch.prior.scale = 0.8

##Prior for Peak Width
peak.sigma.shape = 1000
peak.sigma.size = 0.1


##Vector of bases to choose from
bases.vector = c("A","T","G","C","R","Y","W","S","M","K","B","V","H","D","N")

peak.height.min = 0.4
log.peak.height.prior.int = log(1-pgamma(peak.height.min, 1, scale=1))

#Not tracking peak number anymore, instead switching to TF number
#npeak.lambda=3

#In first implementation, tf.lamda will stand for static number of TFs.
#In later implementations, this will be shifted to a Poisson distribution with 
#tf.lamda as the Poisson variable
tf.lamda = 100

# v1 -- here we assume a fixed number of peaks and just have to fit their locations
# the background is independent student-t distributed, with a common mean for all elements
# and sd varying linearly with divergence (=1-homology)
# I *know* this is a bad model -- just need to implement the sampler first and then
# i can improve on it

mcmc.adam.peak.v1 = function(input.data, n.iter, pars.init=NULL, sparse=0) {
  # do mcmc to find locations of one or more peaks
  # before being passed to this function, data should be regularized to have a mean 0 
  #  and sd of 1
  # IT MUST ALSO BE SORTED IN SPATIAL ORDER ON THE CHROMOSOME
  #  Failure to adhere to these requirements may make your research explode

  # input.data must contain ***
  # n.iter is the number of iterations to do
  # sparse is the number of conformations between each saved conformation

  # we thus model the density at each point as D ~ T(Mu, Sigma, nu) 
  # Mu[i] is defined as mu_0 + sum(peaks=k) peak.heights[k] * dN(i - peak.centers[k], peak.sigmas[k])
  # peak.centers --  a set of peak centers
  # peak.heights --  a set of peak heights (can be positive or negative)
  # peak.sigmas -- a set of peak widths
  # sigma.base -- the intercept of the sigma_probe vs homology regression
  # sigma.divergence.r -- slope of the relationship between divergence and sigma
  # mu_0 is defined as the grand mean of all probe background values
  # nu is the number of degrees of freedom in the student-t background distribution

  # initialize the model parameters
  # we use decent guesses based on the construction of the problem
  # the only key here is that the peak centers are randomized. This is likely to be the
  #  slowest converging part, so I want to make sure I reach a stationary distribution
  #  in independent chains

  probelocs = input.data$loc
  probescores = input.data$data
  probe.homologies = input.data$homology


  # set up an approximation of the overall data
  # this is used in proposing peak jumps
  dens.spline.fit = curfit(c(probelocs, GENOME.LENGTH+1), c(abs(probescores), abs(probescores[1])), s=length(probelocs) - sqrt(2*length(probelocs)), periodic=TRUE) 
  #splyvals = predict.dierckx(spline.fit)[1:GENOME.LENGTH]
  dens.spline.norm =  integral.dierckx(dens.spline.fit)
  #input.data$splinedens = spline.predicted.dens
  
  
  if (is.null(pars.init)) {
    pars.init = list()
      

    ## step id
    pars.init['step'] = 0

    ## for the peaks:
    pars.init['n.peak'] = 10

    ### peak centers
    # we pick these randomly rather than using crude peak finding -- this slows convergence
    #  but also helps evaluate how well converged our model is if we do multiple chains
    #  from different starting points
    pars.init[['peak.mus']] = sample(probelocs, size=10, prob=abs(probescores))
    pars.init['nprobes'] = length(probelocs)
    #pars.init['minloc'] = min(probelocs)
    #pars.init['maxloc'] = max(probelocs)
    #pars.init['peak.mus'] = sample(minloc:maxloc, pars.init$n.peak, replace=FALSE)

    ### peak widths
    pars.init[['peak.sigmas']] = rep(75,pars.init$n.peak)

    ### peak heights
    pars.init[['peak.heights']] = rep(4,pars.init$n.peak)

    ## for the background

    ### mu_0
    pars.init['mu0'] = 0.647576

    ### sigma_probe
    pars.init['sigma.base'] = 0.2

    ### sigma.div.r
    pars.init['sigma.divergence.r'] = 0.1

    ### nu
    pars.init['nu'] = 20

    lik.prv = log.likelihood.main(input.data, pars.init)
    pars.init['loglik'] = lik.prv
  }

  # now do the sampling
  n.accept = 0
  current.pars = pars.init
  #all.step.pars = data.frame(n.peak=current.pars$n.peak, step=0, peak.mus=I(current.pars$peak.mus), nprobes=current.pars$nprobes, peak.sigmas = I(current.pars$peak.sigmas), peak.heights=I(current.pars$peak.heights), mu0 = current.pars$mu0, sigma.base = current.pars$sigma.base, sigma.divergence.r = current.pars$sigma.divergence.r, nu=current.pars$nu, loglik = current.pars$loglik)
  this.step.pars = make.pars.printable(current.pars)
  n.save = n.iter / (sparse+1)
  print(n.iter)
  print(n.iter/(sparse+1))
  all.step.pars = as.list(rep(0,n.iter/(sparse+1)))
  all.step.pars[[1]] =this.step.pars
  print(length(all.step.pars))
  #all.step.pars = list(current.pars)

  n.attempt.movetype = rep(0,NUM.MOVE.TYPES)
  n.accept.movetype = rep(0,NUM.MOVE.TYPES)

  saved.frame.count=2
  sparsecount=0

  png('init.png')
  plot(input.data$loc, input.data$data)
  points(input.data$loc, calc.mode.vals(input.data$loc, pars.init), col='red')
  points(current.pars$peak.mus, current.pars$peak.heights, col='black', bg='green', pch=21)
  dev.off()

  for (iter in 1:n.iter) {
    if (iter %% 50 == 0) {
      png(paste('newstep3_',iter,'.png',sep=""))
      print(paste("working on iter", iter))
      plot(input.data$loc, input.data$data)
      #points(input.data$loc, calc.mode.vals(input.data$loc, pars.init), col='red')
      lines(input.data$loc, calc.mode.vals(input.data$loc, current.pars), col='blue')
      if (current.pars$n.peak > 0) {
        print(current.pars$peak.mus)
        print(current.pars$peak.heights)
        points(current.pars$peak.mus, current.pars$peak.heights, col='black', bg='green', pch=21)
      }
      dev.off()
    }
    new.par.list = do.move(input.data,current.pars, dens.spline.fit, dens.spline.norm)
    current.pars = new.par.list$pars
    move.type = new.par.list$move.type
    n.attempt.movetype[move.type] = n.attempt.movetype[move.type] + 1
    if (new.par.list$accepted) {
      n.accept.movetype[move.type] = n.accept.movetype[move.type] + 1
      n.accept = n.accept + 1
    }

    current.pars$step = current.pars$step + 1
    #all.step.pars = rbind(all.step.pars, data.frame(make.pars.printable(current.pars), stringsAsFactors=FALSE))
    if (sparsecount >= sparse) {
      all.step.pars[[saved.frame.count]] = make.pars.printable(current.pars)
      sparsecount=0
      saved.frame.count = saved.frame.count+1
    } else {
      sparsecount = sparsecount+1
    }
  }

  print(paste("Accepted", n.accept, "out of", n.iter, "steps (", n.accept * 100 / n.iter, "% )"))
  print("Acceptance by move type:")
  print(n.accept.movetype / n.attempt.movetype)
  #return(all.step.pars)
  return(do.call("rbind", lapply(all.step.pars, as.data.frame)))

}

# ----- HELPER FUNCTIONS GO HERE ----------

# define priors for all of the parameters of interest
# to do this, we have to define functions that return the probability density of each parameter having a given value

# n.b. we work with log likelihoods in all cases

# note also that all of the "peak" parameters must have proper priors,
#  since they are subject to birth/death steps


# (alimf) Leaving functions in even if I do not need them per se. Unneeded functions will be marked with a comment
# on top specifying "Not necessary." Prior functions that will be used will be marked "in use"


#not necessary
d.n.peak.prior = function(npeaks) {
  # the number of peaks is poisson distributed with a fixed lambda
  return( dpois(npeaks, lambda=npeak.lambda, log=TRUE) )
}

#not necessary
d.peak.center.prior = function(all.data, this.loc) {
  # the centers have a uniform prior over the set of locations
  if ((this.loc < 0) || this.loc > GENOME.LENGTH) {
    return(NA)
  }

  return(log(1.0/GENOME.LENGTH))
}


#In use with minor mods
d.peak.sigma.prior = function(this.sigma) {
  # the sigmas have a truncated jeffreys prior -- uniform on 1/sigma
  #  we set broad upper and lower bounds to yield a proper prior

  if (this.sigma > max.peak.sigma) {
    return(NA)
  }

  if (this.sigma < min.peak.sigma) {
    return(NA)
  }

  return(log( (1.0/this.sigma) / peak.sigma.prior.int))
}


#To be modified: replace with prior for -deltaG
d.peak.height.prior = function(this.height) {
  # the magnitude of the peak heights  have a gamma prior, with the shape parameter chosen to prefer the peak to be substantially above background
  # we don't care for these purposes what the sign of the peak is -- the prior is symmetric about 0

  mag.height = (this.height)

  if (mag.height < peak.height.min) {
    return(NA)
  }

  return(dgamma(mag.height, height.prior.shape, scale=height.prior.scale, log=TRUE) - log.peak.height.prior.int)
}

#Using
d.energy.prior = function(this.energy) {

mag.energy = (this.energy)

#Log used because obscenely tiny
return(dgamma(mag.energy, energy.prior.shape, scale = energy.prior.scale, log = TRUE))

}


d.mismatch.prior = function(this.mismatch) {

mag.mismatch = (this.mismatch)

return(dgamma(mag.mismatch, mismatch.prior.shape, scale = mismatch.prior.scale, log = TRUE))

}

d.peak.sigma.prior = function(this.sigma) {

mag.sigma = (this.sigma)

return(dgamma(mag.sigma, peak.sigma.shape, scale = peak.sigma.scale, log = TRUE))

}

match = function(data, template) {

if(data = "A") {
return(template %in% c("A","R","W","M","V","H","D","N"))
}
else if(data = "T") {
return(template %in% c("T","Y","W","K","B","H","D","N"))
}
else if(data = "G") {
return(template %in% c("G","R","S","K","B","V","D","N"))
}
else if(data = "C") {
return(template %in% c("C","Y","S","M","B","V","H","N"))
}
else{
return(false)
}

}

  
d.mu0.prior = function(this.mu) {
  # the background mean itself has a N(0,1) prior. It should end up very close to 0,
  #  due to the regularization done before entering this function

  return(dnorm(this.mu, log=TRUE))
}

d.sigma.base.prior = function(this.sigma) {
  # the probe sigmas have an uninformative prior, but in practice they should end up
  #  quite close to 1

  return(log(1.0/this.sigma))
}

#Not being used
d.sigma.divergence.prior = function(this.sigma.div) {
  # the divergence slope gets a noninformative gamma prior -- this keeps it positive,
  #  but doesn't do much else

  return(dgamma(this.sigma.div, 1, scale=1, log=TRUE))
} 

d.nu.prior = function(this.nu) {
  # I need to come up with a good idea here
  # for now I just make it uniform
  return(log(1))
}

# combine all priors conveniently
d.prior = function(all.data, all.pars) {
  retval=0

  if (all.pars$n.peak > 0) {
    for (i in 1:all.pars$n.peak) {
      retval = retval + d.peak.center.prior(all.data, all.pars$peak.mus[i]) + d.peak.sigma.prior(all.pars$peak.sigmas[i]) + d.peak.height.prior(all.pars$peak.heights[i])
    }
  }


  retval = retval + d.n.peak.prior(all.pars$n.peak) + d.mu0.prior(all.pars$mu0) + d.sigma.base.prior(all.pars$sigma.base) + d.sigma.divergence.prior(all.pars$sigma.divergence.r) + d.nu.prior(all.pars$nu)

  return(retval)
}
  


# Now, the likelihood function
# here we plug in the data and current parameters and get back the density
# the likelihood is decomposed into several smaller functions to make it more manageable  
log.likelihood.main = function(all.data, all.pars) {
  # return the log-likelihood of the observed scores given the current parameters
  all.locs = all.data$loc
  all.scores = all.data$data
  all.homologies = all.data$homology

  # calculate some needed quantities using the helper functions below
  # these are the parameters of the distribution we want to model at each probe
  all.modes = calc.mode.vals(all.locs, all.pars)
  all.sds = calc.sd.vals(all.homologies, all.pars)

  # standardize the current values so I can use the standard t distribution
  standard.scores = (all.scores - all.modes)/all.sds

  #return(d.prior(all.data, all.pars) + sum(dt(standard.scores, df=all.pars$nu, log=TRUE)))
  lik.sum = 0
  #for (i in 1:length(all.modes)) {
    #lik.sum = lik.sum + dnorm(all.scores[i], mean=all.modes[i], sd=all.sds[i], log=TRUE)
    #lik.sum = lik.sum + dt.nstd(all.scores[i], means=all.modes[i], stds=all.sds[i], nu=all.pars$nu, log=TRUE)
  #  lik.sum = lik.sum + dt(standard.scores[i], df=all.pars$nu, log=TRUE) - log(all.sds[i])
    #lik.sum = lik.sum + dnorm(standard.scores[i], log=TRUE)
  #}

  lik.sum = sum( dt(standard.scores, df=all.pars$nu, log=TRUE) - log(all.sds) )
  
  return(d.prior(all.data, all.pars) + lik.sum)
}

calc.mode.vals = function(all.locs, all.pars) {
  # calculate the expected value of all probe-level scores based on the current distribution

  mode.vals = rep(all.pars$mu0, length(all.locs))

  peak.mus = all.pars$peak.mus
  peak.sigmas = all.pars$peak.sigmas
  peak.heights = all.pars$peak.heights

  if (all.pars$n.peak > 0) {
    for (i in 1:all.pars$n.peak) {
      this.center = peak.mus[i]
      this.width = peak.sigmas[i]
      this.height = peak.heights[i]
      
      dists.a = abs(this.center-all.locs)
      if (this.center > (GENOME.LENGTH/2)) {
        dists.b = abs((this.center-GENOME.LENGTH)-all.locs)
      } else {
        dists.b = abs((this.center+GENOME.LENGTH)-all.locs)
      }
      all.dists = pmin(dists.a, dists.b)

      normby = dnorm(0, mean=0, sd=this.width)
      test.dens = dnorm(all.dists, mean=0, sd=this.width)
      densities.thispeak = this.height * dnorm(all.dists, mean=0, sd=this.width) / normby
      
      mode.vals = mode.vals + densities.thispeak
    }
  }

  return(mode.vals)

}

calc.sd.vals = function(all.homologies, all.pars, probedists, probe.samestrands) {
  # given the current parameters and the homologies, calculate the 
  #  sd of the background distribution at each probe

  sd.vec = rep(all.pars$sigma.base, length(all.homologies))
  divergence = 1 - all.homologies
  sd.vec = sd.vec + all.pars$sigma.divergence.r * divergence

  return(sd.vec)

}


########################## DIFFERENT MONTE CARLO MOVE TYPES #########################

do.peak.move = function(all.data, current.pars) {
  # attempt a peak move

  old.pars = current.pars

  if (current.pars$n.peak == 0) {
    return(list(pars=old.pars, accepted=FALSE))
  }

  this.peak = sample.int(current.pars$n.peak, size=1)
  move.dist = rnorm(1,sd=peak.move.sigma)
  #1-2*(sample(2,size=1, replace=TRUE) -1) * rnbinom(1,size=1,prob=peak.move.p)
  new.loc = current.pars$peak.mus[this.peak] + move.dist

  # wrap the new peak position around the origin if needed -- we use a circular genome
  if (new.loc > GENOME.LENGTH) {
    new.loc = new.loc - GENOME.LENGTH
  }

  if (new.loc < 0) {
    new.loc = GENOME.LENGTH + new.loc
  }

  current.pars$peak.mus[this.peak] = new.loc
  current.pars$loglik = log.likelihood.main(all.data, current.pars)

  if (accept.test(current.pars$loglik, old.pars$loglik)) {
    return(list(pars = current.pars, accepted = TRUE))
  } else {
    return(list(pars=old.pars, accepted=FALSE))
  }
}

change.peak.width = function(all.data, current.pars) {
  # attempt  to change peak width

  old.pars = current.pars
  if (current.pars$n.peak == 0) {
    return(list(pars=old.pars, accepted=FALSE))
  }

  this.peak = sample.int(current.pars$n.peak, size=1)
  move.dist = rnorm(1,sd=peak.size.change.sigma)
  current.pars$peak.sigmas[this.peak] = (current.pars$peak.sigmas[this.peak]+move.dist)
  current.pars$loglik = log.likelihood.main(all.data, current.pars)

  if (accept.test(current.pars$loglik, old.pars$loglik)) {
    return(list(pars = current.pars, accepted = TRUE))
  } else {
    return(list(pars=old.pars, accepted=FALSE))
  }
}

find.new.trial.loc = function(target.quantile, dens.spline.fit, dens.spline.norm, tol=1e-5, maxiter=1000) {
  # find the location corresponding to a given quantile in the spline fit in
  # dens.spline.fit
  # we act iteratively, taking smaller and smaller steps until convergence

  cur.loc = 0
  cur.int = 0

  movesize = 0.5*GENOME.LENGTH
  n.iter = 0
  while (n.iter < maxiter) {
    newpos = find.new.trial.loc.iter(target.quantile, cur.loc, movesize, cur.int, dens.spline.fit, dens.spline.norm)
    cur.loc = newpos$loc
    cur.int = newpos$int
    if (abs(cur.int - target.quantile) < tol) {
      return(cur.loc)
    }

    n.iter = n.iter+1
    movesize = movesize/2.0
  }

  print("Warning: exceeded max iterations in find.new.trial.loc")
  print(paste("Current vals are", target.quantile, cur.int))
  return(cur.loc)

}

find.new.trial.loc.iter = function(target.quantile, curloc, stepsize, curint, dens.spline.fit, dens.spline.norm) {
  # helper function to find the given quantile in a spline density
  # try a move of the given step size and return the new location
  # the new location is either the current location (if we overshoot)
  # or the current location + step size if that is still not enough
  # also return the new value of the integral up to the current point

  newloc = curloc + stepsize
  newint = curint + integral.dierckx(dens.spline.fit, from=curloc, to=newloc) / dens.spline.norm
  if (newint > target.quantile) {
    return(list(loc=curloc, int=curint))
  } else {
    return(list(loc=newloc, int=newint))
  }
}
  
  

jump.peak = function(all.data, current.pars, dens.spline.fit, dens.spline.norm) {
  # move a peak to a random position dictated by the density of the data
  # this is an asymmetric move, so we adjust the acceptance ratio accordingly
  # this move is designed to aid rapid mixing
  # we randomly flip the sign as well to enable jumping to opposite-signed peaks
  # to simplify the calculations we assume that it must jump to the site of a probe

  old.pars = current.pars
  if (current.pars$n.peak == 0) {
    return(list(pars=old.pars, accepted=FALSE))
  }

  this.peak = sample.int(current.pars$n.peak, size=1)
  this.peak.loc = current.pars$peak.mus[this.peak]
  this.dens = predict.dierckx(dens.spline.fit, this.peak.loc) / dens.spline.norm

  randloc = runif(1)

  new.pos=find.new.trial.loc(randloc, dens.spline.fit, dens.spline.norm)

  new.dens = predict.dierckx(dens.spline.fit, new.pos) / dens.spline.norm
  dens.ratio = this.dens / new.dens
  current.pars$peak.mus[this.peak] = new.pos
  current.pars$peak.heights[this.peak] = (current.pars$peak.heights[this.peak]) * (1-2*(sample(2,size=1, replace=TRUE) -1))
  current.pars$loglik = log.likelihood.main(all.data, current.pars)

  if (accept.test(current.pars$loglik + log(dens.ratio), old.pars$loglik)) {
    return(list(pars = current.pars, accepted = TRUE))
  } else {
    return(list(pars=old.pars, accepted=FALSE))
  }
}

flip.peak = function(all.data, current.pars) {
  # try to change the sign on a peak
  old.pars = current.pars
  if (current.pars$n.peak == 0) {
    return(list(pars=old.pars, accepted=FALSE))
  }

  this.peak = sample.int(current.pars$n.peak, size=1)
  current.pars$peak.heights[this.peak] = -1 * current.pars$peak.heights[this.peak]
  current.pars$loglik = log.likelihood.main(all.data, current.pars)

  if (accept.test(current.pars$loglik, old.pars$loglik)) {
    return(list(pars = current.pars, accepted = TRUE))
  } else {
    return(list(pars=old.pars, accepted=FALSE))
  }
}

rescale.peak = function(all.data, current.pars) {
  # try to change the width of a peak
  old.pars = current.pars
  if (current.pars$n.peak == 0) {
    return(list(pars=old.pars, accepted=FALSE))
  }

  this.peak = sample.int(current.pars$n.peak, size=1)
  move.dist = rnorm(1,sd=peak.height.change.sigma)
  current.pars$peak.heights[this.peak] = (current.pars$peak.heights[this.peak]+move.dist)

  current.pars$loglik = log.likelihood.main(all.data, current.pars)

  if (accept.test(current.pars$loglik, old.pars$loglik)) {
    return(list(pars = current.pars, accepted = TRUE))
  } else {
    return(list(pars=old.pars, accepted=FALSE))
  }
}
    
change.background.mu = function(all.data, current.pars) {
  # try to change the background mean
  old.pars = current.pars

  move.dist = rnorm(1,sd=mu0.change.sigma)
  current.pars$mu0 = current.pars$mu0 + move.dist

  current.pars$loglik = log.likelihood.main(all.data, current.pars)

  if (accept.test(current.pars$loglik, old.pars$loglik)) {
    return(list(pars = current.pars, accepted = TRUE))
  } else {
    return(list(pars=old.pars, accepted=FALSE))
  }
}

change.sigma.base = function(all.data, current.pars) {
  # try to change the background mean
  old.pars = current.pars

  move.dist = rnorm(1,sd=sigma.base.change.sigma)
  current.pars$sigma.base = current.pars$sigma.base + move.dist

  current.pars$loglik = log.likelihood.main(all.data, current.pars)

  if (accept.test(current.pars$loglik, old.pars$loglik)) {
    return(list(pars = current.pars, accepted = TRUE))
  } else {
    return(list(pars=old.pars, accepted=FALSE))
  }
}
  
change.sigma.divergence.r = function(all.data, current.pars) {
  # try to change the dependence of sigma on homology
  old.pars = current.pars

  move.dist = rnorm(1,sd=sigma.divergence.change.sigma)
  current.pars$sigma.divergence.r = current.pars$sigma.divergence.r + move.dist

  current.pars$loglik = log.likelihood.main(all.data, current.pars)

  if (accept.test(current.pars$loglik, old.pars$loglik)) {
    return(list(pars = current.pars, accepted = TRUE))
  } else {
    return(list(pars=old.pars, accepted=FALSE))
  }
}
  
change.nu = function(all.data, current.pars) {
  # try to change the degrees of freedom in the background distribution
  old.pars = current.pars

  move.dist = rnorm(1,sd=nu.change.sigma)
  current.pars$nu = current.pars$nu + move.dist

  current.pars$loglik = log.likelihood.main(all.data, current.pars)

  if (accept.test(current.pars$loglik, old.pars$loglik)) {
    return(list(pars = current.pars, accepted = TRUE))
  } else {
    return(list(pars=old.pars, accepted=FALSE))
  }
}

#### MC moves involving adding/removing peaks
# these are particularly complicated -- we closely follow Green1995 in terms of notation
# and requirements to ensure appropriate behavior

peak.birth = function(all.data, current.pars) {
  # try to add a new peak at a random location
  # the new peak's parameters are chosen randomly from the
  #  corresponding prior distributions

  old.pars= current.pars

  # generate the parameters
  new.peak.mu = runif(1, 0, GENOME.LENGTH)
  new.peak.sigma = 1/runif(1, 1.0/max.peak.sigma, 1.0/min.peak.sigma)
  new.peak.sign = 1-2*(sample(2,size=1, replace=TRUE) -1)
  new.peak.height = rgamma(1,shape=height.prior.shape, scale=height.prior.scale)
  new.peak.height = new.peak.height * new.peak.sign

  # plug in the new peak
  if (current.pars$n.peak == 0) {
    current.pars$n.peak = 1
    current.pars[['peak.mus']] = c(new.peak.mu)
    current.pars[['peak.sigmas']] = c(new.peak.sigma)
    current.pars[['peak.heights']] = c(new.peak.height)
  } else {
    current.pars$n.peak = current.pars$n.peak + 1
    current.pars[['peak.mus']] = c(current.pars[['peak.mus']], new.peak.mu)
    current.pars[['peak.sigmas']] = c(current.pars[['peak.sigmas']], new.peak.sigma)
    current.pars[['peak.heights']] = c(current.pars[['peak.heights']], new.peak.height)
  }

  # calculate components of the acceptance ratio
  lik.portion = log.likelihood.main(all.data, current.pars)
  proposal.par.log.ratio = log( (1.0/GENOME.LENGTH) * ((1.0 / new.peak.sigma) / peak.sigma.prior.int) * (1.0/2.0) * dgamma(abs(new.peak.height), shape=height.prior.shape, scale=height.prior.scale))
  proposal.log.ratio = (-1 * proposal.par.log.ratio) + log(length(all.data$loc) / current.pars$n.peak)
  jacobian.log.ratio = 0
  current.pars$loglik = lik.portion

  #print(paste(lik.portion, proposal.par.log.ratio, proposal.log.ratio, jacobian.log.ratio, old.pars$loglik))

  if (accept.test( lik.portion+proposal.log.ratio+jacobian.log.ratio, old.pars$loglik)) {
    return(list(pars=current.pars, accepted=TRUE))
  } else {
    return(list(pars=old.pars, accepted=FALSE))
  }

}

peak.death = function(all.data, current.pars) {
  # try to delete a random peak

  if (current.pars$n.peak == 0) {
    return(list(pars=current.pars, accepted=FALSE))
  }
  

  # remove all parameters of the peak
  old.pars = current.pars
  del.peak.id = sample.int(current.pars$n.peak, size=1)
  current.pars$n.peak = current.pars$n.peak - 1

  if (del.peak.id == 1) {
    current.pars[['peak.mus']] = current.pars$peak.mus[2:old.pars$n.peak]
    current.pars[['peak.sigmas']] = current.pars$peak.sigmas[2:old.pars$n.peak]
    current.pars[['peak.heights']] = current.pars$peak.heights[2:old.pars$n.peak]
  } else if (del.peak.id == old.pars$n.peak) {
    current.pars[['peak.mus']] = current.pars$peak.mus[1:old.pars$n.peak-1]
    current.pars[['peak.sigmas']] = current.pars$peak.sigmas[1:old.pars$n.peak-1]
    current.pars[['peak.heights']] = current.pars$peak.heights[1:old.pars$n.peak-1]
  } else {
    current.pars[['peak.mus']] = c(current.pars$peak.mus[1:(del.peak.id-1)], current.pars$peak.mus[(del.peak.id+1):old.pars$n.peak])
    current.pars[['peak.sigmas']] = c(current.pars$peak.sigmas[1:(del.peak.id-1)], current.pars$peak.sigmas[(del.peak.id+1):old.pars$n.peak])
    current.pars[['peak.heights']] = c(current.pars$peak.heights[1:(del.peak.id-1)], current.pars$peak.heights[(del.peak.id+1):old.pars$n.peak])
  }

  # calculate components of the acceptance ratio
  lik.portion = log.likelihood.main(all.data, current.pars)
  proposal.par.log.ratio = log( (1.0/GENOME.LENGTH) * ((1.0 / old.pars[['peak.sigmas']][del.peak.id]) / peak.sigma.prior.int) * (1.0/2.0) * dgamma(abs(old.pars[['peak.heights']][del.peak.id]), shape=height.prior.shape, scale=height.prior.scale))
  proposal.log.ratio = proposal.par.log.ratio + log( old.pars$n.peak / length(all.data$loc))
  jacobian.log.ratio = 0
  current.pars$loglik = lik.portion

  if (accept.test( lik.portion+proposal.log.ratio+jacobian.log.ratio, old.pars$loglik)) {
    return(list(pars=current.pars, accepted=TRUE))
  } else {
    return(list(pars=old.pars, accepted=FALSE))
  }
}
  
peak.merge = function(all.data, current.pars) {
  # try to merge two randomly chosen peaks

  # randomly choose two peaks, delete one, and then give the other parameters
  #  that would immitate the combined peak if they were superimposed (or nearly so)

  old.pars = current.pars

  if (old.pars$n.peak < 2) {
    return(list(pars=old.pars, accepted=FALSE))
  }
  old.pars = current.pars
  del.peak.id = sample.int(current.pars$n.peak, size=1)
  merge.peak.id = sample.int(current.pars$n.peak, size=1)
  while (merge.peak.id == del.peak.id) {
    merge.peak.id = sample.int(current.pars$n.peak, size=1)
  }
  current.pars$n.peak = current.pars$n.peak - 1

  current.pars$peak.mus[merge.peak.id] = (current.pars$peak.mus[merge.peak.id] + current.pars$peak.mus[del.peak.id])/2
  current.pars$peak.heights[merge.peak.id] = current.pars$peak.heights[merge.peak.id] + current.pars$peak.heights[del.peak.id]
  current.pars$peak.sigmas[merge.peak.id] = sqrt(current.pars$peak.sigmas[merge.peak.id] * current.pars$peak.sigmas[del.peak.id])

  if (del.peak.id == 1) {
    current.pars[['peak.mus']] = current.pars$peak.mus[2:old.pars$n.peak]
    current.pars[['peak.sigmas']] = current.pars$peak.sigmas[2:old.pars$n.peak]
    current.pars[['peak.heights']] = current.pars$peak.heights[2:old.pars$n.peak]
  } else if (del.peak.id == old.pars$n.peak) {
    current.pars[['peak.mus']] = current.pars$peak.mus[1:old.pars$n.peak-1]
    current.pars[['peak.sigmas']] = current.pars$peak.sigmas[1:old.pars$n.peak-1]
    current.pars[['peak.heights']] = current.pars$peak.heights[1:old.pars$n.peak-1]
  } else {
    current.pars[['peak.mus']] = c(current.pars$peak.mus[1:(del.peak.id-1)], current.pars$peak.mus[(del.peak.id+1):old.pars$n.peak])
    current.pars[['peak.sigmas']] = c(current.pars$peak.sigmas[1:(del.peak.id-1)], current.pars$peak.sigmas[(del.peak.id+1):old.pars$n.peak])
    current.pars[['peak.heights']] = c(current.pars$peak.heights[1:(del.peak.id-1)], current.pars$peak.heights[(del.peak.id+1):old.pars$n.peak])
  }

  # calculate components of the acceptance ratio
  lik.portion = log.likelihood.main(all.data, current.pars)

  # see my notebook for the parameters going into U and the jacobian
  delta = genome.dist(old.pars$peak.mus[merge.peak.id], old.pars$peak.mus[del.peak.id]) / 2
  epsilon = (old.pars$peak.heights[merge.peak.id] - old.pars$peak.heights[del.peak.id])/2
  gamma = old.pars$peak.sigmas[merge.peak.id] / current.pars$peak.sigmas[merge.peak.id]
  #print("===")
  #print(delta)
  #print(epsilon)
  #print('---')
  u.proposal.log.ratio = dnorm(delta, sd=peaksplit.delta.sd, log=TRUE) + dnorm(epsilon, sd=peaksplit.epsilon.sd, log=TRUE) + dnorm(log(gamma), sd=peaksplit.gamma.sd, log=TRUE) 
  proposal.log.ratio = log( current.pars$n.peak ) + u.proposal.log.ratio
  jacobian.log.ratio = log(1/(old.pars$peak.sigmas[del.peak.id]) )
  current.pars$loglik = lik.portion

  #print(u.proposal.log.ratio)
  #print(proposal.log.ratio)
  #print(jacobian.log.ratio)
  #print(lik.portion)
  #print(old.pars$loglik)
  

  #print(lik.portion+proposal.log.ratio+jacobian.log.ratio)
  #print(old.pars$loglik)
  #print('xxx')
  
  if (accept.test( lik.portion+proposal.log.ratio+jacobian.log.ratio, old.pars$loglik)) {
    #print('succ')
    return(list(pars=current.pars, accepted=TRUE))
  } else {
    #print('fail')
    return(list(pars=old.pars, accepted=FALSE))
  }
}

peak.split = function(all.data, current.pars) {
  # try to split one randomly chosen peak into two
  # this is the opposite of peak.merge

  if (current.pars$n.peak == 0) {
    return( list(pars=current.pars, accepted=FALSE))
  }

  old.pars = current.pars
  split.peak.id = sample.int(current.pars$n.peak, size=1)
  new.peak.id = current.pars$n.peak + 1
  current.pars$n.peak = current.pars$n.peak + 1
  new.n.peak = current.pars$n.peak
  
  # generate variables for mapping this direction
  delta = rnorm(1, sd=peaksplit.delta.sd)
  epsilon = rnorm(1, sd=peaksplit.epsilon.sd)
  gamma = exp(rnorm(1, sd=peaksplit.gamma.sd))
  mu0 = current.pars[['peak.mus']][split.peak.id]
  h0 = current.pars[['peak.heights']][split.peak.id]
  sigma0 = current.pars[['peak.sigmas']][split.peak.id]

  # fill in the new parameters as needed
  current.pars[['peak.mus']][split.peak.id] = mu0 + delta
  current.pars[['peak.mus']][new.n.peak] = mu0 - delta
  current.pars[['peak.heights']][split.peak.id] = (h0 + epsilon)
  current.pars[['peak.heights']][new.n.peak] = (h0 - epsilon)
  current.pars[['peak.sigmas']][split.peak.id] = sigma0 * gamma
  current.pars[['peak.sigmas']][new.n.peak] = sigma0/gamma

  # calculate components of the acceptance ratio
  lik.portion = log.likelihood.main(all.data, current.pars)

  # see my notebook for the parameters going into U and the jacobian
  u.proposal.log.ratio = dnorm(delta, sd=peaksplit.delta.sd, log=TRUE) + dnorm(epsilon, sd=peaksplit.epsilon.sd, log=TRUE) + dnorm(log(gamma), sd=peaksplit.gamma.sd, log=TRUE) 
  proposal.log.ratio = -u.proposal.log.ratio - log( current.pars$n.peak ) 
  jacobian.log.ratio = log(  current.pars$peak.sigmas[new.n.peak] )
  current.pars$loglik = lik.portion

#  print("===")
#  print("---")
#  print("===")

  if (accept.test( lik.portion+proposal.log.ratio+jacobian.log.ratio, old.pars$loglik)) {
    return(list(pars=current.pars, accepted=TRUE))
  } else {
    return(list(pars=old.pars, accepted=FALSE))
  }

}


  

accept.test = function(lik.trial, lik.cur, current.pars) {
  # return true if the current step should be accepted, false otherwise

  lhr = lik.trial - lik.cur
  acceptance.r = log(runif(1))

  #print(paste("Old likelihood was", lik.cur))
  #print(paste("New likelihood was", lik.trial))
  #print(paste("log likelihood ratio is", lhr))
  #print(paste("Acceptance random number is", acceptance.r))


  if (is.na(lik.trial)) {
    return(FALSE)
  }

  if (lhr > acceptance.r) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}

sim.data = function(all.locs, all.homologies, all.pars) {
  # generate a simulated data set using the current parameters
  all.modes = calc.mode.vals(all.locs, all.pars)
  all.sds = calc.sd.vals(all.homologies, all.pars)

  sim.standard.scores = rt(length(all.locs), df=all.pars$nu)
  sim.scores = (sim.standard.scores * all.sds) + all.modes

  return(sim.scores)
}
  
dt.nstd = function( locs, means, stds, nu , log) {
  # return the T density function at a given set of points
  gamma.term = gamma( (nu + 1) / 2.0) / gamma( nu / 2.0) 
  var.term = sqrt( stds / (nu * pi) )
  loc.term = (1 + ( stds * (locs - means) * (locs - means) ) / nu) ^ (-1 * (nu+1) / 2.0)
  full.val = (gamma.term * var.term * loc.term)

  if (log) {
    return(log(full.val)) 
  } else {
    return(full.val)
  }
}

# dt' * (tau^-1/2) * 
# tau*(x^2 - 2xu + u^2) = x'
# x' /tau = x^2 - 2xu + u^2
# (x'/tau) - u^2 = x^2 - 2xu
# (x'/ tau) - u^2 = x(x-2u)

stan.z = function(x, mu, tau, nu) {
  return( sqrt( nu * (tau^(1-nu) - 1) + (tau^(-1*(nu+1)) * (x-mu)^2) ) )
}

do.move = function(all.data, current.pars, dens.spline.fit, dens.spline.norm) {
    # this function attempts a move of an appropriately chosen type,
    # tests for acceptance, and then returns the parameter vector for the next
    # set based on this test
    # return a list containing:
    #   pars: the new parameter vector
    #   move.type: the integer code for the move type
    #   accepted: flag of whether or not the move was accepted

    peak.jump.weight=0.2
    peak.birthdeath.weight=0.2
    peak.mergesplit.weight=0.1
    basic.move.weight = (1.0 - (peak.jump.weight + 2*peak.birthdeath.weight + 2*peak.mergesplit.weight))/8.0

    move.weights = c(basic.move.weight, basic.move.weight, basic.move.weight, basic.move.weight, basic.move.weight, basic.move.weight, basic.move.weight, basic.move.weight, peak.jump.weight, peak.birthdeath.weight, peak.birthdeath.weight, peak.mergesplit.weight, peak.mergesplit.weight)

    # first generate a move, and test its acceptance
    this.move = sample.int(NUM.MOVE.TYPES, size=1, prob=move.weights)
    #this.move = 1
    
    # each of these move types returns a list with the updated parameters and
    #  an "accepted" field for whether or not the move was accepted
   
    ## 1: move a peak center 
    if (this.move == 1) {
      new.pars = do.peak.move(all.data, current.pars)
    }

    ## 2: change a peak width
    if (this.move == 2) {
      new.pars = change.peak.width(all.data, current.pars)
    }

    ## 3: flip a peak
    if (this.move == 3) {
      new.pars = flip.peak(all.data, current.pars)
    }

    ## 4: rescale a peak
    if (this.move == 4) {
      new.pars = rescale.peak(all.data, current.pars)
    }

    ## 5: change background mean
    if (this.move == 5) {
      new.pars = change.background.mu(all.data, current.pars)
    }

    ## 6: change background sigma
    if (this.move == 6) {
      new.pars = change.sigma.base(all.data, current.pars)
    }

    ## 7: change divergence r parameter
    if (this.move == 7) {
      new.pars = change.sigma.divergence.r(all.data, current.pars)
    }

    ## 8: change nu
    if (this.move == 8) {
      new.pars = change.nu(all.data, current.pars)
    }

    ## 9: do a peak jump
    if (this.move == 9) {
      new.pars = jump.peak(all.data, current.pars, dens.spline.fit, dens.spline.norm)
    }

    ## 10: try to add a peak
    if (this.move == 10) {
      new.pars = peak.birth(all.data, current.pars)
    }

    ## 11: try to remove a peak
    if (this.move == 11) {
      new.pars = peak.death(all.data, current.pars)
    }

    ## 12: try a peak merge move
    if (this.move == 12) {
      new.pars = peak.merge(all.data, current.pars)
    }

    ## 13: try a peak split move
    if (this.move == 13) {
      new.pars= peak.split(all.data, current.pars)
    }

    #print(this.move)
    #print(new.pars)
    new.pars$move.type = this.move
    return(new.pars)
}

make.pars.printable=function(pars) {
  # create a version of pars that can be included in a data frame
  print.pars = pars
  print.pars[['peak.mus']] = paste(print.pars$peak.mus, collapse=';')
  print.pars[['peak.sigmas']] = paste(print.pars$peak.sigmas, collapse=';')
  print.pars[['peak.heights']] = paste(print.pars$peak.heights, collapse=';')

  return(print.pars)
}

parse.printed.pars = function(parstring) {
  # invert make.pars.printable to recover the pars list at a given step
  listpars = as.list(parstring)
  listpars[['peak.mus']] = as.numeric(strsplit(as.character(listpars[['peak.mus']]), ';')[[1]])
  listpars[['peak.sigmas']] = as.numeric(strsplit(as.character(listpars[['peak.sigmas']]), ';')[[1]])
  listpars[['peak.heights']] = as.numeric(strsplit(as.character(listpars[['peak.heights']]), ';')[[1]])
  return(listpars)
}

genome.dist = function( loc1, loc2 ) {
  # helper function to find distances in the (circular) bacterial genome
  # find the distance between loc1 and loc2 in a genome of length GENOME.LENGTH

  if (loc1 < loc2) {
    tmp=loc1
    loc1=loc2
    loc2=tmp
  }

  dist.basic=abs(loc1-loc2)
  dist.wrap = abs(loc1 - (loc2+GENOME.LENGTH))
  return(min(dist.basic, dist.wrap))
}

calc.posterior.density = function(all.data, par.list) {
  # calculate the posterior mean signal at all points on the trace from a given set of steps
  # par.list should be a slice of the data frame returned by the sampler, containing
  #  the production portion of the mc chain
  # all.data is the original data

  sum.density = rep(0, length(all.data$loc))

  nrows = nrow(par.list)
  for (i in 1:nrows) {
    if (par.list[i,]$n.peak > 0) {
      this.row = parse.printed.pars(par.list[i,])
      these.modes = calc.mode.vals(all.data$loc, this.row)
      sum.density = sum.density + these.modes
    }
  }

  return(sum.density / nrows)
}

find.posterior.peaks = function(probelocs, posterior.means, threshold=2) {
  # find all local minima/maxima in the average posterior signal
  # we do this by taking the discrete differences between posterior means at all probe locations, and finding places
  # where this trace (roughly the derivative) crosses zero
  # the magnitude of the posterior expected signal must be greater than a threshold for a given peak to be returned
  # probelocs is the set of probe locations, and posterior.means the expected values (based on the posterior distribution)
  # at each probe location
  # because we actually define peaks as existing between two probes, rather than within a single probe, we take the average of the 
  # neighboring probes when testing the threshold

  if (length(probelocs) != length(posterior.means) ) {
    stop("Lengths of probelocs and posterior.means must be the same!")
  }

  probelocs.leftof = c(probelocs[length(probelocs)], probelocs[1:(length(probelocs) - 1)])
  probe.middle.mean.locs = (probelocs + probelocs.leftof) / 2

  probevals.leftof = c(posterior.means[length(posterior.means)], posterior.means[1:(length(posterior.means) - 1)])
  probe.middle.mean.vals = (posterior.means + probevals.leftof) / 2


  
  posterior.diffs = diff(posterior.means)
  posterior.diffs.leftof = c(posterior.diffs[length(posterior.diffs)], posterior.diffs)
  posterior.diffs.rightof = c(posterior.diffs, posterior.diffs[1])

  peaklocs= ( sign(posterior.diffs.leftof) != sign(posterior.diffs.rightof) )
  peakloc.inds = (1:length(posterior.diffs))[peaklocs]
  peakloc.inds = peakloc.inds[is.finite(peakloc.inds)]
  peaklocs.good = c()
  for (peak in peakloc.inds) {
    if (abs(probe.middle.mean.vals[peak]) > threshold) {
      peaklocs.good = c(peaklocs.good, peak)
    }
  }

  print("Found peaks at:")
  print(peaklocs.good)
    
 } 
  

  
  
