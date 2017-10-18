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
GENOME.LENGTH = 4085
#GENOME.LENGTH = 1000
NUM.MOVE.TYPES = 7
READING.FRAME.LENGTH = 15
R = 1.98858775 #Uses kcal/(K*mol)
temperature = 298 #Kelvin

## parameters controlling move sizes -- tune for optimal acceptance
deltaG.move.sigma = 5
mismatch.move.sigma = 2
peak.size.change.sigma = 10
mu0.change.sigma=0.05
nu.change.sigma=0.5
sigma.base.change.sigma = 0.01



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




#In first implementation, tf.lamda will stand for static number of TFs.
#In later implementations, this will be shifted to a Poisson distribution with 
#tf.lamda as the Poisson variable
tf.lambda = 100

# v1 -- here we assume a fixed number of peaks and just have to fit their locations
# the background is independent student-t distributed, with a common mean for all elements
# and sd varying linearly with divergence (=1-homology)
# I *know* this is a bad model -- just need to implement the sampler first and then
# i can improve on it

mcmc.adam.peak.v1 = function(input.data, input.seq, n.iter, pars.init=NULL, sparse=0) {
  # do mcmc to find locations of one or more peaks
  # before being passed to this function, data should be regularized to have a mean 0 
  #  and sd of 1
  # IT MUST ALSO BE SORTED IN SPATIAL ORDER ON THE CHROMOSOME
  #  Failure to adhere to these requirements may make your research explode

  # input.data must contain ***
  # input.seq is the name of the fasta file containing the sequence that lines up with the waveform
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

  sequence = sanitize.seq(input.seq)
  
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
   
    ## number of tfs: note that this will be altered when we change over to reversible jump MCMC
    pars.init['tf.num'] = tf.lambda

    ### peak widths
    pars.init[['peak.sigmas']] = rep(75,pars.init$tf.num)

    ### tf binding energies
    pars.init[['deltaG']] = rep(4,pars.init$tf.num)

    ### mismatch energy
    pars.init[['mismatch']] = rep(.4,pars.init$tf.num)

    ### motif
    raw.bases = sample(bases.vector, READING.FRAME.LENGTH*pars.init$tf.num, replace = TRUE)
    for(i in 0:(pars.init$tf.num-1)){
	    pars.init[['motif']][i+1] = paste(raw.bases[(i*READING.FRAME.LENGTH+1):((i+1)*READING.FRAME.LENGTH)], collapse = "")
    }  

    ## for the background

    ### mu_0
    pars.init['mu0'] = 0.647576

    ### sigma_probe
    pars.init['sigma.base'] = 0.2

    ### nu
    pars.init['nu'] = 20

    lik.prv = log.likelihood.main(input.data, sequence, pars.init)
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
  points(input.data$loc, generate.waveform(input.data$loc, sequence, pars.init), col='red')
  #points(current.pars$peak.mus, current.pars$peak.heights, col='black', bg='green', pch=21)
  dev.off()

  for (iter in 1:n.iter) {
    if (iter %% 50 == 0) {
      png(paste('newstep3_',iter,'.png',sep=""))
      print(paste("working on iter", iter))
      plot(input.data$loc, input.data$data)
      #points(input.data$loc, calc.mode.vals(input.data$loc, pars.init), col='red')
      lines(input.data$loc, generate.waveform(input.data$loc, sequence, current.pars), col='blue')
      if (current.pars$tf.num > 0) {
        print(current.pars$deltaG)
        print(current.pars$mismatch)
	print(current.pars$motif)
        #points(current.pars$peak.mus, current.pars$peak.heights, col='black', bg='green', pch=21)
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
#this is the function to sanitize and define the sequence while lining it up to the data
#Note that if a FASTA file contains multiple sequences, this function will ignore the the sequences that come after
#the first one
#Note also that this requires file extension fasta
sanitize.seq = function(input.seq) {

	split.file.name = unlist(strsplit(input.seq,'\\.'))
	if(!(identical(split.file.name[2],"fasta"))){
		stop('This file is not a fasta file')
	}

	con = file(input.seq, open="r")
	raw.read = readLines(con, n = 2)

	raw.sequence = raw.read[2]

	raw.sequence.vector =  unlist(strsplit(raw.sequence,''))

	is.acceptable.bp = raw.sequence.vector %in% c('A','T','G','C')

	if(!all(is.acceptable.bp)){
		stop('Sequences must only have base pairs A, T,G,C. No other base pairs are allowed.')
	}

	#This circularizes the genome. We will edit as neccesary
	sequence.vector = c(raw.sequence.vector, raw.sequence.vector[1:(READING.FRAME.LENGTH-1)])

	return(sequence.vector)

}


#This is the function that generates the peak waveform we need to work with

generate.waveform = function(this.locs, this.sequence, all.pars) {
	
	#frame.matrix generates the list of sequence frame	
	frame.matrix = matrix(this.sequence[1:(length(this.sequence)-(READING.FRAME.LENGTH-1))], ncol = 1)
	for(i in 2:READING.FRAME.LENGTH){
		
		frame.matrix = cbind(frame.matrix,this.sequence[i:(length(this.sequence)-(READING.FRAME.LENGTH-1)+i-1)])
	
	}
	
	mag.pars = all.pars
	
	#Unlists the sequence, splits the individual strings in the resulting vector.
	#This makes a list. This new list is unlisted, which makes a vector of characters
	#We need this vector to be a matrix with 15 columns, and each row contains the current sequnce for the TF.
	#But, r fills matrices down columns, so we make the matrix with reading frames going down, then transpose
	sequence.motifs = t(matrix(unlist(strsplit(unlist(mag.pars$motif),"")), nrow = READING.FRAME.LENGTH))
	

        mismatch = unlist(all.pars$mismatch)

        widths = unlist(all.pars$widths)
	
	waveform = numeric(length(this.locs))


	#Iterates through every readin frame and determines the contribution of that reading frame to the overall wave
	for(i in 1:(length(this.sequence)-READING.FRAME.LENGTH)){
		
		#energy is redefined every iteration because it's subtracted off thanks to mismatches
		energy = unlist(all.pars$deltaG)
		
		current.frame = frame.matrix[i,]
			
		is.mismatch = !(match(current.frame, sequence.motifs))
		
		#Exploits TRUE == 1 to get number of mismatches for each frame comparing to sequence motif
		number.mismatches = rowSums(is.mismatch)

		energy = energy - number.mismatches*mismatch
		
		peak.heights = exp(energy/R*temperature)
		
		#This eliminates any attempt at peaks that don't really contribute to the waveform
		peak.heights[peak.heights <= 1] = 0

		peak.coeffs = peakheights*sqrt(2*pi*widths)
		
		location.means = (1:(length(this.sequence)-(READING.FRAME.LENGTH-1))+READING.FRAME.LENGTH:length(this.sequence))/2
		
		for(j in 1:length(this.locs)){
		
			waveform[j] = waveform[j]+sum(peak.coeffs*dnorm[(i+8)+genome.dist(i+8,this.locs[j]), mean = (i+8), sd = sqrt(widths)])

		}
	}	

	return(allpars['mu0']+waveform)
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

match = function(sequence.frames, motif) {

	isA = (sequence.frames == "A")
	matchA = motif %in% c("A","R","W","M","V","H","D","N")
	Amatches = (isA & matchA)

	isT = (sequence.frames == "T") 
	matchT = motif %in% c("T","Y","W","K","B","H","D","N")
	Tmatches = (isT & matchT)

	isG = (sequence.frames == "G")
	matchG = motif %in% c("G","R","S","K","B","V","D","N")
	Gmatches = (isG & matchG)
	

	isC = (sequence.frames == "C")
	matchC = motif %in% c("C","Y","S","M","B","V","H","N")
	Cmatches = (isC & matchC)

	return (Amatches+Tmatches+Gmatches+Cmatches)	

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

d.nu.prior = function(this.nu) {
  # I need to come up with a good idea here
  # for now I just make it uniform
  return(log(1))
}

# combine all priors conveniently
d.prior = function(all.data, all.pars) {
  retval=0



  retval = retval + d.energy.prior(all.pars$deltaG) + d.mu0.prior(all.pars$mu0) + d.sigma.base.prior(all.pars$sigma.base) + d.peak.sigma.prior(all.pars$peak.sigma) + d.peak.mismatch + d.nu.prior(all.pars$nu)

  return(retval)
}
  


# Now, the likelihood function
# here we plug in the data and current parameters and get back the density
# the likelihood is decomposed into several smaller functions to make it more manageable  
log.likelihood.main = function(all.data, sequence, all.pars) {
  # return the log-likelihood of the observed scores given the current parameters
  all.locs = all.data$loc
  all.scores = all.data$data

  # calculate some needed quantities using the helper functions below
  # these are the parameters of the distribution we want to model at each probe
  all.modes = generate.waveform(all.locs, sequence, all.pars)

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





########################## DIFFERENT MONTE CARLO MOVE TYPES #########################


change.peak.width = function(all.data, sequence, current.pars) {
  # attempt  to change peak width

  old.pars = current.pars

  this.peak = sample.int(current.pars$n.peak, size=1,replace = TRUE)
  move.dist = rnorm(1,sd=peak.size.change.sigma)
  current.pars$peak.sigmas[this.peak] = (current.pars$peak.sigmas[this.peak]+move.dist)
  current.pars$loglik = log.likelihood.main(all.data, sequence, current.pars)

  if (accept.test(current.pars$loglik, old.pars$loglik)) {
    return(list(pars = current.pars, accepted = TRUE))
  } else {
    return(list(pars=old.pars, accepted=FALSE))
  }
}

do.deltaG.move = function(all.data, sequence, current.pars) {

	old.pars = current.pars;
	
	this.deltaG = sample.int(current.pars$num.tfs, size = 1,replace = TRUE)
	

	move.dist = rnorm(1,sd = deltaG.move.sigma)

	new.deltaG = current.pars$deltaG[this.deltaG] + move.dist

	current.pars$deltaG[this.deltaG] = new.deltaG

	current.pars$loglik = log.likelihood.main(all.data, sequence, current.pars)

	if (accept.test(current.pars$loglik, old.pars$loglik)) {
    		return(list(pars = current.pars, accepted = TRUE))
  	} else {
    		return(list(pars=old.pars, accepted=FALSE))
  	}
	
}

do.mismatch.move = function(all.data, sequence, current.pars) {

        old.pars = current.pars;

        this.mismatch = sample.int(current.pars$num.tfs, size = 1, replace = TRUE)


        move.dist = rnorm(1,sd = mismatch.move.sigma)

        new.mismatch = current.pars$deltaG[this.mismatch] + move.dist

        current.pars$mismatch[this.mistmatch] = new.mismatch

        current.pars$loglik = log.likelihood.main(all.data, sequence, current.pars)

        if (accept.test(current.pars$loglik, old.pars$loglik)) {
                return(list(pars = current.pars, accepted = TRUE))
        } else {
                return(list(pars=old.pars, accepted=FALSE))
        }

}

do.motif.move = function(all.data, sequence, current.pars) {

        old.pars = current.pars;

        this.motif = sample.int(current.pars$num.tfs, size = 1, replace = TRUE)
	
	motif.vector = strsplit(old.pars$motif[this.motif])

	base.change = sample.int(READING.FRAME.LENGTH, size = 1, replace = TRUE)
	
	motif.vector[base.change] = sample(bases.vector, size = 1, replace = TRUE)

        current.pars$motif[this.motif] = paste(motif.vector, collapse = "")

        current.pars$loglik = log.likelihood.main(all.data, sequence, current.pars)

        if (accept.test(current.pars$loglik, old.pars$loglik)) {
                return(list(pars = current.pars, accepted = TRUE))
        } else {
                return(list(pars=old.pars, accepted=FALSE))
        }

}

    
change.background.mu = function(all.data, sequence, current.pars) {
  # try to change the background mean
  old.pars = current.pars

  move.dist = rnorm(1,sd=mu0.change.sigma)
  current.pars$mu0 = current.pars$mu0 + move.dist

  current.pars$loglik = log.likelihood.main(all.data, sequence, current.pars)

  if (accept.test(current.pars$loglik, old.pars$loglik)) {
    return(list(pars = current.pars, accepted = TRUE))
  } else {
    return(list(pars=old.pars, accepted=FALSE))
  }
}

change.sigma.base = function(all.data, sequence, current.pars) {
  # try to change the background mean
  old.pars = current.pars

  move.dist = rnorm(1,sd=sigma.base.change.sigma)
  current.pars$sigma.base = current.pars$sigma.base + move.dist

  current.pars$loglik = log.likelihood.main(all.data, sequence, current.pars)

  if (accept.test(current.pars$loglik, old.pars$loglik)) {
    return(list(pars = current.pars, accepted = TRUE))
  } else {
    return(list(pars=old.pars, accepted=FALSE))
  }
}
  
  
change.nu = function(all.data, sequence, current.pars) {
  # try to change the degrees of freedom in the background distribution
  old.pars = current.pars

  move.dist = rnorm(1,sd=nu.change.sigma)
  current.pars$nu = current.pars$nu + move.dist

  current.pars$loglik = log.likelihood.main(all.data, sequence, current.pars)

  if (accept.test(current.pars$loglik, old.pars$loglik)) {
    return(list(pars = current.pars, accepted = TRUE))
  } else {
    return(list(pars=old.pars, accepted=FALSE))
  }
}

#### MC moves involving adding/removing peaks
# these are particularly complicated -- we closely follow Green1995 in terms of notation
# and requirements to ensure appropriate behavior

  

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

do.move = function(all.data, sequence, current.pars, dens.spline.fit, dens.spline.norm) {
    # this function attempts a move of an appropriately chosen type,
    # tests for acceptance, and then returns the parameter vector for the next
    # set based on this test
    # return a list containing:
    #   pars: the new parameter vector
    #   move.type: the integer code for the move type
    #   accepted: flag of whether or not the move was accepted

	
    basic.move.weight = 1/(NUM.MOVE.TYPES)
	
    move.weights = rep(base.move.weight,NUM.MOVE.TYPES)

    # first generate a move, and test its acceptance
    this.move = sample.int(NUM.MOVE.TYPES, size=1, replace = TRUE, prob=move.weights)
    #this.move = 1
    
    # each of these move types returns a list with the updated parameters and
    #  an "accepted" field for whether or not the move was accepted
   
    ## 1: change a TF energy 
    if (this.move == 1) {
      new.pars = do.deltaG.move(all.data, sequence, current.pars)
    }

    ## 2: change a peak width
    if (this.move == 2) {
      new.pars = change.peak.width(all.data, sequence, current.pars)
    }

    ## 3: flip a peak
    if (this.move == 3) {
      new.pars = do.mismatch.move(all.data, sequence, current.pars)
    }

    ## 4: rescale a peak
    if (this.move == 4) {
      new.pars = do.motif.move(all.data, sequence, current.pars)
    }

    ## 5: change background mean
    if (this.move == 5) {
      new.pars = change.background.mu(all.data, sequence, current.pars)
    }

    ## 6: change background sigma
    if (this.move == 6) {
      new.pars = change.sigma.base(all.data, sequence, current.pars)
    }

    ## 7: change nu
    if (this.move == 7) {
      new.pars = change.nu(all.data, sequence, current.pars)
    }


    #print(this.move)
    #print(new.pars)
    new.pars$move.type = this.move
    return(new.pars)
}

make.pars.printable=function(pars) {
  # create a version of pars that can be included in a data frame
  print.pars = pars
  print.pars[['deltaG']] = paste(print.pars$deltaG, collapse=';')
  print.pars[['peak.sigmas']] = paste(print.pars$peak.sigmas, collapse=';')
  print.pars[['mismatch']] = paste(print.pars$mismatch, collapse=';')

  return(print.pars)
}

parse.printed.pars = function(parstring) {
  # invert make.pars.printable to recover the pars list at a given step
  listpars = as.list(parstring)
  listpars[['deltaG']] = as.numeric(strsplit(as.character(listpars[['deltaG']]), ';')[[1]])
  listpars[['peak.sigmas']] = as.numeric(strsplit(as.character(listpars[['peak.sigmas']]), ';')[[1]])
  listpars[['mismatch']] = as.numeric(strsplit(as.character(listpars[['mismatch']]), ';')[[1]])
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
  

  
  
