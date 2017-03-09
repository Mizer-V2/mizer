# Methods used for projecting for the size based modelling package

# Copyright 2012 Finlay Scott and Julia Blanchard. Distributed under the GPL 2 or later
# Maintainer: Finlay Scott, CEFAS

# Calculate the amount of food exposed to each predator by predator size

#' getPhiPrey method for the size based model
#' 
#' Calculates the amount \eqn{E_{a,i}(w)} of food exposed to each predator by
#' predator size. This method is used by the \code{\link{project}} method for
#' performing simulations.
#' @param object An \linkS4class{MizerParams} object
#' @param n A matrix of species abundances (species x size)
#' @param n_pp A vector of the background abundance by size
#' @param tau The acitivity level
#' 
#' @return A two dimensional array (predator species x predator size)
#' @seealso \code{\link{project}}
#' @export
#' @examples
#' \dontrun{
#' data(NS_species_params_gears)
#' data(inter)
#' params <- MizerParams(NS_species_params_gears, inter)
#' # With constant fishing effort for all gears for 20 time steps
#' sim <- project(params, t_max = 20, effort = 0.5)
#' n <- sim@n[21,,]
#' n_pp <- sim@@n_pp[21,]
#' getPhiPrey(params,n,n_pp)
#' }

getPhiPrey <- function(object, n, n_pp, tau){
  #        cat("In getPhiPrey\n")
  # Check n dims
  if(dim(n)[1] != dim(object@interaction)[1])
    stop("n does not have the right number of species (first dimension)")
  if(dim(n)[2] != length(object@w))
    stop("n does not have the right number of size groups (second dimension)")
  if(length(n_pp) != length(object@w_full))
    stop("n_pp does not have the right number of size groups")
  # n_eff_prey is the total prey abundance by size exposed to each predator (prey
  # not broken into species - here we are just working out how much a predator
  # eats - not which species are being eaten - that is in the mortality calculation
  n_eff_prey <- sweep((object@interaction %*% n) * tau, 2, object@w * object@dw, "*") 
  # Quick reference to just the fish part of the size spectrum
  idx_sp <- (length(object@w_full) - length(object@w) + 1):length(object@w_full)
  # pred_kernel is predator x predator size x prey size
  # So multiply 3rd dimension of pred_kernel by the prey abundance
  # Then sum over 3rd dimension to get total eaten by each predator by predator size
  phi_prey_species <- rowSums(sweep(object@pred_kernel[,,idx_sp,drop=FALSE],c(1,3),n_eff_prey,"*"),dims=2)
  # Eating the background
  phi_prey_background <- rowSums(sweep(object@pred_kernel,3,object@dw_full*object@w_full*n_pp,"*"),dims=2)
  return(phi_prey_species+phi_prey_background)
}


# Feeding level
# The amount of food consumed by a predator, by each predator size

#' getFeedingLevel method for the size based model
#' 
#' Calculates the amount of food \eqn{f_i(w)} consumed by a predator by predator
#' size based on food availability, search volume and maximum intake. This
#' method is used by the \code{\link{project}} method for performing
#' simulations.
#' @param object A \code{MizerParams} or \code{MizerSim} object
#' @param n A matrix of species abundance (species x size). Only used if 
#'   \code{object} argument is of type \code{MizerParams}.
#' @param n_pp A vector of the background abundance by size. Only used if 
#'   \code{object} argument is of type \code{MizerParams}.
#' @param phi_prey The PhiPrey matrix (optional) of dimension no. species x no. 
#'   size bins. If not passed in, it is calculated internally using the 
#'   \code{\link{getPhiPrey}} method. Only used if \code{object} argument is of type 
#'   \code{MizerParams}.
#' @param time_range Subset the returned fishing mortalities by time. The time 
#'   range is either a vector of values, a vector of min and max time, or a 
#'   single value. Default is the whole time range. Only used if the 
#'   \code{object} argument is of type \code{MizerSim}.
#' @param drop should extra dimensions of length 1 in the output be dropped, 
#'   simplifying the output. Defaults to TRUE
#' @param tau The acitivity level
#'   
#' @note If a \code{MizerParams} object is passed in, the method returns a two 
#' dimensional array (predator species x predator size) based on the abundances 
#' also passed in.
#' 
#' If a \code{MizerSim} object is passed in, the method returns a three
#' dimensional array (time step x predator species x predator size) with the
#' feeding level calculated at every time step in the simulation.
#' @seealso \code{\link{getPhiPrey}}
#' @export
#' @examples
#' \dontrun{
#' data(NS_species_params_gears)
#' data(inter)
#' params <- MizerParams(NS_species_params_gears, inter)
#' # With constant fishing effort for all gears for 20 time steps
#' sim <- project(params, t_max = 20, effort = 0.5)
#' # Get the feeding level at one time step
#' n <- sim@@n[21,,]
#' n_pp <- sim@@n_pp[21,]
#' fl <- getFeedingLevel(params,n,n_pp)
#' # Get the feeding level at all saved time steps
#' fl <- getFeedingLevel(sim)
#' # Get the feeding level for time 15 - 20
#' fl <- getFeedingLevel(sim, time_range = c(15,20))
#' }

getFeedingLevel <- function(object, n, n_pp, phi_prey, tau){
  if (missing(phi_prey)) phi_prey <- getPhiPrey(object, n=n, n_pp=n_pp, tau=tau)
  # Check dims of phi_prey
  if (!all(dim(phi_prey) == c(nrow(object@species_params),length(object@w)))){
    stop("phi_prey argument must have dimensions: no. species (",nrow(object@species_params),") x no. size bins (",length(object@w),")")
  }
  # encountered food = available food * search volume
 
  encount <- tau * object@search_vol * phi_prey
  # calculate feeding level
  f <- encount/(encount + object@intake_max)
  return(f)
}

# Predation rate
# Soundtrack: Nick Drake - Pink Moon

#' getPredRate method for the size based model
#' 
#' Calculates the predation rate of each predator species at size on prey size. 
#' In formulas \deqn{\phi_i(w_p/w) (1-f_i(w)) \gamma_i w^q N_i(w) dw}
#' This method is used by the \code{\link{project}} method for performing
#' simulations. In the simulations, it is combined with the interaction matrix
#' (see \code{\link{MizerParams}}) to calculate the realised predation mortality
#' (see \code{\link{getM2}}).
#' @param object A \code{MizerParams} object.
#' @param n A matrix of species abundance (species x size).
#' @param n_pp A vector of the background abundance by size.
#' @param feeding_level The current feeding level (optional). A matrix of size
#'   no. species x no. size bins. If not supplied, is calculated internally
#'   using the \code{getFeedingLevel()} method.
#' @param tau The acitivity level
#'   
#' @return A three dimensional array (predator species x predator size x prey size), 
#'   where the predator size runs over the community size range only and prey size
#'   runs over community plus background spectrum.
#' @export
#' @seealso \code{\link{project}}, \code{\link{getM2}}, \code{\link{getFeedingLevel}} and \code{\link{MizerParams}}
#' @examples
#' \dontrun{
#' data(NS_species_params_gears)
#' data(inter)
#' params <- MizerParams(NS_species_params_gears, inter)
#' # With constant fishing effort for all gears for 20 time steps
#' sim <- project(params, t_max = 20, effort = 0.5)
#' # Get the feeding level at one time step
#' n <- sim@@n[21,,]
#' n_pp <- sim@@n_pp[21,]
#' getPredRate(params,n,n_pp)
#' }
getPredRate <- function(object, n, n_pp, feeding_level, tau){
  if(missing(feeding_level)) feeding_level <- getFeedingLevel(object, n=n, n_pp=n_pp, tau=tau)
  if (!all(dim(feeding_level) == c(nrow(object@species_params),length(object@w)))){
    stop("feeding_level argument must have dimensions: no. species (",nrow(object@species_params),") x no. size bins (",length(object@w),")")
  }
  n_total_in_size_bins <- sweep(n, 2, object@dw, '*') # N_i(w)dw
  pred_rate <- sweep(object@pred_kernel,c(1,2),(1-feeding_level)*tau*object@search_vol*n_total_in_size_bins,"*")
  return(pred_rate)
}


# getM2
# This uses the predation rate which is also used in M2background
# Too much overlap? Inefficient? Same thing is calculated twice

#' getM2 method for the size based model
#'
#' Calculates the total predation mortality \eqn{\mu_{p,i}(w_p)} on each prey
#' species by prey size. This method is used by the \code{\link{project}} method
#' for performing simulations.
#' @param object A \code{MizerParams} or \code{MizerSim} object.
#' @param n A matrix of species abundance (species x size). Only used if
#'   \code{object} argument is of type \code{MizerParams}.
#' @param n_pp A vector of the background abundance by size. Only used if
#'   \code{object} argument is of type \code{MizerParams}.
#' @param pred_rate An array of predation rates of dimension no. sp x no.
#'   community size bins x no. of size bins in whole spectra (i.e. community +
#'   background, the w_full slot). The array is optional. If it is not provided
#'   it is calculated by the \code{getPredRate()} method.
#' @param time_range Subset the returned fishing mortalities by time. The time
#'   range is either a vector of values, a vector of min and max time, or a
#'   single value. Default is the whole time range. Only used if the
#'   \code{object} argument is of type \code{MizerSim}.
#' @param drop Only used when object is of type \code{MizerSim}. Should
#'   dimensions of length 1 in the output be dropped, simplifying the output.
#'   Defaults to TRUE
#' @param tau The acitivity level
#'
#' @return
#'   If a \code{MizerParams} object is passed in, the method returns a two
#'   dimensional array (prey species x prey size) based on the abundances also
#'   passed in. If a \code{MizerSim} object is passed in, the method returns a
#'   three dimensional array (time step x prey species x prey size) with the
#'   predation mortality calculated at every time step in the simulation.
#' @seealso \code{\link{getPredRate}} and \code{\link{project}}.
#' @export
#' @examples
#' \dontrun{
#' data(NS_species_params_gears)
#' data(inter)
#' params <- MizerParams(NS_species_params_gears, inter)
#' # With constant fishing effort for all gears for 20 time steps
#' sim <- project(params, t_max = 20, effort = 0.5)
#' # Get M2 at one time step
#' n <- sim@@n[21,,]
#' n_pp <- sim@@n_pp[21,]
#' getM2(params,n,n_pp)
#' # Get M2 at all saved time steps
#' getM2(sim)
#' # Get M2 over the time 15 - 20
#' getM2(sim, time_range = c(15,20))
#' }
getM2 <- function(object, n, n_pp, pred_rate, tau){
  if(missing(pred_rate)) pred_rate <- getPredRate(object,n=n,n_pp=n_pp, tau=tau)
  if ((!all(dim(pred_rate) == c(nrow(object@species_params),length(object@w),length(object@w_full)))) | (length(dim(pred_rate))!=3)){
    stop("pred_rate argument must have 3 dimensions: no. species (",nrow(object@species_params),") x no. size bins (",length(object@w),") x no. size bins in community + background (",length(object@w_full),")")
  }
  # get the element numbers that are just species
  idx_sp <- (length(object@w_full) - length(object@w) + 1):length(object@w_full)
  # Interaction is predator x prey so need to transpose so it is prey x pred
  # Sum pred_kernel over predator sizes to give total predation rate of
  # each predator on each prey size
  m2 <- t(object@interaction) %*% colSums(aperm(pred_rate, c(2,1,3)),dims=1)[,idx_sp]
  return(m2)
}


#' getM2Background method for the size based model
#'
#' Calculates the predation mortality \eqn{\mu_p(w)} on the background spectrum
#' by prey size. Used by the \code{project} method for running size based
#' simulations.
#' @param object A \code{MizerParams} object.
#' @param n A matrix of species abundance (species x size).
#' @param n_pp A vector of the background abundance by size.
#' @param pred_rate An array of predation rates of dimension no. sp x no.
#'   community size bins x no. of size bins in whole spectra (i.e. community +
#'   background, the w_full slot). The array is optional. If it is not provided
#'   it is calculated by the \code{getPredRate()} method.
#' @param tau The acitivity level
#' 
#' @return A vector of predation mortalities by background prey size.
#' @seealso \code{\link{project}} and \code{\link{getM2}}.
#' @export
#' @examples
#' \dontrun{
#' data(NS_species_params_gears)
#' data(inter)
#' params <- MizerParams(NS_species_params_gears, inter)
#' # With constant fishing effort for all gears for 20 time steps
#' sim <- project(params, t_max = 20, effort = 0.5)
#' # Get M2 of the background spectrum at one time step
#' n <- sim@@n[21,,]
#' n_pp <- sim@@n_pp[21,]
#' getM2Background(params,n,n_pp)
#' }

getM2Background <- function(object, n, n_pp, pred_rate, tau){
  if(missing(pred_rate)) pred_rate <- getPredRate(object,n=n,n_pp=n_pp,tau=tau)
  if ((!all(dim(pred_rate) == c(nrow(object@species_params),length(object@w),length(object@w_full)))) | (length(dim(pred_rate))!=3)){
    stop("pred_rate argument must have 3 dimensions: no. species (",nrow(object@species_params),") x no. size bins (",length(object@w),") x no. size bins in community + background (",length(object@w_full),")")
  }
  
  M2background <- colSums(pred_rate,dims=2)
  return(M2background)
}


# getFMortGear
#' Get the fishing mortality by time, gear, species and size
#'
#' Calculates the fishing mortality by gear, species and size at each time step
#' in the \code{effort} argument. Used by the \code{project} method to perform
#' simulations.
#' 
#' @param object A \code{MizerParams} object or a \code{MizerSim} object.
#' @param effort The effort of each fishing gear. Only needed if the object
#'   argument is of class \code{MizerParams}. See notes below.
#' @param time_range Subset the returned fishing mortalities by time. The time
#'   range is either a vector of values, a vector of min and max time, or a
#'   single value. Default is the whole time range. Only used if the
#'   \code{object} argument is of type \code{MizerSim}.
#'   
#' @return An array. If the effort argument has a time dimension, or a
#'   \code{MizerSim} is passed in, the output array has four dimensions (time x
#'   gear x species x size). If the effort argument does not have a time
#'   dimension (i.e. it is a vector or a single numeric), the output array has
#'   three dimensions (gear x species x size).
#' @note Here: fishing mortality = catchability x selectivity x effort.
#' 
#' The \code{effort} argument is only used if a \code{MizerParams} object is
#' passed in. The \code{effort} argument can be a two dimensional array (time x
#' gear), a vector of length equal to the number of gears (each gear has a
#' different effort that is constant in time), or a single numeric value (each
#' gear has the same effort that is constant in time). The order of gears in the
#' \code{effort} argument must be the same the same as in the \code{MizerParams}
#' object.
#' 
#' If the object argument is of class \code{MizerSim} then the effort slot of
#' the \code{MizerSim} object is used and the \code{effort} argument is not
#' used.
#' @export
#' @examples
#' \dontrun{
#' data(NS_species_params_gears)
#' data(inter)
#' params <- MizerParams(NS_species_params_gears, inter)
#' # Get the fishing mortality when effort is constant
#' # for all gears and time:
#' getFMortGear(params, effort = 1)
#' # Get the fishing mortality when effort is different
#' # between the four gears but constant in time:
#' getFMortGear(params, effort = c(0.5,1,1.5,0.75))
#' # Get the fishing mortality when effort is different
#' # between the four gears and changes with time:
#' effort <- array(NA, dim = c(20,4))
#' effort[,1] <- seq(from=0, to = 1, length=20)
#' effort[,2] <- seq(from=1, to = 0.5, length=20)
#' effort[,3] <- seq(from=1, to = 2, length=20)
#' effort[,4] <- seq(from=2, to = 1, length=20)
#' getFMortGear(params, effort=effort)
#' # Get the fishing mortality using the effort already held in a MizerSim object.
#' sim <- project(params, t_max = 20, effort = 0.5)
#' getFMortGear(sim)
#' getFMortGear(sim, time_range=c(10,20))
#' }
getFMortGear <- function(object, effort, ...){
  no_gear <- dim(object@catchability)[1]
  # If a single value, just repeat it for all gears
  if(length(effort) == 1){
    effort <- rep(effort, no_gear)
  }
  if(!is.matrix(effort)) effort <- array(effort,dim=c(1,no_gear))
  if ((length(effort) != no_gear & !is.matrix(effort)) || (!is.null(dim(effort)) && dim(effort)[2] != no_gear))
    stop("Effort must be a single value or a vector as long as the number of gears\n")
  # F = sel * q * effort
  sel_q <- sweep(object@selectivity, c(1,2), object@catchability, "*")
  # Kinda nasty! ends up with 4D array 
  fmort_gear <- aaply(effort, 1, function(x,sel_q) sweep(sel_q, c(1), x, "*"), sel_q=sel_q, .drop=FALSE)
  # fmort_gear is 4D, and first D is time with length 1
  # Drop time dimension - bit annoying because we want to keep the other dims even if they have length 1
  if(dim(fmort_gear)[1]==1){
    out <- array(fmort_gear, dim=dim(fmort_gear)[2:4])
    dimnames(out) <- dimnames(fmort_gear)[2:4]
    fmort_gear <- out
  }
  return(fmort_gear)
}


# Total fishing mortality from all gears
# species x size and maybe also by time if effort is time based

#' Get the total fishing mortality from all fishing gears by time, species and
#' size.
#' 
#' Calculates the fishing mortality from all gears by species and size at each
#' time step in the \code{effort} argument.
#' The total fishing mortality is just the sum of the fishing mortalities
#' imposed by each gear.
#' 
#' @param object A \code{MizerParams} object or a \code{MizerSim} object
#' @param effort The effort of each fishing gear. Only needed if the object
#'   argument is of class \code{MizerParams}. See notes below.
#' @param time_range Subset the returned fishing mortalities by time. The time
#'   range is either a vector of values, a vector of min and max time, or a
#'   single value. Default is the whole time range. Only used if the
#'   \code{object} argument is of type \code{MizerSim}.
#' @param drop Only used when object is of type \code{MizerSim}. Should
#'   dimensions of length 1 be dropped, e.g. if your community only has one
#'   species it might make presentation of results easier. Default is TRUE
#'
#' @return An array. If the effort argument has a time dimension, or object is
#'   of class \code{MizerSim}, the output array has three dimensions (time x
#'   species x size). If the effort argument does not have a time dimension, the
#'   output array has two dimensions (species x size).
#' @note Here: fishing mortality = catchability x selectivity x effort.
#'
#' The \code{effort} argument is only used if a \code{MizerParams} object is
#' passed in. The \code{effort} argument can be a two dimensional array (time x
#' gear), a vector of length equal to the number of gears (each gear has a
#' different effort that is constant in time), or a single numeric value (each
#' gear has the same effort that is constant in time). The order of gears in the
#' \code{effort} argument must be the same the same as in the \code{MizerParams}
#' object.
#'
#' If the object argument is of class \code{MizerSim} then the effort slot of the \code{MizerSim} object is used and the \code{effort} argument is not used.
#' @export
#' @seealso \code{getFMortGear}, \code{project}
#' @examples
#' \dontrun{
#' data(NS_species_params_gears)
#' data(inter)
#' params <- MizerParams(NS_species_params_gears, inter)
#' # Get the total fishing mortality when effort is constant for all gears and time:
#' getFMort(params, effort = 1)
#' # Get the total fishing mortality when effort is different
#' # between the four gears but constant in time:
#' getFMort(params, effort = c(0.5,1,1.5,0.75))
#' # Get the total fishing mortality when effort is different
#' # between the four gears and changes with time:
#' effort <- array(NA, dim = c(20,4))
#' effort[,1] <- seq(from=0, to = 1, length=20)
#' effort[,2] <- seq(from=1, to = 0.5, length=20)
#' effort[,3] <- seq(from=1, to = 2, length=20)
#' effort[,4] <- seq(from=2, to = 1, length=20)
#' getFMort(params, effort=effort)
#' # Get the total fishing mortality using the effort already held in a MizerSim
#' object.
#' sim <- project(params, t_max = 20, effort = 0.5)
#' getFMort(sim)
#' getFMort(sim, time_range = c(10,20))
#' }
getFMort <- function(object, effort, ...){
  fMortGear <- getFMortGear(object, effort, ...)
  if(is.matrix(effort)) marg <-c(1,3,4) else marg <- c(2,3)
  fMort <- apply(fMortGear, marg, sum)
  return(fMort)
}


# get total Z
#' getZ method for the size based model
#'
#' Calculates the total mortality \eqn{\mu_i(w)} on each species by size from
#' predation mortality (M2), background mortality (M) and fishing mortality for
#' a single time step.
#' @param object A \code{MizerParams} object.
#' @param n A matrix of species abundance (species x size).
#' @param n_pp A vector of the background abundance by size.
#' @param effort A numeric vector of the effort by gear or a single numeric
#'   effort value which is used for all gears.
#' @param m2 A two dimensional array of predation mortality (optional). Has
#'   dimensions no. sp x no. size bins in the community. If not supplied is
#'   calculated using the \code{getM2()} method.
#' @param tau The acitivity level
#'
#' @return A two dimensional array (prey species x prey size). 
#'
#' @export
#' @seealso \code{\link{getM2}}, \code{\link{getFMort}}
#' @examples
#' \dontrun{
#' data(NS_species_params_gears)
#' data(inter)
#' params <- MizerParams(NS_species_params_gears, inter)
#' # Project with constant fishing effort for all gears for 20 time steps
#' sim <- project(params, t_max = 20, effort = 0.5)
#' # Get the total mortality at a particular time step
#' getZ(params,sim@@n[21,,],sim@@n_pp[21,],effort=0.5)
#' }
getZ <- function(object, n, n_pp, effort, m2, tau){
  if (!all(dim(m2) == c(nrow(object@species_params),length(object@w)))){
    stop("m2 argument must have dimensions: no. species (",nrow(object@species_params),") x no. size bins (",length(object@w),")")
  }
  f_mort <- getFMort(object, effort = effort)
  if(missing(m2)) m2 <- getM2(object, n=n, n_pp=n_pp, tau=tau)
  zz <- m2 + f_mort
  z = sweep(zz,1,object@species_params$z0,"+")
  return(z)
}

# get tau
#' getTau method for the size based model
#'
#' Calculates the activity level tau at tau_{t+dt} based on e/mu.
#' @param object A \code{MizerParams} object.
#' @param n A matrix of species abundance (species x size).
#' @param n_pp A vector of the background abundance by size.
#' @param tau The acitivity level at t
#' @param delta The time-scale of the adaptive dynamics.
#' @param dt Time step of the model
#'
#' @return A two dimensional array (prey species x prey size). 
#'
#' @export
#' @seealso \code{\link{getEReproAndGrowth}}, \code{\link{getFeedingLevel}}
#' @examples
#' \dontrun{
#' data(NS_species_params_gears)
#' data(inter)
#' params <- MizerParams(NS_species_params_gears, inter)
#' # Project with constant fishing effort for all gears for 20 time steps
#' sim <- project(params, t_max = 20, effort = 0.5)
#' # Get tau
#' getTau(params,sim@@n[21,,],sim@@n_pp[21,])
#' }
#' 
getTau <- function(object, n, n_pp, tau, delta, dt,temp_correct=F,temp=15){
  
  #get gradient and calculate tau + dt
  F_prime <- grad(getFitness, x = tau, object=object, n=n, n_pp=n_pp, method="simple")     
  expp <- exp(dt*delta*F_prime) 
  tau_new <- tau*expp/(tau*(expp-1)+1)
  
  # manual limits
  tau_new[is.infinite(expp)] <- 0.999
  tau_new[F_prime==0 ] <- tau[F_prime==0 ]
  
  tau_new[tau_new<0.001] <- 0.001
  
  tau_new[tau_new>0.999] <- 0.999
  
  if(temp_correct == T){
    
    phi_prey <- getPhiPrey(object, n=n, n_pp=n_pp, tau=tau_new)
    O_2 <- object@O2_supply(temp,
                            O2crit = object@species_params$O2_crit,
                            P50 = object@species_params$O_2_P50,
                            Tmax = object@species_params$Tmax,
                            Topt = object@species_params$Topt,
                            om = object@species_params$om,
                            de = object@species_params$de,
                            dt = dt
    )
    tcorr <- object@get_tc(temp = temp)
    tau_max <- eval_tau_max_temp(f = O_2,
                                 tc = tcorr,
                                 omega = object@species_params$omega,
                                 gamma = object@search_vol*phi_prey,
                                 delta = object@species_params$k,
                                 h = object@species_params$h,
                                 phi = object@species_params$phi,
                                 alpha = object@species_params$alpha,
                                 k = object@species_params$ks,
                                 p = object@species_params$p,
                                 q = object@species_params$q,
                                 n = object@species_params$n,
                                 m = object@w
    )

    tau_new[tau_new>tau_max] <- tau_max[tau_new>tau_max]
    
  }

  tau_new
}

eval_tau_max_temp <- function(f=object@O2_supply(seq(1,30,l=100)),
                              tc = object@get_tc(seq(1,30,l=100)),
                              omega = 0.4,
                              gamma=50,
                              delta=2,
                              h=30,
                              phi=0.15,
                              alpha=0.25,
                              k=2,
                              p=0.8,
                              q=0.9,
                              n=0.8,
                              m=100){
  
  m <- matrix(m,dim(gamma)[1], dim(gamma)[2], byrow=T)
  tau_max <- -1/2*(delta*h*k*m^(n - p + q)*omega*tc^2 - f*gamma*m^n + (alpha*gamma*h*m^q + gamma*k*m^n)*omega*tc - sqrt(delta^2*h^2*k^2*m^(2*n - 2*p + 2*q)*omega^2*tc^4 + 2*(alpha*delta*gamma*h^2*k*m^(n - p + 2*q) - delta*gamma*h*k^2*m^(2*n - p + q))*omega^2*tc^3 + f^2*gamma^2*m^(2*n) - 2*(alpha*f*gamma^2*h*m^(n + q) + f*gamma^2*k*m^(2*n))*omega*tc + (2*delta*f*gamma*h*k*m^(2*n - p + q)*omega + (alpha^2*gamma^2*h^2*m^(2*q) + 2*alpha*gamma^2*h*k*m^(n + q) + gamma^2*k^2*m^(2*n))*omega^2)*tc^2))/(delta*gamma*k*m^n*omega*tc)  
  tau_max
}


#' getFitness method for the size based model
#'
#' Calculates the Fitness as e/mu.
#' @param object A \code{MizerParams} object.
#' @param n A matrix of species abundance (species x size).
#' @param n_pp A vector of the background abundance by size.
#' @param tau The acitivity level
#' @return A two dimensional array (prey species x prey size). 
#'
#' @export
#' @seealso \code{\link{getEReproAndGrowth}}, \code{\link{getFeedingLevel}}
#' @examples
#' \dontrun{
#' data(NS_species_params_gears)
#' data(inter)
#' params <- MizerParams(NS_species_params_gears, inter)
#' # Project with constant fishing effort for all gears for 20 time steps
#' sim <- project(params, t_max = 20, effort = 0.5)
#' # Get tau
#' getFitness(sim@@tau,params,sim@@n[21,,],sim@@n_pp[21,])
#' }
getFitness <- function(tau, object, n, n_pp){
  
  # energy in
  feeding_level <- getFeedingLevel(object, n=n, n_pp=n_pp, tau = tau)
  e <- feeding_level * object@intake_max * (1-object@species_params$alpha-object@species_params$phi)
  # Subtract basal metabolism and activity 
  e <- e - object@std_metab - tau*object@activity*object@std_metab
  e <- sweep(e,2,object@w,'/')
  # mortality
  P <- getPredRate(object, n=n, n_pp=n_pp, feeding_level=feeding_level, tau = tau)
  mu <- getM2(object, n=n, n_pp=n_pp, pred_rate = P)
  
  # Fitness
  F <- matrix(NA,dim(e)[1],dim(e)[2])
  F[!(mu==0)] <- e[!(mu==0)]/mu[!(mu==0)]
  F[mu==0] <- e[mu==0]
  F[is.nan(F)] <- 0
  F
}

# Energy after metabolism and movement
#' getEReproAndGrowth method for the size based model
#'
#' Calculates the energy available by species and size for reproduction and
#' growth after metabolism and movement have been accounted for. Used by the
#' \code{project} method for performing simulations.
#' @param object A \code{MizerParams} object.
#' @param n A matrix of species abundance (species x size).
#' @param n_pp A vector of the background abundance by size.
#' @param feeding_level The current feeding level (optional). A matrix of size
#'   no. species x no. size bins. If not supplied, is calculated internally
#'   using the \code{getFeedingLevel()} method.
#' @param tau The acitivity level
#'
#' @return A two dimensional array (species x size) 
#' @export
#' @seealso \code{\link{project}} and \code{\link{getFeedingLevel}}.
#' @examples
#' \dontrun{
#' data(NS_species_params_gears)
#' data(inter)
#' params <- MizerParams(NS_species_params_gears, inter)
#' # Project with constant fishing effort for all gears for 20 time steps
#' sim <- project(params, t_max = 20, effort = 0.5)
#' # Get the energy at a particular time step
#' getEReproAndGrowth(params,sim@@n[21,,],sim@@n_pp[21,])
#' }
getEReproAndGrowth <- function(object, n, n_pp, feeding_level, tau){
  if (!all(dim(feeding_level) == c(nrow(object@species_params),length(object@w)))){
    stop("feeding_level argument must have dimensions: no. species (",nrow(object@species_params),") x no. size bins (",length(object@w),")")
  }
  if(missing(feeding_level)) feeding_level <- getFeedingLevel(object, n=n, n_pp=n_pp, tau = tau)
  # assimilated intake
  e <- feeding_level * object@intake_max * (1-object@species_params$alpha-object@species_params$phi)
  # Subtract basal metabolism and activity 
  e <- e - (1+tau*object@activity)*object@std_metab
  e[e<0] <- 0 # Do not allow negative growth
  return(e)
}



# Energy left for reproduction
# assimilated food intake, less metabolism and activity, split between reproduction and growth

#' getESpawning method for the size based model
#'
#' Calculates the energy available by species and size for reproduction after
#' metabolism and movement have been accounted for.
#' Used by the \code{project} method for performing simulations.
#' @param object A \code{MizerParams} object.
#' @param n A matrix of species abundance (species x size).
#' @param n_pp A vector of the background abundance by size.
#' @param e The energy available for reproduction and growth (optional). A
#'   matrix of size no. species x no. size bins. If not supplied, is calculated
#'   internally using the \code{getEReproAndGrowth()} method.
#' @param tau The acitivity level
#' 
#' 
#' @return A two dimensional array (prey species x prey size) 
#' @export
#' @seealso \code{\link{project}} and \code{\link{getEReproAndGrowth}}.
#' @examples
#' \dontrun{
#' data(NS_species_params_gears)
#' data(inter)
#' params <- MizerParams(NS_species_params_gears, inter)
#' # Project with constant fishing effort for all gears for 20 time steps
#' sim <- project(params, t_max = 20, effort = 0.5)
#' # Get the energy at a particular time step
#' getESpawning(params,sim@@n[21,,],sim@@n_pp[21,])
#' }
getESpawning <- function(object, n, n_pp, e, tau=tau){
  if (!all(dim(e) == c(nrow(object@species_params),length(object@w)))){
    stop("e argument must have dimensions: no. species (",nrow(object@species_params),") x no. size bins (",length(object@w),")")
  }
  if(missing(e)) e <- getEReproAndGrowth(object,n=n,n_pp=n_pp, tau=tau)
  e_spawning <- object@psi * e 
  return(e_spawning)
}

#' getEGrowth method for the size based model
#'
#' Calculates the energy \eqn{g_i(w)} available by species and size for growth
#' after metabolism, movement and reproduction have been accounted for. Used by
#' the \code{\link{project}} method for performing simulations.
#' @param object A \linkS4class{MizerParams} object.
#' @param n A matrix of species abundance (species x size).
#' @param n_pp A vector of the background abundance by size.
#' @param e The energy available for reproduction and growth (optional, although
#'   if specified, e_spawning must also be specified). A matrix of size no.
#'   species x no. size bins. If not supplied, is calculated internally using
#'   the \code{\link{getEReproAndGrowth}} method.
#' @param e_spawning The energy available for spawning (optional, although if
#'   specified, e must also be specified). A matrix of size no. species x no.
#'   size bins. If not supplied, is calculated internally using the
#'   \code{\link{getESpawning}} method.
#' @param tau The acitivity level
#'  
#' @return A two dimensional array (prey species x prey size) 
#' @export
#' @seealso \code{\link{project}}
#' @examples
#' \dontrun{
#' data(NS_species_params_gears)
#' data(inter)
#' params <- MizerParams(NS_species_params_gears, inter)
#' # Project with constant fishing effort for all gears for 20 time steps
#' sim <- project(params, t_max = 20, effort = 0.5)
#' # Get the energy at a particular time step
#' getEGrowth(params,sim@@n[21,,],sim@@n_pp[21,])
#' }
getEGrowth <- function(object, n, n_pp, e_spawning, e, tau){
  if (!all(dim(e_spawning) == c(nrow(object@species_params),length(object@w)))){
    stop("e_spawning argument must have dimensions: no. species (",nrow(object@species_params),") x no. size bins (",length(object@w),")")
  }
  if (!all(dim(e) == c(nrow(object@species_params),length(object@w)))){
    stop("e argument must have dimensions: no. species (",nrow(object@species_params),") x no. size bins (",length(object@w),")")
  }
  # Assimilated intake less activity and metabolism
  # energy for growth is intake - energy for growth
  if(missing(e_spawning)) e_spawning <- getESpawning(object, n=n, n_pp=n_pp)
  if(missing(e)) e <- getEReproAndGrowth(object,n=n,n_pp=n_pp, tau=tau)
  e_growth <- e - e_spawning
  return(e_growth)
}

#' getRDI method for the size based model
#'
#' Calculates the density independent recruitment (total egg production)
#' \eqn{R_{p,i}} before density dependence, by species. Used by the
#' \code{project} method for performing simulations.
#' @param object A \code{MizerParams} object.
#' @param n A matrix of species abundance (species x size).
#' @param n_pp A vector of the background abundance by size.
#' @param e_spawning The energy available for spawning (optional). A matrix of
#'   size no. species x no. size bins. If not supplied, is calculated internally
#'   using the \code{\link{getESpawning}} method.
#' @param sex_ratio Proportion of the population that is female. Default value
#'   is 0.5.
#'   
#' @return A numeric vector the length of the number of species 
#' @export
#' @seealso \code{\link{project}}
#' @examples
#' \dontrun{
#' data(NS_species_params_gears)
#' data(inter)
#' params <- MizerParams(NS_species_params_gears, inter)
#' # Project with constant fishing effort for all gears for 20 time steps
#' sim <- project(params, t_max = 20, effort = 0.5)
#' # Get the recruitment at a particular time step
#' getRDI(params,sim@@n[21,,],sim@@n_pp[21,])
#' }
getRDI <- function(object, n, n_pp, e_spawning, sex_ratio = 0.5){
  if (!all(dim(e_spawning) == c(nrow(object@species_params),length(object@w)))){
    stop("e_spawning argument must have dimensions: no. species (",nrow(object@species_params),") x no. size bins (",length(object@w),")")
  }
  # Should we put this in the class as part of species_params?
  # Index of the smallest size class for each species
  #w0_idx <- as.vector(tapply(object@species_params$w_min,1:length(object@species_params$w_min),function(w_min,wx) max(which(wx<=w_min)),wx=params@w))
  if(missing(e_spawning)) e_spawning <- getESpawning(object, n=n, n_pp=n_pp)
  e_spawning_pop <- (e_spawning*n) %*% object@dw
  rdi <- sex_ratio*(e_spawning_pop * object@species_params$erepro)/object@w[object@species_params$w_min_idx] 
  return(rdi)
}


#' getRDD method for the size based model
#'
#' Calculates the density dependent recruitment (total egg production) \eqn{R_i}
#' for each species. This is the flux entering the smallest size class of each
#' species. The density dependent recruiment is the density independent
#' recruitment after it has been put through the density dependent
#' stock-recruitment relationship function. This method is used by the
#' \code{project} method for performing simulations.
#' @param object An \code{MizerParams} object
#' @param n A matrix of species abundance (species x size)
#' @param n_pp A vector of the background abundance by size
#' @param rdi A matrix of density independent recruitment (optional) with
#'   dimensions no. sp x 1. If not specified rdi is calculated internally using
#'   the \code{\link{getRDI}} method.
#' @param sex_ratio Proportion of the population that is female. Default value
#'   is 0.5
#'   
#' @return A numeric vector the length of the number of species. 
#' @export
#' @examples
#' \dontrun{
#' data(NS_species_params_gears)
#' data(inter)
#' params <- MizerParams(NS_species_params_gears, inter)
#' # Project with constant fishing effort for all gears for 20 time steps
#' sim <- project(params, t_max = 20, effort = 0.5)
#' # Get the energy at a particular time step
#' getRDD(params,sim@@n[21,,],sim@@n_pp[21,])
#' }
getRDD <- function(object, n, n_pp, rdi, sex_ratio = 0.5){
  if (!all(dim(rdi) == c(nrow(object@species_params),1))){
    stop("rdi argument must have dimensions: no. species (",nrow(object@species_params),") x 1")
  }
  if(missing(rdi)) rdi <- getRDI(object, n=n, n_pp=n_pp, sex_ratio = sex_ratio)
  rdd <- object@srr(rdi = rdi, species_params = object@species_params)
  return(rdd)
}


# get_time_elements
# internal function to get the array element references of the time dimension
# for the time based slots of a MizerSim object
# time_range can be character or numeric
# Necessary to include a slot_name argument because the effort and abundance
# slots have different time dimensions
get_time_elements <- function(sim,time_range,slot_name="n"){
  if (!(slot_name %in% c("n","effort")))
    stop("'slot_name' argument should be 'n' or 'effort'")
  if (!is(sim,"MizerSim"))
    stop("First argument to get_time_elements function must be of class MizerSim")
  time_range <- range(as.numeric(time_range))
  # Check that time range is even in object
  sim_time_range <- range(as.numeric(dimnames(slot(sim,slot_name))$time))
  if ((time_range[1] < sim_time_range[1]) | (time_range[2] > sim_time_range[2]))
    stop("Time range is outside the time range of the modell")
  time_elements <- (as.numeric(dimnames(slot(sim,slot_name))$time) >= time_range[1]) & (as.numeric(dimnames(slot(sim,slot_name))$time) <= time_range[2])
  names(time_elements) <- dimnames(slot(sim,slot_name))$time
  return(time_elements)
}

