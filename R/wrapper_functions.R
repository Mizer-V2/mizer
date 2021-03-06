#' Sets up parameters for a community-type model
#' 
#' This functions creates a \code{\linkS4class{MizerParams}} object so that
#' community-type models can be easily set up and run. A community model has
#' several features that distinguish it from the food-web type models. Only one
#' 'species' is resolved, i.e. one 'species' is used to represent the whole
#' community. The resource spectrum only extends to the start of the community
#' spectrum. Recruitment to the smallest size in the community spectrum is
#' constant and set by the user. As recruitment is constant, the proportion of
#' energy invested in reproduction (the slot \code{psi} of the returned 
#' \code{MizerParams} object) is set to 0. Standard metabolism has been turned 
#' off (the parameter \code{ks} is set to 0). Consequently, the growth rate is 
#' now determined solely by the assimilated food (see the package vignette for 
#' more details).
#' 
#' The function has many arguments, all of which have default values. The main 
#' arguments that the users should be concerned with are \code{z0}, 
#' \code{recruitment}, \code{alpha} and \code{f0} as these determine the average
#' growth rate of the community.
#' 
#' Fishing selectivity is modelled as a knife-edge function with one parameter, 
#' \code{knife_edge_size}, which determines the size at which species are 
#' selected.
#' 
#' The resulting \code{MizerParams} object can be projected forward using 
#' \code{project()} like any other \code{MizerParams} object. When projecting 
#' the community model it may be necessary to reduce \code{dt} to 0.1 to avoid 
#' any instabilities with the solver. You can check this by plotting the biomass
#' or abundance through time after the projection.
#' 
#' @param z0 The background mortality of the community. Default value is 0.1.
#' @param alpha The assimilation efficiency of the community. Default value 0.2
#' @param f0 The average feeding level of individuals who feed mainly on the 
#'   resource. This value is used to calculate the search rate parameter 
#'   \code{gamma} (see the package vignette). Default value is 0.7.
#' @param h The maximum food intake rate. Default value is 10.
#' @param beta The preferred predator prey mass ratio. Default value is 100.
#' @param sigma The width of the prey preference. Default value is 2.0.
#' @param q The search volume exponent. Default value is 0.8.
#' @param n The scaling of the intake. Default value is 2/3.
#' @param kappa The carrying capacity of the background spectrum. Default value
#'   is 1000.
#' @param lambda The exponent of the background spectrum. Default value is 2 + q
#'   - n.
#' @param r_pp Growth rate of the primary productivity. Default value is 10.
#' @param gamma Volumetric search rate. Estimated using \code{h}, \code{f0} and 
#'   \code{kappa} if not supplied.
#' @param recruitment The constant recruitment in the smallest size class of the
#'   community spectrum. This should be set so that the community spectrum 
#'   continues the background spectrum. Default value = \code{kappa} * 
#'   \code{min_w}^-\code{lambda}.
#' @param rec_mult Additional multiplier for the constant recruitment. Default 
#'   value is 1.
#' @param knife_edge_size The size at the edge of the knife-selectivity 
#'   function. Default value is 1000.
#' @param knife_is_min Is the knife-edge selectivity function selecting above 
#'   (TRUE) or below (FALSE) the edge. Default is TRUE.
#' @param max_w The maximum size of the community. The \code{w_inf} of the 
#'   species used to represent the community is set to 0.9 * this value. The 
#'   default value is 1e6.
#' @param min_w The minimum size of the community. Default value is 1e-3.
#' @param ... Other arguments to pass to the \code{MizerParams} constructor.
#' @export
#' @return An object of type \code{\linkS4class{MizerParams}}
#' @references K. H. Andersen,J. E. Beyer and P. Lundberg, 2009, Trophic and 
#'   individual efficiencies of size-structured communities, Proceedings of the 
#'   Royal Society, 276, 109-114
#' @examples
#' \dontrun{
#' params <- set_community_model(f0=0.7, z0=0.2, recruitment=3e7)
#' sim <- project(params, effort = 0, t_max = 100, dt=0.1)
#' plotBiomass(sim)
#' plotSpectra(sim)
#' }
set_community_model <- function(max_w = 1e6,
                                min_w = 1e-3,
                                z0 = 0.1,
                                alpha = 0.2,
                                h = 10,
                                beta = 100,
                                sigma = 2.0,
                                q = 0.8,
                                n = 2/3,
                                kappa = 1000,
                                lambda = 2+q-n,
                                f0 = 0.7,
                                r_pp = 10,
                                gamma = NA,
                                knife_edge_size = 1000,
                                knife_is_min = TRUE,
                                recruitment = kappa * min_w^-lambda,
                                rec_mult = 1,
                                ...
                                ){
    w_inf <- max_w * 0.9
    w_pp_cutoff <- min_w
    ks <- 0 # Turn off standard metabolism
    p <- n # But not used as ks = 0
    # Estimate gamma if not supplied
    if (is.na(gamma)){
        gamma <- (f0 * h * beta^(2-lambda)) / ((1-f0)*sqrt(2*pi)*kappa*sigma)
    }
    # Make the species data.frame
    com_params_df <- data.frame(
        species = "Community",
        w_inf = w_inf,
        w_mat = 1e12, # Has no affect as psi set to 0 but we set it to something 
                      # to help the constructor
        h = h, # max food intake
        gamma = gamma,# vol. search rate,
        ks = ks,# standard metabolism coefficient,
        beta = beta,
        sigma = sigma,
        z0 = z0, # background mortality
        alpha = alpha,
        erepro = 1, # not used
        sel_func = "knife_edge",
        knife_edge_size = knife_edge_size,
        knife_is_min = knife_is_min,
        constant_recruitment = recruitment * rec_mult # to be used in the SRR
    )
    # Set the recruitment function for constant recruitment
    constant_recruitment <- function(rdi, species_params){
        return(species_params$constant_recruitment)
    }
    com_params <- MizerParams(com_params_df, p=p, n=n,q=q, lambda = lambda, 
                              kappa = kappa, min_w = min_w, max_w = max_w, 
                              w_pp_cutoff = w_pp_cutoff, r_pp = r_pp, ...)
    com_params@srr <- constant_recruitment
    com_params@psi[] <- 0 # Need to force to be 0. Can try setting w_mat but 
                          # due to slope still not 0
    # Set w_mat to NA for clarity - it is not actually being used
    com_params@species_params$w_mat[] <- NA
    return(com_params)
}





#' Sets up parameters for a trait-based model
#' 
#' This functions creates a \code{MizerParams} object so that trait-based-type 
#' models can be easily set up and run. The trait-based size spectrum model can
#' be derived as a simplification of the general size-based model used in
#' \code{mizer}. All the species-specific parameters are the same, except for
#' the asymptotic size, which is considered the most important trait
#' characterizing a species. Other parameters are related to the asymptotic
#' size. For example, the size at maturity is given by w_inf * eta, where eta is
#' the same for all species. For the trait-based model the number of species is
#' not important. For applications of the trait-based model see Andersen &
#' Pedersen (2010). See the \code{mizer} vignette for more details and examples
#' of the trait-based model.
#' 
#' The function has many arguments, all of which have default values. Of
#' particular interest to the user are the number of species in the model and
#' the minimum and maximum asymptotic sizes. The asymptotic sizes of the species
#' are spread evenly on a logarithmic scale within this range.
#' 
#' The stock recruitment relationship is the default Beverton-Holt style. The
#' maximum recruitment is calculated using equilibrium theory (see Andersen &
#' Pederson, 2010) and a multiplier, \code{k0}. Users should adjust \code{k0} to
#' get the spectra they want.
#' 
#' The factor for the search volume, \code{gamma}, is calculated using the
#' expected feeding level, \code{f0}.
#' 
#' Fishing selectivity is modelled as a knife-edge function with one parameter,
#' \code{knife_edge_size}, which is the size at which species are selected. Each
#' species can either be fished by the same gear (\code{knife_edge_size} has a
#' length of 1) or by a different gear (the length of \code{knife_edge_size} has
#' the same length as the number of species and the order of selectivity size is
#' that of the asymptotic size).
#' 
#' The resulting \code{MizerParams} object can be projected forward using
#' \code{project()} like any other \code{MizerParams} object. When projecting
#' the community model it may be necessary to reduce \code{dt} to 0.1 to avoid
#' any instabilities with the solver. You can check this by plotting the biomass
#' or abundance through time after the projection.
#' @param no_sp The number of species in the model. The default value is 10. The
#'   more species, the longer takes to run.
#' @param min_w_inf The asymptotic size of the smallest species in the
#'   community.
#' @param max_w_inf The asymptotic size of the largest species in the community.
#' @param no_w The number of size bins in the community spectrum.
#' @param min_w The smallest size of the community spectrum.
#' @param max_w The largest size of the community spectrum. Default value is the
#'   largest w_inf in the community x 1.1.
#' @param min_w_pp The smallest size of the background spectrum.
#' @param no_w_pp The number of the extra size bins in the background spectrum
#'   (i.e. the difference between the number of sizes bins in the community
#'   spectrum and the full spectrum).
#' @param w_pp_cutoff The cut off size of the background spectrum. Default value
#'   is 1.
#' @param k0 Multiplier for the maximum recruitment. Default value is 50.
#' @param n Scaling of the intake. Default value is 2/3.
#' @param p Scaling of the standard metabolism. Default value is 0.75.
#' @param q Exponent of the search volume. Default value is 0.9.
#' @param eta Factor to calculate \code{w_mat} from asymptotic size.
#' @param r_pp Growth rate of the primary productivity. Default value is 4.
#' @param kappa Carrying capacity of the resource spectrum. Default value is
#'   0.005.
#' @param lambda Exponent of the resource spectrum. Default value is (2+q-n).
#' @param alpha SDA
#' @param phi excretion
#' @param omega oxygen required per g of cabon metabolised
#' @param ks Standard metabolism coefficient. Default value is 4.
#' @param k Activity multiplyer. Default value is 4.
#' @param delta Time scale of adaptive foraging dynamics. Default value is 1.
#' @param z0pre The coefficient of the background mortality of the community. z0
#'   = z0pre * w_inf ^ (n-1). The default value is 0.6.
#' @param h Maximum food intake rate. Default value is 30.
#' @param beta Preferred predator prey mass ratio. Default value is 100.
#' @param sigma Width of prey size preference. Default value is 1.3.
#' @param f0 Expected average feeding level. Used to set \code{gamma}, the
#'   factor for the search volume. The default value is 0.5.
#' @param gamma Volumetric search rate. Estimated using \code{h}, \code{f0} and
#'   \code{kappa} if not supplied.
#' @param knife_edge_size The minimum size at which the gear or gears select
#'   species. Must be of length 1 or no_sp.
#' @param gear_names The names of the fishing gears. A character vector, the
#'   same length as the number of species. Default is 1 - no_sp.
#' @param ... Other arguments to pass to the \code{MizerParams} constructor.
#' @export
#' @return An object of type \code{MizerParams}
#' @seealso \linkS4class{MizerParams}
#' @references K. H. Andersen and M. Pedersen, 2010, Damped trophic cascades
#'   driven by fishing in model marine ecosystems. Proceedings of the Royal
#'   Society V, Biological Sciences, 1682, 795-802.
#' @examples
#' \dontrun{
#' trait_params <- set_trait_model(no_sp = 15)
#' init_pop <- get_initial_n(trait_params, n0_mult = 0.001)
#' sim <- project(trait_params, effort = 0, t_max = 50, dt=0.2,
#'     initial_n = init_pop, t_save = 1)
#' plot(sim)
#' ## Set up industrial fishery that only fishes on species with w_inf <= 500 g
#' ## And where the selectivity of the industrial fishery = w_inf * 0.05
#' no_sp <- 10
#' min_w_inf <- 10
#' max_w_inf <- 1e5
#' w_inf <- 10^seq(from=log10(min_w_inf), to = log10(max_w_inf), length=no_sp)
#' knife_edges <- w_inf * 0.05
#' industrial_gears <- w_inf <= 500
#' other_gears <- w_inf > 500
#' gear_names <- rep("Industrial", no_sp)
#' gear_names[other_gears] <- "Other"
#' params_gear <- set_trait_model(no_sp = no_sp, min_w_inf = min_w_inf,
#'     max_w_inf = max_w_inf, knife_edge_size = knife_edges, gear_names = gear_names)
#' ## Only turn on Industrial fishery. Set effort of the Other gear to 0
#' sim <- project(params_gear, t_max = 20, effort = c(Industrial = 1, Other = 0))
#' }
set_trait_model <- function(no_sp = 10,
                            min_w_inf = 10,
                            max_w_inf = 1e5,
                            no_w = 100,
                            min_w = 0.001,
                            max_w = max_w_inf * 1.1,
                            min_w_pp = 1e-10,
                            no_w_pp = round(no_w)*0.3,
                            w_pp_cutoff = 1,
                            k0 = 50, # recruitment adjustment parameter
                            n = 2/3,
                            p = 0.75,
                            q = 0.9, 
                            eta = 0.25,
                            r_pp = 4,
                            kappa = 0.005,
                            omega=0.4,
                            lambda = 2+q-n,
                            alpha = 0.15,
                            phi=0.25,
                            ks = 4,
                            k = 4,
                            delta = 1,
                            z0pre = 0.6,
                            h = 30,
                            beta = 100,
                            sigma = 1.3,
                            f0 = 0.5,
                            gamma = NA,
                            knife_edge_size = 1000,
                            gear_names = "knife_edge_gear",
                            O2_crit = 2,
                            O_2_P50 = 8,
                            Tmax = 30,
                            Topt = 15,
                            om = 1.87,
                            de = 1038,
                            Ea = 0.52,
                            tref = 15,
                            ...){
    # If not supplied, calculate gamma using equation 2.1 in A&P 2010
    if(is.na(gamma)){
        alpha_e <- sqrt(2*pi) * sigma * beta^(lambda-2) * exp((lambda-2)^2 * sigma^2 / 2) # see A&P 2009
        gamma <- h * f0 / (alpha_e * kappa * (1-f0)) # see A&P 2009 
    }
    w_inf <- 10^seq(from=log10(min_w_inf), to = log10(max_w_inf), length=no_sp)
    w_mat <- w_inf * eta

    # Check gears
    if (length(knife_edge_size) > no_sp){
        stop("There cannot be more gears than species in the model")
    }
    if ((length(knife_edge_size) > 1) & (length(knife_edge_size) != no_sp)){
        warning("Number of gears is less than number of species so gear information is being recycled. Is this what you want?")
    }
    if ((length(gear_names) != 1) & (length(gear_names) != no_sp)){
        stop("Length of gear_names argument must equal the number of species.")
    }

    # Make the species parameters data.frame
    trait_params_df <- data.frame(
            species = 1:no_sp,
            w_inf = w_inf,
            w_mat = w_mat,
            h = h, # max food intake
            gamma = gamma, # vol. search rate,
            ks = ks,# standard metabolism coefficient,
            k = k,# standard metabolism coefficient,
            omega= omega,
            n=n,
            p=p,
            q=q,
            delta = delta,
            beta = beta,
            sigma = sigma,
            z0 = z0pre * w_inf^(n-1), # background mortality
            alpha = alpha,
            phi = phi,
            #r_max = r_max,
            sel_func = "knife_edge",
            knife_edge_size = knife_edge_size,
            gear = gear_names,
            O2_crit = O2_crit,
            O_2_P50 = O_2_P50,
            Tmax = Tmax,
            Topt = Topt,
            om = om,
            de = de,
            erepro = 1 # not used but included out of necessity
    )
    # Make the MizerParams
    trait_params <- MizerParams(trait_params_df, 
                                min_w = min_w, 
                                max_w=max_w, 
                                no_w = no_w, 
                                min_w_pp = min_w_pp, 
                                w_pp_cutoff = w_pp_cutoff, 
                                n = n, 
                                p=p, 
                                q=q, 
                                r_pp=r_pp,
                                Ea=Ea,
                                tref=tref, 
                                kappa=kappa, 
                                lambda = lambda) 
    
    # Sort out maximum recruitment - see A&P 2009
    # Get max flux at recruitment boundary, R_max
    # R -> | -> g0 N0
    # R is egg flux, in numbers per time
    # Actual flux at recruitment boundary = RDD = NDD * g0 (where g0 is growth rate)
    # So in our BH SRR we need R_max comparable to RDI (to get RDD)
    # R_max = N0_max * g0 (g0 is the average growth rate of smallest size, i.e. at f0 = 0.5)
    # N0 given by Appendix A of A&P 2010 - see Ken's email 12/08/13
    # Taken from Ken's code 12/08/13 - equation in paper is wrong!
    alpha_p <- f0 * h * beta^(2 * n - q - 1) * exp((2 * n * (q - 1) - q^2 + 1) * sigma^2 / 2)
    alpha_rec <- alpha_p / (alpha * h * f0 - ks)
    # Calculating dw using Ken's code - see Ken's email 12/08/13
    tmpA <- w_inf[1]
    tmpB <- (log10(w_inf[length(w_inf)]) - log10(w_inf[1])) / (no_sp - 1) # Difference between logged w_infs, fine
    dw_winf <- tmpB * tmpA *10^(tmpB*((1:no_sp)-1)) # ?
    N0_max <- k0 * w_inf^(n*2-q-3+alpha_rec) * dw_winf  # Why * dw_winf, not / ? Ken confirms * in email
    # No need to include (1 - psi) in growth equation because allocation to reproduction at this size = 0, so 1 - psi = 1
    g0 <- (alpha * f0 * h * trait_params@w[1]^n - ks * trait_params@w[1]^p)
    r_max <- N0_max * g0

    trait_params@species_params$r_max <- r_max

    return(trait_params)
}

