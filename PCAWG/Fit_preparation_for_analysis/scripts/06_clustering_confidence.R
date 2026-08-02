# =============================================================================
# deltatau_min.R  --  resolution limit for tickTack timing clusters
#
# Given per-sample / per-segment info (mutations per segment M, number of
# segments S, purity, coverage, and optionally mutation rate mu), returns the
# smallest clock separation Delta-tau_min that the data could have resolved.
#
# Calibrated on simulation study EXP3 (sparse finite mixture, 2:1 segments,
# purity 0.8, coverage 60, S=20, 252 configs). See SUPP_identifiability.md.
#   P(correct K) = sigma(-1.64 + 4.12 ln dtau + 1.51 ln M)
#   single-axis invariant  z* = dtau * M^0.37
#   z*(50%) = 1.49,  z*(80%) = 2.08
#
# A single-cluster (K=1) call is read as an EQUIVALENCE statement:
#   observing one clock at this M/purity/coverage bounds any true gap below
#   Delta-tau_min --> "synchronous to within Delta-tau_min of molecular time".
# =============================================================================

# ---- calibration constants (from EXP3) --------------------------------------
.DTM <- list(
  a       = 0.367,   # empirical exponent  z* = dtau * M^a
  z50     = 1.49,    # invariant at 50% correct-K
  z80     = 2.08,    # invariant at 80% correct-K
  kappa_ref = 0.852, # SNR term at the calibration point (purity 0.8, cov 60)
  purity_ref = 0.8, coverage_ref = 60, S_ref = 20
)

# ---- 2:1 information geometry ------------------------------------------------
# clonal mutations split into multiplicity-2 (pre-gain, duplicated) and
# multiplicity-1 (post-gain). Expected mult-2 fraction:  phi(tau)=tau/(3-tau).
.phi_of_tau <- function(tau) tau / (3 - tau)
.dtau_dphi  <- function(tau) { phi <- .phi_of_tau(tau); 3 / (1 + phi)^2 }

# information-floor resolution (perfect multiplicity assignment):
#   sigma_tau = |dtau/dphi| * sqrt(phi(1-phi)/N)   ->  kappa = sigma_tau*sqrt(N)
.kappa_floor <- function(tau = 0.5) {
  phi <- .phi_of_tau(tau)
  abs(.dtau_dphi(tau)) * sqrt(phi * (1 - phi))
}

# VAF peaks for a 2:1 (total CN 3) segment at purity rho:
#   mult-1 reads at v1 = rho/(2+rho);  mult-2 at v2 = 2 rho/(2+rho)
.vpeaks <- function(rho) c(v1 = rho / (2 + rho), v2 = 2 * rho / (2 + rho))

# multiplicity mis-assignment rate at coverage DP (nearest-peak rule):
.misclass_rate <- function(rho, DP) {
  v <- .vpeaks(rho); thr <- floor(DP * (v["v1"] + v["v2"]) / 2)
  e2 <- pbinom(thr, DP, v["v2"])          # mult-2 read misassigned to mult-1
  e1 <- 1 - pbinom(thr, DP, v["v1"])      # mult-1 read misassigned to mult-2
  as.numeric(0.5 * (e1 + e2))
}

# effective SNR term: floor inflated by the binary-channel mis-assignment.
# info degrades as (1-2e)^2, so sigma (hence kappa) scales as 1/(1-2e).
kappa_eff <- function(tau = 0.5, purity = 0.8, coverage = 60) {
  e <- .misclass_rate(purity, coverage)
  .kappa_floor(tau) / max(1e-6, (1 - 2 * e))
}

# =============================================================================
# dtau_min(): the resolution limit
#
#   M         expected clonal mutations timing ONE segment (mult-1 + mult-2).
#             This is the PER-SEGMENT burden -- the axis swept in EXP3.
#   purity    tumour purity rho in (0,1)
#   coverage  mean sequencing depth on the segment
#   tau       molecular-time location of the cluster (default 0.5; the gauge
#             point where information is weakest -- a conservative choice)
#   S         number of segments supporting the cluster. EXP3 was calibrated at
#             S = 20. Two DISTINCT channels set resolution:
#               (1) per-segment tau precision, set by M  -> the fitted score
#               (2) cross-segment pooling of the shared clock,
#                   sigma_clock ~ sigma_seg / sqrt(S)     -> correction below
#             We evaluate the fitted score at per-segment M (its S=20 world),
#             then rescale by sqrt(S_ref / S) for the actual segment count. At
#             S=20 the correction is 1; fewer segments -> larger (more censored).
#   target    "80" (default, recommended for headline claims) or "50"
#
# Returns Delta-tau_min in molecular-time units [0,1].
#
# CAVEAT: the pooling correction is the theoretical 1/sqrt(S) law anchored at
# the single calibrated point S=20. For S far from 20 (e.g. 1-3 segments, as in
# many real single-cluster samples) treat the value as an approximation to be
# confirmed by a targeted simulation at the sample's own (S, M, purity).
# =============================================================================
dtau_min <- function(M, purity = 0.8, coverage = 60, tau = 0.5, S = 20,
                     target = c("80", "50")) {
  target <- match.arg(target)
  z  <- if (target == "80") .DTM$z80 else .DTM$z50
  snr_scale  <- kappa_eff(tau, purity, coverage) / .DTM$kappa_ref  # SNR vs ref
  pool_scale <- sqrt(.DTM$S_ref / max(1, S))                        # pooling vs ref
  z / M^.DTM$a * snr_scale * pool_scale
}

# invert: per-segment M needed to resolve a target gap dtau (at given S)
M_needed <- function(dtau, purity = 0.8, coverage = 60, tau = 0.5, S = 20,
                     target = c("80", "50")) {
  target <- match.arg(target)
  z <- if (target == "80") .DTM$z80 else .DTM$z50
  snr_scale  <- kappa_eff(tau, purity, coverage) / .DTM$kappa_ref
  pool_scale <- sqrt(.DTM$S_ref / max(1, S))
  ((z * snr_scale * pool_scale) / dtau)^(1 / .DTM$a)
}

# raw probability a gap dtau is resolved at burden M (2-D confidence score)
p_correct <- function(dtau, M) {
  eta <- -1.636 + 4.123 * log(dtau) + 1.513 * log(M)
  1 / (1 + exp(-eta))
}
