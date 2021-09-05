#include <iostream>
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

// [[Rcpp::export]]
arma::umat OPS(int n, arma::vec prop) {

  int n_splits = prop.n_rows;
  int n_split0 = arma::sum(prop);
  int n_samples_split0 = ceil(n / double (n_split0));
  int n_samples_max = n_samples_split0 * arma::max(prop);
  arma::umat IDX(n_splits, n_samples_max);
  IDX.fill(n);
  arma::uvec count = arma::zeros<arma::uvec>(n_splits);
  for (int k = 0; k < n_samples_split0; k++) {
    arma::uvec temp = n_split0 * k + arma::randperm(n_split0);
    int count_prop = 0;
    for (int l = 0; l < n_splits; l++) {
      arma::uvec temp_l = temp(arma::span(count_prop, count_prop + prop(l) - 1));
      IDX(l, arma::span(count(l), count(l) + prop(l) - 1)) = temp_l.t();
      count_prop += prop(l);
      count(l) += prop(l);
    }
  }

  IDX = IDX + 1;
  IDX.elem(arma::find(IDX > n)).zeros();

  return(IDX);
}

// [[Rcpp::export]]
double ker_r(double u) {
  double out = 0;
  if (u >= 0 & u <= 1) {
    out = 1.5 * (1 - std::pow(u, 2));
  }
  return(out);
}

// [[Rcpp::export]]
double ker_l(double u) {
  double out = 0;
  if (u >= -1 & u < 0) {  // avoid ties
    out = 1.5 * (1 - std::pow(u, 2));
  }
  return(out);
}

// [[Rcpp::export]]
double mfit_r(double x, arma::vec Y, arma::vec X, double h) {

  int n = X.n_rows;
  double eps = arma::datum::eps;

  double s0 = 0;
  double s1 = 0;
  double s2 = 0;
  double t0 = 0;
  double t1 = 0;
  for (int i = 0; i < n; i++) {
    double Z = X(i) - x;
    double u = Z / h;
    double K = ker_r(u);
    s0 += K;
    double KZ = K * Z;
    s1 += KZ;
    s2 += KZ * Z;
    t0 += K * Y(i);
    t1 += KZ * Y(i);
  }

  double mfit = (t0 * s2 - t1 * s1) / (s0 * s2 - std::pow(s1, 2) + eps);

  return(mfit);
}

// [[Rcpp::export]]
double mfit_l(double x, arma::vec Y, arma::vec X, double h) {

  int n = X.n_rows;
  double eps = arma::datum::eps;

  double s0 = 0;
  double s1 = 0;
  double s2 = 0;
  double t0 = 0;
  double t1 = 0;
  for (int i = 0; i < n; i++) {
    double Z = X(i) - x;
    double u = Z / h;
    double K = ker_l(u);
    s0 += K;
    double KZ = K * Z;
    s1 += KZ;
    s2 += KZ * Z;
    t0 += K * Y(i);
    t1 += KZ * Y(i);
  }

  double mfit = (t0 * s2 - t1 * s1) / (s0 * s2 - std::pow(s1, 2) + eps);

  return(mfit);
}

// [[Rcpp::export]]
arma::vec delta(arma::vec x_seq, arma::vec Y, arma::vec X, double h) {
  int nx = x_seq.n_rows;
  arma::vec ts(nx);
  ts.fill(NA_REAL);
  for (int ix = 0; ix < nx; ix++) {
    double x = x_seq(ix);
    if (x > h & x < (1 - h)) {
      double mfit_r_x = mfit_r(x, Y, X, h);
      double mfit_l_x = mfit_l(x, Y, X, h);
      ts(ix) = mfit_r_x - mfit_l_x;
    }
  }
  return(ts);
}

// [[Rcpp::export]]
List detect(arma::vec x_seq, arma::vec Y, arma::vec X, double h, double c_h, int L) {

  List out;

  arma::vec jump_loc = arma::zeros(L);
  arma::vec jump_size = arma::zeros(L);

  arma::vec x_sur = x_seq;
  arma::vec delta_sur = delta(x_seq, Y, X, h);
  arma::vec delta_abs_sur = arma::abs(delta_sur);
  for (int l = 0; l < L; l++) {
    if (x_sur.n_rows > 0) {
      double delta_abs_sur_max = delta_abs_sur.max();
      for (int i = 0; i < x_sur.n_rows; i++) {
        if (delta_abs_sur(i) == delta_abs_sur_max) {
          jump_loc(l) = x_sur(i);
          jump_size(l) = delta_sur(i);
        }
      }
      arma::uvec idx_sur = arma::find(arma::abs(x_sur - jump_loc(l)) >= c_h * h);
      x_sur = x_sur(idx_sur);
      delta_sur = delta_sur(idx_sur);
      delta_abs_sur = delta_abs_sur(idx_sur);
    }
  }

  arma::uvec idx = arma::find(jump_loc != 0);
  jump_loc = jump_loc(idx);
  jump_size = jump_size(idx);
  out["jump_loc"] = jump_loc;
  out["jump_size"] = jump_size;

  return(out);
}

// [[Rcpp::export]]
double ker_c(double u) {
  double out = 0;
  if (std::abs(u) <= 1) {
    out = 0.75 * (1 - std::pow(u, 2));
  }
  return(out);
}

// [[Rcpp::export]]
double mfit_c(double x, arma::vec Y, arma::vec X, double h) {

  int n = X.n_rows;
  double eps = arma::datum::eps;

  double s0 = 0;
  double s1 = 0;
  double s2 = 0;
  double t0 = 0;
  double t1 = 0;
  for (int i = 0; i < n; i++) {
    double Z = X(i) - x;
    double u = Z / h;
    double K = ker_c(u);
    s0 += K;
    double KZ = K * Z;
    s1 += KZ;
    s2 += KZ * Z;
    t0 += K * Y(i);
    t1 += KZ * Y(i);
  }

  double mfit = (t0 * s2 - t1 * s1) / (s0 * s2 - std::pow(s1, 2) + eps);

  return(mfit);
}

// [[Rcpp::export]]
List JIC(arma::vec x_seq, arma::vec Y, arma::vec X, double h, double c_h, int L, double gamma, arma::vec adj_fac_JIC, arma::vec adj_fac_BIC) {

  List out;

  List out_detect = detect(x_seq, Y, X, h, c_h, L);
  arma::vec jump_loc = out_detect["jump_loc"];
  arma::vec jump_size = out_detect["jump_size"];
  L = jump_loc.n_rows;

  int n = X.n_rows;
  arma::mat mfit_J = arma::zeros(n, L);
  arma::mat mfit_C = arma::zeros(n, L);

  for (int l = 0; l < L; l++) {  // l = 1, ..., L, USE l + 1 if l is needed

    arma::vec jump_loc_l = jump_loc(arma::span(0, l));
    arma::vec jump_size_l = jump_size(arma::span(0, l));

    for (int i = 0; i < n; i++) {
      double x = X(i);
      for (int j = 0; j <= l; j++) {  // j = 1, ..., l
        if (x > jump_loc_l(j)) {
          mfit_J(i, l) += jump_size_l(j);
        }
      }
    }

    arma::vec Y_adj = Y - mfit_J.col(l);
    for (int i = 0; i < n; i++) {
      double x = X(i);
      mfit_C(i, l) = mfit_c(x, Y_adj, X, h);
    }
  }

  arma::mat mfit = mfit_C + mfit_J;

  arma::mat JIC = arma::zeros(L, adj_fac_JIC.n_rows);
  arma::mat BIC = arma::zeros(L, adj_fac_BIC.n_rows);
  for (int l = 0; l < L; l++) {
    double rss_l = arma::sum(arma::pow(Y - mfit.col(l), 2));
    arma::vec jump_size_l = jump_size(arma::span(0, l));
    double pen_l = arma::sum(1 / arma::pow(arma::abs(jump_size_l), gamma));
    double obj_l = n * std::log(rss_l / n);
    for (int i = 0; i < adj_fac_JIC.n_rows; i++) {
      JIC(l, i) = obj_l + adj_fac_JIC(i) * pen_l;
    }
    for (int i = 0; i < adj_fac_BIC.n_rows; i++) {
      BIC(l, i) = obj_l + adj_fac_BIC(i) * (l + 1);
    }
  }

  arma::urowvec Kn_JIC = arma::index_min(JIC, 0) + 1;
  arma::urowvec Kn_BIC = arma::index_min(BIC, 0) + 1;

  out["Kn_JIC"] = Kn_JIC;
  out["Kn_BIC"] = Kn_BIC;
  out["JIC"] = JIC;
  out["BIC"] = BIC;
  out["jump_loc"] = jump_loc;
  out["mfit"] = mfit;
  out["mfit_J"] = mfit_J;

  return(out);
}

// [[Rcpp::export]]
List val(arma::vec x_seq, arma::vec Y_O, arma::vec X_O, arma::vec Y_E, arma::vec X_E, double h, double c_h, int L) {

  List out;

  List out_detect = detect(x_seq, Y_O, X_O, h, c_h, L);
  arma::vec th_O = out_detect["jump_loc"];
  arma::vec dh_O = out_detect["jump_size"];
  arma::vec dh_E = delta(th_O, Y_E, X_E, h);
  L = th_O.n_rows;

  arma::vec val_E = arma::zeros(L);
  for (int l = 0; l < L; l++) {
    arma::vec dh_O_zero = dh_O;
    if (l < L - 1) {
      dh_O_zero(arma::span(l + 1, L - 1)).fill(0);
    }
    val_E(l) = arma::sum(arma::pow(dh_O_zero - dh_E, 2));
  }

  out["val_E"] = val_E;
  out["th_O"] = th_O;
  out["dh_O"] = dh_O;
  out["dh_E"] = dh_E;

  return(out);
}

// [[Rcpp::export]]
List RSS_h(
    arma::vec jump_loc, arma::vec jump_size,
    arma::vec Y_O, arma::vec X_O,
    arma::vec Y_H, arma::vec X_H,
    double h
) {

  List out;

  int L = jump_loc.n_rows;

  int n_O = X_O.n_rows;
  arma::vec mfit_J_O = arma::zeros(n_O);
  for (int i = 0; i < n_O; i++) {
    double x = X_O(i);
    for (int l = 0; l < L; l++) {
      if (x > jump_loc(l)) {
        mfit_J_O(i) += jump_size(l);
      }
    }
  }

  int n_H = X_H.n_rows;
  arma::vec mfit_J_H = arma::zeros(n_H);
  for (int i = 0; i < n_H; i++) {
    double x = X_H(i);
    for (int l = 0; l < L; l++) {
      if (x > jump_loc(l)) {
        mfit_J_H(i) += jump_size(l);
      }
    }
  }

  arma::vec Y_O_adj = Y_O - mfit_J_O;
  arma::vec Y_H_adj = Y_H - mfit_J_H;
  arma::vec mfit_C_H = arma::zeros(n_H);
  for (int i = 0; i < n_H; i++) {
    double x = X_H(i);
    mfit_C_H(i) = mfit_c(x, Y_O_adj, X_O, h);
  }
  double RSS = arma::sum(arma::pow(Y_H_adj -  mfit_C_H, 2));

  out["RSS"] = RSS;
  out["Y_H_adj"] = Y_H_adj;
  out["mfit_C_H"] = mfit_C_H;
  out["Y_O_adj"] = Y_O_adj;

  return(out);
}

// [[Rcpp::export]]
List sFDP(arma::vec W, arma::vec q_all) {

  List out;

  arma::vec W_ext = arma::sort(arma::abs(W));
  arma::mat FDP = arma::zeros(W_ext.n_rows, 2);
  for (int i = 0; i < W_ext.n_rows; i++) {
    double s = W_ext(i);
    arma::uvec temp_D = arma::find(W >= s);
    int D0 = temp_D.n_rows;
    int D = std::max(D0, 1);
    arma::uvec temp_FD_app = arma::find(W <= -s);
    int FD_app = temp_FD_app.n_rows;
    int FD_app_mod = 1 + FD_app;
    FDP(i, 0) = 1.0 * FD_app / D;
    FDP(i, 1) = 1.0 * FD_app_mod / D;
  }

  arma::mat thld = arma::zeros(q_all.n_rows, FDP.n_cols);
  for (int iq = 0; iq < q_all.n_rows; iq++) {
    double q = q_all(iq);
    for (int iW = 0; iW < FDP.n_cols; iW++) {
      arma::uvec ok = arma::find(FDP.col(iW) <= q);
      if (ok.is_empty()) {
        thld(iq, iW) = arma::datum::inf;
      } else {
        thld(iq, iW) = W_ext(ok(0));
      }
    }
  }

  out["thld"] = thld;
  out["FDP"] = FDP;
  out["s"] = W_ext;

  return(out);
}
