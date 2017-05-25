#include <Rcpp.h>
#include <queue>
using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//
//

// [[Rcpp::export]]
double nwestimator(const double & inputval, const NumericVector & xvals, const NumericVector & yvals, const double & h) {
  NumericVector kernel_val = dnorm(inputval - xvals, 0, h);
  double est_value = sum(kernel_val*yvals)/sum(kernel_val);
  return est_value;
}

// [[Rcpp::export]]
NumericVector nwvector(const NumericVector & x, const NumericVector & y, const double & h) {
  int n = x.size();
  NumericVector result(n);
  for(int i = 0; i < n; ++i) {
    result[i] = nwestimator(x[i], x, y, h);
  }
  return result;
}

// LOOCV
double LOOCV_loop(const int & input, const IntegerVector & inds, const NumericVector & x, const NumericVector & y, const double & h) {
  double est_val = nwestimator(x[input], x[inds != input], y[inds != input], h);
  double sq_error = std::pow(y[input] - est_val, 2.0);
  return sq_error;
}

// [[Rcpp::export]]
double LOOCV(const NumericVector & x, const NumericVector & y, const double & h, const IntegerVector & ind_keep) {
  int lenx = x.size();
  IntegerVector inds = seq(0, lenx - 1);
  NumericVector cv(lenx);
  for(int i = 0; i < lenx; ++i) {
    cv[i] = LOOCV_loop(i, inds, x, y, h);
  }
  cv = cv[ind_keep];
  double ret_val = sum(cv)/lenx;
  return ret_val;
}


// [[Rcpp::export]]
double biasnwestimator(const double & inputval, const NumericVector & xvals, const NumericVector & yvals,
                       const double & h, const double & inputnw, const NumericVector & nwvals, bool shift_sq = false) {
  NumericVector kernel_val = dnorm(inputval - xvals, 0, h);
  NumericVector shift = yvals - nwvals + inputnw;
  if (shift_sq == TRUE) {
    shift = pow(shift, 2.0);
  }
  double est_value = sum(kernel_val * shift)/sum(kernel_val);
  return est_value;
}

// [[Rcpp::export]]
NumericVector biasnwvector(const NumericVector & x, const NumericVector & y, const NumericVector & nwvals, const double & h) {
  int n = x.size();
  NumericVector result(n);
  for(int i = 0; i < n; ++i) {
    result[i] = biasnwestimator(x[i], x, y, h, nwvals[i], nwvals);
  }
  return result;
}


// functions for index sorting
template <int RTYPE>
class IndexComparator {
public:
  typedef typename Rcpp::traits::storage_type<RTYPE>::type STORAGE;

  IndexComparator(const Vector<RTYPE>& data_) : data(data_.begin()) {}

  inline bool operator()(int i, int j) const {
    return data[i] > data[j] || (data[i] == data[j] && j > i);
  }

private:
  const STORAGE* data;
};

template <>
class IndexComparator<STRSXP> {
public:
  IndexComparator( const CharacterVector& data_ ) : data(data_.begin()) {}

  inline bool operator()(int i, int j) const {
    return (String)data[i] > (String)data[j] || (data[i] == data[j] && j > i );
  }

private:
  const Vector<STRSXP>::const_iterator data;
};

template <int RTYPE>
class IndexQueue {
public:
  typedef std::priority_queue<int, std::vector<int>, IndexComparator<RTYPE> > Queue;

  IndexQueue(const Vector<RTYPE>& data_) : comparator(data_), q(comparator), data(data_) {}

  inline operator IntegerVector() {
    int n = q.size();
    IntegerVector res(n);
    for( int i=0; i<n; i++) {
      // +1 for 1-based R indexing
      res[i] = q.top() + 1;
      q.pop();
    }
    return res;
  }
  inline void input( int i) {
    // if( data[ q.top() ] < data[i] ) {
    if (comparator(i, q.top())) {
      q.pop();
      q.push(i);
    }
  }
  inline void pop() { q.pop(); }
  inline void push(int i) { q.push(i); }

private:
  IndexComparator<RTYPE> comparator;
  Queue q;
  const Vector<RTYPE>& data;
};


template <int RTYPE>
IntegerVector top_index(Vector<RTYPE> v, int n) {
  int size = v.size();

  // not interesting case. Less data than n
  if( size < n){
    return seq( 0, n-1 );
  }

  IndexQueue<RTYPE> q( v );
  for( int i=0; i<n; i++) q.push(i);
  for( int i=n; i<size; i++) q.input(i);
  return q;
}

// [[Rcpp::export]]
IntegerVector top_index(SEXP x, int n) {
  switch (TYPEOF(x)) {
  case INTSXP: return top_index<INTSXP>(x, n);
  case REALSXP: return top_index<REALSXP>(x, n);
  case STRSXP: return top_index<STRSXP>(x, n);
  default: stop("type not handled");
  }
  return IntegerVector() ; // not used
}


// [[Rcpp::export]]
NumericVector quantile_thresh(NumericVector x, const int & thresh_val) {
  NumericVector zs(x.size());
  std::fill(zs.begin(), zs.end(), 0.0);
  NumericVector x_abs = abs(x);
  IntegerVector keep_index = top_index(x_abs, thresh_val) - 1;
  zs[keep_index] = x[keep_index];
  return zs;
}

// [[Rcpp::export]]
List NOIS_loop(const NumericVector & xx, const NumericVector & yy, const double & first_h,
               const double & local_q, const double & tol, const int & maxit) {

  int nn = xx.size();
  //declaring type
  double xx_inner;
  NumericVector yy_adj(nn);
  NumericVector kernlist(nn);
  NumericVector kernlist_nzind;
  NumericVector kern_nzsqrt(nn);
  NumericVector kern_nzsqrt_sub;
  NumericVector kern_nzsqrtinv(nn);
  NumericVector kern_nzsqrtinv_sub;
  IntegerVector nz_ind;
  IntegerVector nz_tmp = seq_len(nn) - 1;
  IntegerVector qq(nn);
  int qq_inner;
  NumericVector rr(nn);
  NumericVector cond_abs(nn);
  NumericVector gamma_diff(nn);
  NumericVector thresh(nn);
  double cond_inner = 0.0;
  NumericVector cond_check(nn);
  IntegerVector iter(nn);
  int i;
  IntegerVector converged(nn);

  List ret;

  double local_inner = 0.0;
  NumericVector local_fit(nn);
  NumericMatrix gamma_curr(nn, nn);
  NumericVector gamma_inner(nn);
  NumericVector gamma_next(nn);

  for(int jj = 0; jj < nn; ++jj) {

    xx_inner = xx[jj];
    kernlist = dnorm(xx_inner - xx, 0.0, first_h);
    nz_ind = nz_tmp;
    nz_ind = nz_ind[(kernlist != 0.0) & (kernlist >= 1e-20)];
    kernlist_nzind = kernlist[nz_ind];
    std::fill(kern_nzsqrt.begin(), kern_nzsqrt.end(), 0.0);
    kern_nzsqrt_sub = sqrt(kernlist_nzind);
    kern_nzsqrt[nz_ind] = kern_nzsqrt_sub;
    std::fill(kern_nzsqrtinv.begin(), kern_nzsqrtinv.end(), 0.0);
    kern_nzsqrtinv_sub = 1.0/kern_nzsqrt_sub;
    kern_nzsqrtinv[nz_ind] = kern_nzsqrtinv_sub;
    qq[jj] = ceil(local_q * nz_ind.size());
    qq_inner = qq[jj];

    for(i = 0; i < maxit; ++i) {
      yy_adj = yy - gamma_curr(_, jj);
      local_inner = nwestimator(xx_inner, xx, yy_adj, first_h);
      rr = kern_nzsqrt * (yy - local_inner);
      gamma_next = kern_nzsqrtinv * quantile_thresh(rr, qq_inner);
      gamma_diff = gamma_next - gamma_curr(_, jj);
      cond_abs = abs(gamma_diff);
      cond_inner = max(cond_abs);
      gamma_curr(_, jj) = gamma_next;
      if (cond_inner <= tol) {
        converged[jj] = 1;
        break;
      }
    }
    local_fit[jj] = local_inner;
    cond_check[jj] = cond_inner;
    iter[jj] = i;
  }

  ret["local_fit"] = local_fit;
  ret["gamma_curr"] = gamma_curr;
  ret["iter"] = iter + 1;
  ret["cond_check"] = cond_check;
  ret["converged"] = converged;
  ret["qq_inner"] = qq;
  return ret;

}

