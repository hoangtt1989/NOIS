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

// // not currently implemented/used
// // [[Rcpp::export]]
// int qdet(const double & local_q, NumericVector nz) {
//   int ret = max(1, ceil(local_q * nz.size()));
//   return ret;
// }

// // [[Rcpp::export]]
// List NOIS_inner(NumericVector yy,
//                          double xx_inner, NumericVector xx,
//                          double first_h, NumericVector kern_nzsqrt,
//                          NumericVector kern_nzsqrtinv,
//                          int qq_inner, double tol, int maxit) {
//   int n = yy.size();
//   // declaring type
//   NumericVector yy_adj(n);
//   double cond_inner;
//   NumericVector rr(n);
//   NumericVector cond_abs(n);
//   NumericVector gamma_diff(n);
//   NumericVector thresh(n);
//   int i;
//   bool converged;
//   List ret;
//
//   // outside the loop
//   NumericMatrix gamma_curr(n, n);
//   NumericVector gamma_inner(n);
//   NumericVector gamma_next(n);
//   // std::fill(gamma_next.begin(), gamma_next.end(), 0.0);
//   // std::fill(gamma_inner.begin(), gamma_inner.end(), 0.0);
//
//
//   // inside the loop
//   for(i = 0; i < maxit; ++i) {
//     yy_adj = yy - gamma_curr(_, i);
//     rr = kern_nzsqrt*(yy - nwestimator(xx_inner, xx, yy_adj, first_h));
//     gamma_next = kern_nzsqrtinv*quantile_thresh(rr, qq_inner);
//     gamma_diff = gamma_next - gamma_curr(_ ,i);
//     cond_abs = abs(gamma_diff);
//     cond_inner = max(cond_abs);
//     gamma_curr(_, i+1) = gamma_next;
//     if (cond_inner <= tol) {
//       break;
//     }
//   }
//   // return gamma_diff;
//   ret["gamma"] = gamma_curr(_, i);
//   ret["iter"] = i;
//   ret["converged"] = converged;
//   return ret;
// }


