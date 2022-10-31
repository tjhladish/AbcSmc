#ifndef PLS_H
#define PLS_H

#include "AbcUtil.h"

//using namespace std;
using namespace Eigen;
using namespace ABC;

using std::cout;
using std::cerr;
using std::complex;
using std::endl;

typedef enum { KERNEL_TYPE1, KERNEL_TYPE2 } METHOD;
typedef enum { PRESS, RMSEP } VALIDATION_OUTPUT;
typedef enum { LOO, NEW_DATA } VALIDATION_METHOD;

/*
 *   Variable definitions from source paper:
 *     X     : predictor variables matrix (N × K)
 *     Y     : response variables matrix (N × M)
 *     B_PLS : PLS regression coefficients matrix (K × M)
 *     W     : PLS weights matrix for X (K × A)
 *     P     : PLS loadings matrix for X (K × A)
 *     Q     : PLS loadings matrix for Y (M × A)
 *     R     : PLS weights matrix to compute scores T directly from original X (K × A)
 *     T     : PLS scores matrix of X (N × A)
 *     w_a   : a column vector of W
 *     p_a   : a column vector of P
 *     q_a   : a column vector of Q
 *     r_a   : a column vector of R
 *     t_a   : a column vector of T
 *     K     : number of X-variables
 *     M     : number of Y-variables
 *     N     : number of objects
 *     A     : number of components in PLS model
 *     a     : integer counter for latent variable dimension
 */

class PLS_Model {
  public:
    Mat2Dc P, W, R, Q, T;
    int A;
    METHOD method;
    void initialize(int num_predictors, int num_responses, int num_components) {
        A = num_components;
        P.setZero(num_predictors, num_components);
        W.setZero(num_predictors, num_components);
        R.setZero(num_predictors, num_components);
        Q.setZero(num_responses, num_components);
        // T will be initialized if needed
        return;
    }

    //"Modified kernel algorithms 1 and 2"
    //from Dayal and MacGregor (1997) "Improved PLS Algorithms" J. of Chemometrics. 11,73-85.
    void plsr(const Mat2D X, const Mat2D Y, METHOD algorithm) {
        method = algorithm;
        int M = Y.cols(); // Number of response variables == columns in Y

        if (algorithm == KERNEL_TYPE1) T.setZero(X.rows(), A);

        Mat2D XY = X.transpose() * Y;
        Mat2D XX;
        if (algorithm == KERNEL_TYPE2) XX = X.transpose() * X;

        for (int i=0; i<A; i++) {
            Colc w, p, q, r, t;
            complex<float_type> tt;
            if (M==1) {
                w = XY.cast<complex<float_type> >();
            } else {
                EigenSolver<Mat2D> es( (XY.transpose() * XY) );
                q = dominant_eigenvector(es);
                w = (XY*q);
            }

            w /= sqrt((w.transpose()*w)(0,0)); // use normalize function from eigen?
            r=w;
            for (int j=0; j<=i-1; j++) {
                r -= (P.col(j).transpose()*w)(0,0)*R.col(j);
            }
            if (algorithm == KERNEL_TYPE1) {
                t = X*r;
                tt = (t.transpose()*t)(0,0);
                p.noalias() = (X.transpose()*t);
            } else if (algorithm == KERNEL_TYPE2) {
                tt = (r.transpose()*XX*r)(0,0);
                p.noalias() = (r.transpose()*XX).transpose();
            }
            p /= tt;
            q.noalias() = (r.transpose()*XY).transpose(); q /= tt;
            XY -= ((p*q.transpose())*tt).real(); // is casting this to 'real' always safe?
            W.col(i)=w;
            P.col(i)=p;
            Q.col(i)=q;
            R.col(i)=r;
            if (algorithm == KERNEL_TYPE1) T.col(i) = t;
        }
        return;
    }


    // latent X values, i.e. the orthogonal metrics you wish you could measure
    const Mat2Dc scores(Mat2D X_new) { return scores(X_new, A); }
    const Mat2Dc scores(Mat2D X_new, int comp) {
        assert (A >= comp);
        return X_new * R.leftCols(comp);
    }

    // compute the regression coefficients (aka 'beta')
    const Mat2Dc coefficients() { return coefficients(A); }
    const Mat2Dc coefficients(int comp) {
        assert (A >= comp);
        return R.leftCols(comp)*Q.leftCols(comp).transpose();
    }

    // predicted Y values, given X values and pls model
    const Mat2D fitted_values(const Mat2D& X) { return fitted_values(X, A); }
    const Mat2D fitted_values(const Mat2D& X, int comp) {
        assert (A >= comp);
        return X*coefficients(comp).real();
    }

    // unexplained portion of Y values
    const Mat2D residuals(const Mat2D& X, const Mat2D& Y) { return residuals(X, Y, A); }
    const Mat2D residuals(const Mat2D& X, const Mat2D& Y, int comp) {
        assert (A >= comp);
        return Y - fitted_values(X, comp);
    }

    // Sum of squared errors
    Row SSE(const Mat2D& X, const Mat2D& Y) { return this->SSE(X, Y, A); }
    Row SSE(const Mat2D& X, const Mat2D& Y, int comp) {
        return residuals(X, Y, comp).colwise().squaredNorm();
    }

    // Total sum of squares
    Row SST(const Mat2D& Y) {
        Row sst(Y.cols());
        for (int c = 0; c < Y.cols(); c++) {
            sst(c) = (Y.col(c).array() - (Y.col(c).sum()/Y.rows())).square().sum();
        }
        return sst;
    }

    // fraction of explainable variance
    Row explained_variance(const Mat2D& X, const Mat2D& Y) { return explained_variance(X, Y, A); }
    Row explained_variance(const Mat2D& X, const Mat2D& Y, int comp) {
        assert (A >= comp);
    //    cerr << "ev: " << this->SSE(X, Y, comp).cwiseQuotient( SST(Y) ) << endl;
        return (1.0 - this->SSE(X, Y, comp).cwiseQuotient( SST(Y) ).array()).matrix();
    }

    // leave-one-out validation of model (i.e., are we overfitting?)
    Mat2D loo_validation(const Mat2D& X, const Mat2D& Y, VALIDATION_OUTPUT out_type) {
        Mat2D Xv = X.bottomRows(X.rows()-1);
        Mat2D Yv = Y.bottomRows(Y.rows()-1);

        Mat2D SSEv = Mat2D::Zero(Y.cols(), this->A);

        PLS_Model plsm_v;
        plsm_v.initialize(Xv.cols(), Yv.cols(), this->A);
        for (int i = 0; i < X.rows()-1; i++) {
            // run pls for the data, less one observation
            plsm_v.plsr(Xv, Yv, this->method);
            for (int j = 1; j <= this->A; j++) {
                // now see how well the data predict the missing observation
                Row res = plsm_v.residuals(X.row(i), Y.row(i), j).row(0);
                // tally the squared errors
                SSEv.col(j-1) += res.cwiseProduct(res).transpose();
            }
            Xv.row(i) = X.row(i);
            Yv.row(i) = Y.row(i);
        }
        if ( out_type == PRESS ) {
            return SSEv;
        } else if ( out_type == RMSEP ) {
            SSEv /= X.rows();
            return SSEv.cwiseSqrt();
        }
    }


    std::vector<Mat2D> _loo_cv_error_matrix(const Mat2D& X, const Mat2D& Y) {
        Mat2D Xv = X.bottomRows(X.rows()-1);
        Mat2D Yv = Y.bottomRows(Y.rows()-1);

        // vector of error matrices(rows=Y.rows(), cols=Y.cols())
        // col = component #, row = obs #, tier = Y category
        std::vector<Mat2D> Ev(Y.cols(), Mat2D::Zero(X.rows(), this->A));

        PLS_Model plsm_v;
        plsm_v.initialize(Xv.cols(), Yv.cols(), this->A);
        for (int i = 0; i < X.rows(); i++) {
            plsm_v.plsr(Xv, Yv, this->method);
            //cerr << "Using matrix:\n" << Xv << "\n";
            //cerr << "Excluding row " << i << ":\n" << X.row(i) << "\n\n";
            for (int j = 1; j <= this->A; j++) {
                Row res = plsm_v.residuals(X.row(i), Y.row(i), j).row(0);
                for (int k = 0; k < res.size(); k++) Ev[k](i,j-1) = res(k);
            }
            if (i < Xv.rows()) {
                Xv.row(i) = X.row(i);
                Yv.row(i) = Y.row(i);
            }
        }
        return Ev;
    }


    std::vector<Mat2D> _new_data_cv_error_matrix(const Mat2D& X_new, const Mat2D& Y_new) {
        // vector of error matrices(rows=Y.rows(), cols=Y.cols())
        // col = component #, row = obs #, tier = Y category
        std::vector<Mat2D> Ev(Y_new.cols(), Mat2D::Zero(X_new.rows(), this->A));

        for (int j = 1; j <= this->A; j++) { // j is component #
            Mat2D res = residuals(X_new, Y_new, j);
            for (int k = 0; k < res.cols(); k++) { // k is Y category
                Ev[k].col(j-1) = res.col(k);
            }
        }
        return Ev;
    }

    // if val_method is LOO, X and Y should be original data
    // if val_method is NEW_DATA, X and Y should be observations not included in the original model
    Rowi optimal_num_components(const Mat2D& X, const Mat2D& Y, VALIDATION_METHOD val_method) {
        // col = component #, row = obs #, tier = Y category

        std::vector<Mat2D> errors;
        if (val_method == LOO) {
            errors = _loo_cv_error_matrix(X, Y);
        } else {
            errors = _new_data_cv_error_matrix(X, Y);
        }

        Mat2D press = Mat2D::Zero(Y.cols(), A);
        Rowi min_press_idx = Rowi::Zero(Y.cols());
        Row  min_press_val(Y.cols());
        Rowi best_comp(Y.cols());

        // Determine PRESS values
        for (unsigned int i=0; i<errors.size(); i++) {    // for each Y category
            for (int j=0; j<errors[i].rows(); j++) {      // for each observation
                for (int k=0; k<errors[i].cols(); k++) {  // for each component
                    press(i,k) += pow(errors[i](j,k), 2);
                }
            }
        }

        min_press_val = press.col(0);
        // Find the component number that minimizes PRESS for each Y category
        for (int i=0; i<press.rows(); i++) {              // for each Y category
            for (int j=0; j<press.cols(); j++) {          // for each component
                if (press(i,j) < min_press_val(i)) {
                    min_press_val(i) = press(i,j);
                    min_press_idx(i) = j;
                }
            }
        }

        best_comp = min_press_idx.array() + 1; // +1 to convert from index to component number
        // Find the min number of components that is not significantly
        // different from the min PRESS at alpha = 0.1 for each Y category
        const float ALPHA = 0.1;
        for (int i=0; i<press.rows(); i++) {              // for each Y category
            for (int j=0; j<min_press_idx(i); j++) {      // for each smaller number of components
                Col err1 = errors[i].col(min_press_idx(i));
                Col err2 = errors[i].col(j);
                float p = wilcoxon(err1, err2);
                if (p > ALPHA) {
                    best_comp(i) = j+1; // +1 to convert from index to component number
                    break;
                }
            }
        }

        return best_comp;
    }

};


#endif
