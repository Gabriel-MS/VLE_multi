#ifndef GIBBS_H_INCLUDED
#define GIBBS_H_INCLUDED

#define EIGEN_NO_DEBUG
#include <Eigen/Dense>
#include <unsupported/Eigen/MatrixFunctions>

using namespace std;
using namespace Eigen;

//Arquivo dedicado a conter funções que calculem a energia livre de Gibbs em excesso e o coeficiente de atividade

double gibbs_excess_NRTL(int nc, double T, double R, VectorXd x, MatrixXd A, MatrixXd alfa_NRTL)
{
    MatrixXd LAMBDA(nc,nc), ln_G(nc,nc), G(nc,nc), E(nc,nc), Gt_x(nc,nc);
    double G_ex_RT, G_ex;

    ln_G = (-alfa_NRTL/T).cwiseProduct(A);
    G = ln_G.array().exp();
    LAMBDA = (A/T).cwiseProduct(G);
    Gt_x = G.transpose()*x;
    E = LAMBDA*(Gt_x.asDiagonal().inverse());
    G_ex_RT = x.transpose()*E*x;
    G_ex = G_ex_RT*R*T;

    return G_ex;
}

VectorXd gamma_NRTL(int nc, double T, double R, VectorXd x, MatrixXd A, MatrixXd alfa_NRTL)
{
    MatrixXd LAMBDA(nc,nc), ln_G(nc,nc), G(nc,nc), E(nc,nc), L(nc,nc), Gt_x(nc,nc);
    VectorXd ln_gamma(nc), gamma(nc);
    double G_ex_RT, G_ex;

    ln_G = (-alfa_NRTL/T).cwiseProduct(A);
    G = ln_G.array().exp();
    LAMBDA = (A/T).cwiseProduct(G);
    Gt_x = G.transpose()*x;
    E = LAMBDA*(Gt_x.asDiagonal().inverse());
    L = G*(Gt_x.asDiagonal().inverse());
    ln_gamma = (E + E.transpose()-L*x.asDiagonal()*E.transpose())*x;
    gamma = ln_gamma.array().exp();

    return gamma;
}

VectorXd gamma_UNIQUAC(int nc, VectorXd q, VectorXd q_, VectorXd r, VectorXd x, MatrixXd Aij)
{
    VectorXd xr(nc), xq(nc), xq_(nc), PHI(nc), theta(nc), theta_(nc), SUMxr(nc), SUMxq(nc), SUMxq_(nc);
    VectorXd PHI_x(nc), ln_PHI_x(nc), theta_PHI(nc), ln_theta_PHI(nc), qln_theta_PHI(nc), L(nc), xL(nc), SUMxL(nc);
    VectorXd PHI_x_sumxL(nc), theta_A(nc), SUMtheta_A(nc), ln_SUMtheta_A(nc), q_lnt(nc), ln_gamma(nc), gamma(nc);
    RowVectorXd thetat(nc);
    MatrixXd ln_PHI_xdiag(nc,nc), PHI_xdiag(nc,nc);


    VectorXd uni(nc);
    int i;

    for(i=0; i<nc; i++)
    {
        uni(i) = 1;
    }

    //Calculating parameters
    xr = x.cwiseProduct(r);
    xq = x.cwiseProduct(q);
    xq_ = x.cwiseProduct(q_);

    PHI = xr/(xr.sum());
    theta = xq/(xq.sum());
    theta_ = xq_/(xq_.sum());

    //calculating gamma
    //1º termo
    PHI_x = PHI.cwiseQuotient(x);
    PHI_xdiag = PHI_x.asDiagonal();
    ln_PHI_xdiag = PHI_xdiag.array().log();
    ln_PHI_x = ln_PHI_xdiag.diagonal();

    //2º termo
    MatrixXd theta_PHIdiag(nc,nc), ln_theta_PHIdiag(nc,nc);
    theta_PHI = theta.cwiseQuotient(PHI);
    theta_PHIdiag = theta_PHI.asDiagonal();
    ln_theta_PHIdiag = theta_PHIdiag.array().log();
    ln_theta_PHI = ln_theta_PHIdiag.diagonal();
    qln_theta_PHI = 5*(q.cwiseProduct(ln_theta_PHI));

    //3º termo
    L = 5*(r-q)-(r-uni);

    //4º termo
    xL = x.cwiseProduct(L);
    PHI_x_sumxL = PHI_x*(xL.sum());

    //5ºtermo
    MatrixXd theta_Adiag(nc,nc), ln_theta_Adiag(nc,nc), theta_m(nc,nc), Aijt(nc,nc), theta_mt(nc,nc), theta_Amt(nc,nc), theta_Am(nc,nc);
    MatrixXd theta_Af(nc,nc), theta_Bm(nc,nc), theta_B(nc,nc);
    VectorXd ln_theta_A(nc), theta_At(nc);
    int j;

    thetat = theta_.transpose();
    for(j=0; j<nc; j++)
    {
        for(i=0; i<nc; i++)
        {
            theta_m(i,j) = thetat[i];
        }
    }
    theta_mt = theta_m.transpose();
    Aijt = Aij.transpose();
    theta_Am = theta_m.cwiseProduct(Aij);
    theta_Bm = theta_m.cwiseProduct(Aijt);
    theta_Amt = theta_Am.transpose();

    RowVectorXd unit(nc);
    for(i=0; i<nc; i++)
    {
        unit(i) = 1;
    }

    theta_A = unit*theta_Am;
    theta_At = theta_A.transpose();

    for(j=0; j<nc; j++)
    {
        for(i=0; i<nc; i++)
        {
            theta_Af(i,j) = theta_At[i];
        }
    }
    theta_B = theta_Bm.cwiseQuotient(theta_Af);

    theta_A = thetat*Aij;
    theta_Adiag = theta_A.asDiagonal();
    ln_theta_Adiag = theta_Adiag.array().log();
    ln_theta_A = ln_theta_Adiag.diagonal();
    q_lnt = q_.cwiseProduct(ln_theta_A);

    //6º termo
    q_ = q_;

    //7º termo
    VectorXd qtheta_B(nc,nc);
    VectorXd theta_Bf(nc);

    theta_Bf = unit*theta_B;
    qtheta_B = q.cwiseProduct(theta_Bf);


    ln_gamma = ln_PHI_x + qln_theta_PHI + L - PHI_x_sumxL - q_lnt + q_ - qtheta_B;

    MatrixXd ln_gammadiag(nc,nc), gammadiag(nc,nc);

    ln_gammadiag = ln_gamma.asDiagonal();
    gammadiag = ln_gammadiag.array().exp();
    gamma = gammadiag.diagonal();


    return gamma;
}

double gibbs_excess_UNIQUAC(int nc, double T, VectorXd q, VectorXd r, VectorXd x, MatrixXd A, double R)
{
    VectorXd epsilon(nc), m(nc), phi(nc), theta(nc), log_q(nc), log_phi(nc), log_theta(nc), log_lambda_x(nc), lambda_x(nc);
    MatrixXd A_T(nc,nc), lambda(nc,nc);
    double G_ex_RT, G_ex;


    log_q = q.array().log();


    epsilon = q + q.asDiagonal()*log_q;
    m = 1 - q.array()*5;
    phi = (pow((r.transpose()*x),-1))*r;
    log_phi = phi.array().log();
    theta = (pow((q.transpose()*x),-1))*q;
    log_theta = theta.array().log();
    A_T = A/T;
    lambda = q.asDiagonal()*((-A_T).exp());
    lambda_x = lambda.transpose()*x;
    log_lambda_x = lambda_x.array().log();
    G_ex_RT = x.transpose()*(m.asDiagonal()*(log_phi) + 4*(q.asDiagonal()*(log_theta)) + epsilon - q.asDiagonal()*(log_lambda_x) - q);
    G_ex = G_ex_RT*(R*T);

    return G_ex;
}



#endif // GIBBS_H_INCLUDED
