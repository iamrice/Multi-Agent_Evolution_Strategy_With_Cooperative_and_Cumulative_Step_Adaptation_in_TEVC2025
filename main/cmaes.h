#ifndef CMAES_H
#define CMAES_H

#include <iostream>
#include <Eigen/Dense>
#include <vector>
#include <cmath>
#include <random>
#include "../Benchmarks/Benchmarks.h"
#include "../Benchmarks/WSN_Benchmarks.h"


using namespace Eigen;
static std::default_random_engine e(time(0));
static std::normal_distribution<double> n(0,1);

class semi_cmaes{
public:
    Benchmarks* pFunc;
    WSN_Benchmarks* wsn_func;
    int node_index;
    int N;
    // VectorXd xmean;
    double sigma;
    double stopfitness;
    double stopeval;

    // Strategy parameter setting: Selection
    int lambda;
    double mu;
    VectorXd weights;
    double mueff;

    // Strategy parameter setting: Adaptation
    double cc,cs,c1,cmu,damps;

    // Initialize dynamic (internal) strategy parameters and constants
    VectorXd pc,ps;
    MatrixXd B;
    VectorXd D;
    MatrixXd C;
    MatrixXd invsqrtC;
    int eigeneval;
    double chiN;
    double range;

    int counteval;


    semi_cmaes(WSN_Benchmarks* func,int node_index=-1){
        pFunc = nullptr;
        wsn_func = func;
        this->node_index = node_index;
        N=wsn_func->getDimension();
        range = wsn_func->getMaxX()-wsn_func->getMinX();
        para_init();
    }
    semi_cmaes(Benchmarks* func,int node_index=-1){
        pFunc = func;
        wsn_func = nullptr;
        this->node_index = node_index;
        N=pFunc->getDimension();
        range = pFunc->getMaxX()-pFunc->getMinX();
        para_init();
    }

    void para_init(){

        // xmean =  VectorXd::Random(N); // objective variables initial point
        sigma =  0.5; // coordinate wise standard deviation (step size)
        stopfitness =  1e-10; // stop if fitness < stopfitness (minimization)
        stopeval =   N * N; // stop after stopeval number of function evaluations

        // Strategy parameter setting: Selection
        lambda =  4 + floor(3 * log(N)); // population size, offspring number
        mu =  lambda / 2.0; // number of parents/points for recombination
        int weight_length = floor(mu);
        weights = VectorXd::Zero(weight_length);
        for (int i = 0; i < weight_length; ++i) {
            weights(i) = log(mu + 0.5) - log(i + 1); // muXone array for weighted recombination
        }
        mu = floor(mu);
        weights = weights / weights.sum(); // normalize recombination weights array
        mueff =  weights.sum() * weights.sum() / weights.cwiseProduct(weights).sum(); // variance-effectiveness of sum w_i x_i

        // Strategy parameter setting: Adaptation
        cc =  (4 + mueff / N) / (N + 4 + 2 * mueff / N); // time constant for cumulation for C
        cs =  (mueff + 2) / (N + mueff + 5); // t-const for cumulation for sigma control
        c1 =  2 / ((N + 1.3) * (N + 1.3) + mueff); // learning rate for rank-one update of C
        cmu =  std::min(1 - c1, 2 * (mueff - 2 + 1 / mueff) / ((N + 2) * (N + 2) + mueff)); // and for rank-mu update
        damps =  1 + 2 * std::max(0.0, sqrt((mueff - 1) / (N + 1)) - 1) + cs; // damping for sigma

        // Initialize dynamic (internal) strategy parameters and constants
        pc =  VectorXd::Zero(N); // evolution paths for C
        ps =  VectorXd::Zero(N); // evolution paths for sigma
        B =  MatrixXd::Identity(N, N); // B defines the coordinate system
        D =  VectorXd::Ones(N)*100; // diagonal D defines the scaling
        C =  B * D.array().square().matrix().asDiagonal() * B.transpose(); // covariance matrix C
        invsqrtC =  B * D.array().inverse().matrix().asDiagonal() * B.transpose(); // C^-1/2
        eigeneval =  0; // track update of B and D
        chiN =  sqrt(N) * (1 - 1.0 / (4 * N) + 1.0 / (21 * N * N)); // expectation of ||N(0,I)|| == norm(randn(N,1))

        counteval = 0;
    }

    VectorXd step(VectorXd xmean){
        //std::cout<<"flag06"<<std::endl;
        MatrixXd arx(N, lambda);
        VectorXd arfitness(lambda);
        for (int k = 0; k < lambda; ++k) {
            VectorXd rand_vec = VectorXd::Zero(N).unaryExpr([](double dummy){return n(e);});
            // arx.col(k) = xmean + sigma * (B * (D.asDiagonal() * rand_vec)); 
            arx.col(k) = xmean + sigma * rand_vec;
            // arfitness(k) = sphereFunction(arx.col(k));
            VectorXd x_to_eva = arx.col(k);
            if(pFunc!=nullptr){
                if(node_index == -1)
                    arfitness(k) = pFunc->global_eva(x_to_eva.data());
                else
                    arfitness(k) = pFunc->local_eva(x_to_eva.data(),node_index);
            }
            else{
                if(node_index == -1)
                    arfitness(k) = wsn_func->global_eva(x_to_eva.data());
                else
                    arfitness(k) = wsn_func->local_eva(x_to_eva.data(),node_index);
            }
        }
        //std::cout<<"flag07"<<std::endl;

        // Sort by fitness and obtain sorted indices
        std::vector<std::pair<double, int>> sorted_fitness_indices;
        sorted_fitness_indices.reserve(lambda);
        for (int i = 0; i < lambda; ++i) {
            sorted_fitness_indices.emplace_back(arfitness(i), i);
        }
        std::sort(sorted_fitness_indices.begin(), sorted_fitness_indices.end());

        //std::cout<<"flag08"<<std::endl;
        // Extract indices in sorted order
        VectorXi arindex(lambda);
        MatrixXd sorted_arx = arx.array();
        for (int i = 0; i < lambda; ++i) {
            arindex(i) = sorted_fitness_indices[i].second;
            sorted_arx.col(i) = arx.col(arindex(i));
        }

        //std::cout<<"flag09"<<std::endl;
        VectorXd xold = xmean.array();
        xmean.setZero();
        for(int i=0;i<mu;i++){
            xmean += arx.col(arindex(i)) * weights(i);
        }
        // xmean = arx * arindex.head(mu).cast<double>().cwiseProduct(weights.transpose()); 

        //std::cout<<"flag10"<<std::endl;
        // Cumulation: Update evolution paths
        // VectorXd ps_temp = (1 - cs) * ps + sqrt(cs * (2 - cs) * mueff) * invsqrtC * (xmean - xold) / sigma;
        // ps = ps_temp;
        // // Adapt step size sigma
        // sigma = sigma * exp((cs / damps) * (ps.norm() / chiN - 1));

        return xmean;
    }
};

#endif