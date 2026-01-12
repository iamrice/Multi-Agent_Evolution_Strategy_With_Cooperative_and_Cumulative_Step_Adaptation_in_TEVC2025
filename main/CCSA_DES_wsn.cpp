#include "./cmaes.h"
#include "../Benchmarks/WSN_Benchmarks.h"
#include <fstream>
#include <sys/time.h>
#include <time.h>
#include <unistd.h>
#include <regex>
#include <math.h>

using namespace std;

long getCurrentTime()
{
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return tv.tv_sec * 1000 + tv.tv_usec / 1000;
}
int main(int argc, char* argv[]){
    srand(getCurrentTime());

    int run_times = 25;
    string func_name,exid;
    string log_filename;
    int max_eva,run_idx=0;
    int interval = 1;
    double w1=0.9,w2=0.38,rate=0.01, local_rate = 0.001,angle = 90;
    double lr = 1;
    double fixed_step = 1e-6;

    func_name = argv[argc-1];
    
    exid = "ex02";
    max_eva = 1500000;
    log_filename = "../output/"+exid+"_CCSA-DES_"+func_name+"_"+to_string(run_idx);

    cout<<log_filename<<endl;
    
    WSN_Benchmarks* pFunc = new WSN_Benchmarks(func_name);
    int nodenum = pFunc->getNodeNum();
    pFunc->max_eva_times = max_eva * nodenum;

    MatrixXd Ajacent_mtx(nodenum,nodenum);
    double** W = pFunc->getNetworkGraph();
    for(int i=0;i<nodenum;i++){
        Ajacent_mtx.row(i) = VectorXd::Map(W[i],nodenum);
    }
    MatrixXd binary_net_mtx = (Ajacent_mtx.array()>0).cast<double>();
    MatrixXd average_mtx = MatrixXd::Ones(nodenum,nodenum) / nodenum;

    vector<vector<int>> neighbor_list;
    for(int i=0;i<nodenum;i++){
        vector<int> nei = pFunc->getNeighbors(i);
        nei.push_back(i);
        neighbor_list.push_back(nei);
    }

    int dim = pFunc->getDimension();
    // vector<cmaes*> optSet;
    vector<semi_cmaes*> optSet;
    for(int i=0;i<nodenum;i++){
        semi_cmaes* local_opt = new semi_cmaes(pFunc,i);
        optSet.push_back(local_opt);
    }
    double ub = pFunc->getMaxX();
    double lb = pFunc->getMinX();
    MatrixXd direction_momentum = MatrixXd::Zero(nodenum,dim);
    VectorXd gradient_momentum = VectorXd::Zero(dim);
    VectorXd grad_size = VectorXd::Zero(nodenum);
    VectorXd step_size = VectorXd::Zero(nodenum);

    MatrixXd system = MatrixXd::Random(nodenum,dim).array()*(ub-lb) + lb;
    MatrixXd grad_direction = MatrixXd::Zero(nodenum,dim);

    VectorXd sys_alpha = VectorXd::Ones(nodenum) * 1e-6;
    MatrixXd sys_ps = MatrixXd::Zero(nodenum,dim);


    int node_degree = 4;
    VectorXd rank_based_weights = VectorXd::Zero(node_degree);
    for (int i = 0; i < node_degree; ++i) {
        rank_based_weights(i) = log(node_degree + 0.5) - log(i + 1); // muXone array for weighted recombination
    }
    rank_based_weights /= rank_based_weights.sum(); 
    
    double commu_count=0;
    int iter=0;
    while(!pFunc->reachMaxEva()){
        for(int i=0;i<nodenum;i++){
            optSet[i]->sigma = sys_alpha(i);
            VectorXd origin_xmean = system.row(i);
            VectorXd xmean = system.row(i);

            optSet[i]->ps = VectorXd::Zero(dim);
            for(int g=0;g<interval;g++){
                VectorXd new_xmean = optSet[i]->step(xmean);
                optSet[i]->ps = (1 - optSet[i]->cs) * optSet[i]->ps + sqrt(optSet[i]->cs * (2 - optSet[i]->cs) * optSet[i]->mueff) * optSet[i]->invsqrtC * (new_xmean - xmean) / optSet[i]->sigma;                    
                if(g>0){
                    bool tau = false;
                    if((optSet[i]->ps.norm() / optSet[i]->chiN - 1) * (direction_momentum.row(i).norm() - 1) > 0){
                        tau = true;
                    }
                    optSet[i]->sigma = optSet[i]->sigma * exp(local_rate * (optSet[i]->ps.norm() / optSet[i]->chiN - 1) * tau);                        
                }
                xmean = new_xmean;
            }
            
            step_size(i) = (xmean - origin_xmean).norm();
            grad_size(i) = (pFunc->local_eva(origin_xmean.data(),i) - pFunc->local_eva(xmean.data(),i)) / step_size(i);
            grad_direction.row(i) = (xmean-origin_xmean).array() / step_size(i);
            system.row(i) = xmean;
            sys_alpha(i) = optSet[i]->sigma;
        }      

        MatrixXd merged_direction = MatrixXd::Zero(nodenum,dim);
        for(int i=0;i<nodenum;i++){
            for(int j=0;j<nodenum;j++){
                if(Ajacent_mtx(i,j) > 0){
                    merged_direction.row(i) += grad_size(j) * grad_direction.row(j);
                }
            }
            merged_direction.row(i) = merged_direction.row(i) / merged_direction.row(i).norm();
        }
        Eigen::Array<bool, Eigen::Dynamic, Eigen::Dynamic> nan_mask = merged_direction.array().isNaN();
        Eigen::Array<bool, Eigen::Dynamic, Eigen::Dynamic> inf_mask = merged_direction.array().isInf();
        if (!nan_mask.any() && !inf_mask.any()) {
            angle = 90 - 90 * (pFunc->eva_count*1.0 / pFunc->max_eva_times);
            w2 = (-2*w1*cos(angle/180*M_PI)+sqrt(4*w1*w1*cos(angle/180*M_PI)*cos(angle/180*M_PI)-4*w1*w1+4))/2;    

            direction_momentum = Ajacent_mtx * direction_momentum;
            if(iter==0){
                direction_momentum = merged_direction; 
            }else{
                direction_momentum = w1 * direction_momentum + w2 * merged_direction;
            }
            for(int i=0;i<nodenum;i++){
                sys_alpha(i) *= exp(rate*(direction_momentum.row(i).norm() - 1));
            }
        }
        system = Ajacent_mtx * system; 
        
        commu_count += system.size() + direction_momentum.size() + grad_direction.size();
        
        if (iter % 100 == 0){
            VectorXd xmean = system.colwise().mean();//v1
            double curr_fit = pFunc->global_eva(xmean.data());
            double disagreement = (system.rowwise() - xmean.transpose()).rowwise().norm().mean();
            fstream ft(log_filename, ios::app | ios::out);
            ft << curr_fit<<" "<<iter <<" " <<pFunc->eva_count*1.0/nodenum <<" "<<disagreement<<" "<<sys_alpha.mean()<<" "<< direction_momentum.rowwise().norm().mean()<<" "<<commu_count<<endl;
            ft.close();
            if(disagreement < 1e-10)
                break;
        }

        iter++;
    }
    VectorXd xmean = system.colwise().mean();
    fstream ft(log_filename+"_res", ios::app | ios::out);
    ft<<xmean.transpose()<<endl;
    ft.close();
}