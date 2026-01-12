#include <sstream>
#include <fstream>
#include <iostream>
#include<cstring>
#include<cstdio>
#include<cmath>
#include<algorithm>
#include "../Benchmarks/WSN_Benchmarks.h"
#include "../Benchmarks/Benchmarks.h"
using namespace std;
using json = nlohmann::json;

double measurement_2d(double x0, double y0, double x1, double y1, string type){
    double res = 0;
    if(type == "distance"){
        res = sqrt(pow(x0 - x1,2) + pow(y0  - y1,2) );
    }else if (type == "TDOA"){
        double send_source[3] = {85.3031,   62.2055 ,  35.0952};
        res = sqrt(pow(x0 - send_source[0],2) + pow(y0 - send_source[1],2)) + sqrt(pow(x1 - x0,2) + pow(y1 - y0,2));
    }else if (type == "RSS"){
        double dis = sqrt(pow(x0 - x1,2) + pow(y0  - y1,2) );
        res = 100-10*2*log10(dis);
    }
    return res;
}

double measurement_3d(double x0, double y0, double z0, double x1, double y1, double z1, string type){
    double res = 0;
    if(type == "distance"){
        res = sqrt(pow(x0 - x1,2) + pow(y0  - y1,2) + pow(z0 - z1,2));
    }else if (type == "TDOA"){
        double send_source[2] = {85.3031,   62.2055};
        res = sqrt(pow(x0 - send_source[0],2) + pow(y0 - send_source[1],2) + pow(z0 - send_source[2],2)) + sqrt(pow(x1 - x0,2) + pow(y1 - y0,2) + pow(z1 - z0,2));
    }else if (type == "RSS"){
        double dis = sqrt(pow(x0 - x1,2) + pow(y0  - y1,2) + pow(z0 - z1,2));
        res = 100-10*2*log10(dis);
    }
    return res;
}

WSN_Benchmarks::WSN_Benchmarks(string ID,int max_eva_times){    
    string config_path = "../Benchmarks/default_config.json";
    string data_path = "../Benchmarks/data/";
    ifstream file_config(config_path);
    json config;
    file_config >> config;
    file_config.close();
    this->func_config = config["benchmarks"][ID];
    this->funcID = ID;
    this->node_num = func_config["node_num"];

    this->max_eva_times = max_eva_times;
    this->eva_count = 0;
    this->reach_max_eva_times = false;

    measurement_type = func_config["measurement"];
    target_num = func_config["target_num"];
    source_num = func_config["node_num"];
    coordinate_dim = func_config["coordinate"];
    func_config["dimension"] = target_num * coordinate_dim;
    W=read_data<double>(data_path+"W_"+to_string(source_num)+"_WSN_location");
    load_and_init(data_path+"target_WSN_location", data_path+"source_WSN_location");
}


void WSN_Benchmarks::load_and_init(string target_fname, string source_fname){
    source = read_data<double>(source_fname);
    target = read_data<double>(target_fname);
    if(target == nullptr||source == nullptr){
        cout<<target_fname<<" "<<source_fname<<endl;
        return;
    }
    bool if_noisy = func_config["noisy"];
    double degree=0;
    if(if_noisy == true){
        degree = func_config["noisy_degree"];
    }

    std::random_device rd;
    std::mt19937 gen(rd());
    std::normal_distribution<> d(0, degree);

    noisy_dis = new double*[source_num];
    for(int i=0;i<source_num;i++){
        noisy_dis[i] = new double[target_num];
        for(int j=0;j<target_num;j++){                
            if(coordinate_dim == 3){
                noisy_dis[i][j] = measurement_3d(target[j][0],target[j][1],target[j][2],source[i][0],source[i][1],source[i][2],func_config["measurement"]);
            }else if(coordinate_dim == 2){
                noisy_dis[i][j] = measurement_2d(target[j][0],target[j][1],source[i][0],source[i][1],func_config["measurement"]); 
            }                
            if(if_noisy == true){
                noisy_dis[i][j] += d(gen);
            }
        }
    }
    
}

WSN_Benchmarks::~WSN_Benchmarks() {

}

bool WSN_Benchmarks::reachMaxEva(){
    if(eva_count>=max_eva_times){
        if(!reach_max_eva_times){
            // cout<<"The time of evaluation has reached the maximum bound. Later evaluation results would not be recorded.\n"; 
            reach_max_eva_times = true;
        }
        return true;
    }
    return false;
}


double WSN_Benchmarks::local_eva(double* x, int groupIndex){
    double res = 0;
    for (size_t t = 0; t < target_num; t++){
        double est_dis;
        if(coordinate_dim == 3)
            est_dis = measurement_3d(x[t*3],x[t*3+1],x[t*3+2],source[groupIndex][0],source[groupIndex][1],source[groupIndex][2],measurement_type);
        else if(coordinate_dim == 2)
            est_dis = measurement_2d(x[t*2],x[t*2+1],source[groupIndex][0],source[groupIndex][1],measurement_type);
        res += pow(noisy_dis[groupIndex][t] - est_dis,2);
    }
    this->eva_count += 1;
    return res;
}


double WSN_Benchmarks::global_eva(double* x) {
	double res = 0;
    for(int i=0;i<node_num;i++){
        double r = local_eva(x,i);
        res += r;
    }
    res/=node_num;
	return res;
}

double WSN_Benchmarks::getMinX() {
	return func_config["lower_bound"];
}
double WSN_Benchmarks::getMaxX() {
	return func_config["upper_bound"];
}
int WSN_Benchmarks::getNodeNum() {
	return func_config["node_num"];
}
int WSN_Benchmarks::getDimension(){
    return func_config["dimension"];
}
double** WSN_Benchmarks::getNetworkGraph(){
    if(W==nullptr){
        string data_path = "../Benchmarks/data/";
        W=read_data<double>(data_path+"W_"+to_string(this->node_num)+"n");
    }
    return W;
}


vector<int> WSN_Benchmarks::getNeighbors(int groupIndex){
    vector<int> groups;
    for(int i=0;i<node_num;i++){
        if(i == groupIndex)
            continue;
        if(W[groupIndex][i]>0)
            groups.push_back(i);
    }    
    return groups;
}

template<typename T>
T** WSN_Benchmarks::read_data(string fileName) {
	// cout << fileName<<endl;
    T** res = nullptr;
    ifstream file(fileName);
	if (file.is_open()) {
		// cout << " is opened;\n";
        vector<vector<T>> data;
		string line;
        while(getline(file,line)){
    // cout<<"1"<<endl;
            stringstream ss(line);
            vector<T> rowData;
            T x;
            while(ss>>x){
                rowData.push_back(x);
            }
            data.push_back(rowData);
        }
		file.close();
    // cout<<"1"<<endl;

        res = new T*[data.size()];
        int count = 0;
        for(auto rowData:data){
    // cout<<"1"<<endl;
            // cout<<rowData.size()<<endl;
            T* row = new T[rowData.size()];
            memcpy(row,&rowData[0],rowData.size()*sizeof(T));
            res[count]=row;
            count++;
        }
	}
	else {
		// cout << " can not be opened;\n";
	}
	return res;
}

