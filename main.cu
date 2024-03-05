#include <iostream>
#include <cuda_runtime.h>

#include "Cluster.h"
#include "Point.h"
#include <cmath>
#include <vector>
#include <string>
#include <fstream>
#include <chrono>
#include <filesystem>
#include "SequentialKMeans.h"
#include <algorithm>


#define N 400000
#define TPB 128
static std::vector<int> used_pointIds;
using namespace std;

__host__ void readPoints(double *h_xval, double *h_yval, int* h_clusters) {

    std::string tmp = "";
    string line;
    std::cout << std::filesystem::current_path().string() << "Ciao" << std::endl;
    ifstream infile("../cmake-build-debug/input1.txt");
    if (!infile.is_open()) {
        cout << "Error: Failed to open file." << endl;
        return;
    }
    int j = 0;
    while (getline(infile, line)) {
        for (int i = 0; i < static_cast<int>(line.length()); i++) {
            if ((48 <= static_cast<int>(line[i]) && static_cast<int>(line[i]) <= 57) || line[i] == '.' ||
                line[i] == '+' || line[i] == '-' || line[i] == 'e') {
                tmp += line[i];

            } else if (!tmp.empty()) {


                h_xval[i] = std::stod(tmp);
                //xval.push_back(std::stod(tmp));
                tmp = "";
            }
        }
        if (!tmp.empty()) {

           h_yval[j] = std::stod(tmp);
            //yval.push_back(std::stod(tmp));
            tmp = "";
        }
       h_clusters[j] = -1; //inizializzo anche il vettore clusters a -1
       j++;
    }
    //std::cout << "dim clusters: " << h_clusters << std::endl;

    cout << "\nDataset fetched!" << endl
         << endl;
    //std::cout << "Punti totali letti : " << xval.size() << std::endl;
}



static bool first = true;
std::vector<int>indexGenerator(int K, int total_points){
    //mem punti gia usati per init cluster
    if (first) {
        for (int i = 0; i < K; i++) {
            while (true) {
                int index = rand() % total_points;

                //check se index gia usato per un altro cluster

                if (std::find(used_pointIds.begin(), used_pointIds.end(), index) ==
                    used_pointIds.end()) {
                    used_pointIds.push_back(index);
                    //creo un cluster con avete centroide il punto attuale
                    break;
                }
            }
        }
        first = false;
    }
    return used_pointIds;
}
__global__ void clusterAssignment(const double *d_xval, const double *d_yval, int *d_clusterval, const double *d_centroidX, const double *d_centroidY, int K){
    const int idx = blockIdx.x*blockDim.x + threadIdx.x;
    if (idx >= N) return;
    float min_dist = INFINITY;
    int closest_centroid = 0;

    double dist;
    for(int i = 0; i < K; i++)
    {
        double sum = 0.0;
        sum += pow(d_centroidX[i] - d_xval[idx], 2.0);
        sum += pow(d_centroidY[i] - d_yval[idx], 2.0);
        dist = sqrt(sum);
        if(dist < min_dist){
            min_dist = dist;
            closest_centroid = i;
        }
    }
    //assegno id-cluster al thread corrente
    d_clusterval[idx] = closest_centroid;
}

__global__ void clusterPointsSum(double* d_xval, double* d_yval, int* d_clusterVal, double* d_clusterSumX, double* d_clusterSumY, int* d_clusterSize){
    //indice del thread a livello grid

    const int idx = blockIdx.x * blockDim.x + threadIdx.x;

    if (idx >= N) return;

    int clusterId = d_clusterVal[idx];
    //sommo tutti i punti appartenenti ai clusters
    atomicAdd(&(d_clusterSumX[clusterId]), d_xval[idx]);
    atomicAdd(&(d_clusterSumY[clusterId]), d_yval[idx]);
    atomicAdd(&(d_clusterSize[clusterId]), 1);

}


//////
//////Definizione classe
//////
class ParallelKMeans{
private:
    int K, iters, dimensions{}, total_points{};
    std::vector<Cluster> clusters;
    std::string input_dir;
    std::string output_dir;
    Point all_points;


    int getNearestClusterId(int pointIdx);
    int getNearestClusterId_Parallel(Point point);

public:
    ParallelKMeans(int K, int iterations, std::string output_dir, std::string input_dir);
    void run();
    void run_parallel(std::vector<Point> &all_points);
    void run_parallel2(std::vector<Point> &all_points);
    ~ParallelKMeans() {
        // Releases all the remaining resources allocated on the GPU
        cudaDeviceReset();
    }
};
/////
/////Implementazione costruttore
/////
ParallelKMeans::ParallelKMeans(int K, int iterations, std::string output_dir, std::string input_dir):all_points(input_dir) {
    this->K = K;
    this->iters = iterations;
    this->output_dir = output_dir;
    //this->all_points = Point(input_dir);
}
void ParallelKMeans::run() {
    total_points = all_points.getDimensions();  //num totale di punti

    //alloco memoria Host

    double *h_xval = (double *) malloc(N * sizeof(double));
    double *h_yval = (double *) malloc(N * sizeof(double));
    int *h_clusterVal = (int *) malloc(N * sizeof(int));
    //double*h_centroids =
    double *h_centroidX = (double *) malloc(K * sizeof(double));
    double *h_centroidY = (double *) malloc(K * sizeof(double));
    int *h_clusterSize = (int *) malloc(K * sizeof(int));
    double *h_clusterSumX = (double *) malloc(K * sizeof(double));
    double *h_clusterSumY = (double *) malloc(K * sizeof(double));

    //alloco memoria Device
    double *d_xval;
    double *d_yval;
    //double*h_centroids =
    int *d_clusterVal;
    double *d_centroidX;
    double *d_centroidY;
    int *d_clusterSize;
    double *d_clusterSumX;
    double *d_clusterSumY;

    cudaMalloc(&d_xval, total_points * sizeof(double));
    cudaMalloc(&d_yval, total_points * sizeof(double));
    cudaMalloc(&d_centroidX, K * sizeof(double));
    cudaMalloc(&d_centroidY, K * sizeof(double));
    cudaMalloc(&d_clusterVal, N * sizeof(int));
    cudaMalloc(&d_clusterSize, K * sizeof(int));
    cudaMalloc(&d_clusterSumX, K * sizeof(double));
    cudaMalloc(&d_clusterSumY, K * sizeof(double));

    readPoints(h_xval, h_yval, h_clusterVal);       //leggo i punti dataset

    dimensions = 2;   //dimensione punto

    //std::cout<<"dimensione all_points: "<<total_points<<" dimensione clusters: "<<all_points.getDimClusters()<<" punto 1: "<<all_points.getYval(1)<< std::endl;
    // Inizializzo Clusters
    std::vector<int> used_pointIds;  //mem punti gia usati per init cluster

    bool exit = false;
    double x, y;

    used_pointIds = indexGenerator(K, N);
    for (int i = 0; i < K; i++) {

        h_centroidX[i] = h_xval[used_pointIds[i]];
        h_centroidY[i] = h_yval[used_pointIds[i]];
        //x = all_points.getXval(used_pointIds[i]);
        //y = all_points.getYval(used_pointIds[i]);
        h_clusterVal[used_pointIds[i]] = i;
        h_clusterSize[i] = 0;
        //all_points.setCluster(used_pointIds[i], i);
        //Cluster cluster(i,used_pointIds[i], x, y);  //creo un cluster avente centroide il punto attuale
        //clusters.push_back(cluster);
    }
    std::cout << "Clusters Inizializzati = " << std::endl
              << std::endl;

    std::cout << "Eseguo Clustering K-Means..." << std::endl;

    int iter = 1;

    //copia host -> device

    cudaMemcpy(d_xval, h_xval, N * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_yval, h_yval, N * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_centroidX, h_centroidX, K * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_centroidY, h_centroidY, K * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_clusterSize, h_clusterSize, K * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_clusterSumX, h_centroidX, K * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_clusterSumY, h_centroidY, K * sizeof(double), cudaMemcpyHostToDevice);


    while (iter <= iters) {
        std::cout << "Iter - " << iter << "/" << iters << std::endl;
        bool done = true;
        clusterAssignment<<<(N + TPB - 1) / TPB, TPB>>>(d_xval, d_yval, d_clusterVal, d_centroidX, d_centroidY, K);

        //setto a 0 la sommma dei punti appartenenti ai cluster, per iniziare aggiornamento centroide
        cudaMemset(d_clusterSumX, 0.0, K * sizeof(double));
        cudaMemset(d_clusterSumY, 0.0, K * sizeof(double));

        clusterPointsSum<<<(N + TPB - 1) / TPB, TPB>>>(d_xval, d_yval, d_clusterVal, d_clusterSumX, d_clusterSumY,
                                                       d_clusterSize);

        //copia variabili aggiornate da device -> host

        cudaMemcpy(h_clusterSumX, d_clusterSumX, K * sizeof(double), cudaMemcpyDeviceToHost);
        cudaMemcpy(h_clusterSumY, d_clusterSumY, K * sizeof(double), cudaMemcpyDeviceToHost);
        cudaMemcpy(h_clusterSize, d_clusterSize, K * sizeof(int), cudaMemcpyDeviceToHost);

        cudaMemset(d_clusterSize, 0, K * sizeof(int));

        //Ricalcolo nuovi centroidi per ogni cluster
        for (int i = 0; i < K; i++) {
            h_centroidX[i] = h_clusterSumX[i] / h_clusterSize[i];
            h_centroidY[i] = h_clusterSumY[i] / h_clusterSize[i];
        }

        cudaMemcpy(d_centroidX, h_centroidX, K * sizeof(double), cudaMemcpyHostToDevice);
        cudaMemcpy(d_centroidY, h_centroidY, K * sizeof(double), cudaMemcpyHostToDevice);

        iter++;

    }
    cudaMemcpy(h_clusterVal, d_clusterVal, N * sizeof(int), cudaMemcpyDeviceToHost);


    ///scrivo i risultati clustering su file
    //scrittura centroidi
    std::ofstream outfile;
    std::cout << std::filesystem::current_path().string() << std::endl;
    std::cout << output_dir + "/" + "clusters.txt" << std::endl;
    outfile.open(output_dir + "/" + "clusters.txt");
    if (outfile.is_open()) {
        for (int i = 0; i < K; i++) {
            //std::cout <<  i << " cluster contiene: "<< clusters[i].getSize() <<std::endl;
            std::cout << "Cluster " << i << " centroid : ";

            std::cout << h_centroidX[i] << " ";    // Output console
            std::cout << h_centroidY[i] << " ";
            outfile << h_centroidX[i] << " " << h_centroidY[i]; // Output file
        }
        std::cout << std::endl;
        outfile << std::endl;
        outfile.close();
    } else {
        std::cout << "Error: Unable to write to clusters.txt";
    }
    //scrittura punti e cluster di appartenenza
    outfile.open(output_dir + "/" + "clustering.txt");
    for (int i = 0; i < N; i++) {
        // indice cluster è ID-1
        if (outfile.is_open()) {

            outfile << h_xval[i] << " " <<  h_yval[i] << " "
                    << h_clusterVal[i]  ;// Output to file
        }

        outfile << std::endl;

    }
    outfile.close();
}



double averageParallelExecutions(int K, int iters, std::string output_dir, std::string input_dir)
{
    double mediaS, mediaP;
    double sum;

    for(int i  = 0; i < 2; i++ ) {
        auto start = std::chrono::high_resolution_clock::now();
        ParallelKMeans kmeans(K, iters, output_dir, input_dir);
        //kmeans.run_parallel2(all_points);
        kmeans.run();
        //double end = omp_get_wtime( );
        auto end = std::chrono::high_resolution_clock::now();

        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
        std::cout << "Tempo di esecuzione parallela: " << duration.count() << " millisecondi" << std::endl;
        sum += static_cast<double>(duration.count());
        std::cout << "<<------------------------------>>" << std::endl;

    }
    mediaP = static_cast<double>(sum)/2;
    return mediaP;
}

double averageSeqExecutions(int K, int iters, std::string output_dir, std::string input_dir)
{
    double mediaS;
    double sum;

    for(int i  = 0; i < 2; i++ ) {
        auto start = std::chrono::high_resolution_clock::now();
        SequentialKMeans kmeans(K, iters, output_dir, input_dir);
        kmeans.run();
        //double end = omp_get_wtime( );
        auto end = std::chrono::high_resolution_clock::now();

        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
        std::cout << "Tempo di esecuzione Sequenziale: " << duration.count() << " millisecondi" << std::endl;
        sum += static_cast<double>(duration.count());
        std::cout << "<<------------------------------>>" << std::endl;

    }
    mediaS = static_cast<double>(sum)/2;
    return mediaS;
}

int main() {
    std::string output_dir = "../cmake-build-debug/cluster_details";   //dir output
    int K = 4;                               //numero cluster
    std::string input_dir= "input1.txt";

    // Avvio il clustering
    int iters = 100;


    //auto mediaS = averageSeqExecutions(K, iters, output_dir, input_dir);
    auto mediaP = averageParallelExecutions(K, iters, output_dir, input_dir);
    //std::cout << "Media esecuzione Sequenziale : " << mediaS << std::endl;
    std::cout << "Media esecuzione Parallela : " << mediaP << std::endl;

    //double speedup = static_cast<double>(mediaS) / static_cast<double>(mediaP);
    //std::cout << "Speedup: " << speedup << std::endl;
    return 0;
}
