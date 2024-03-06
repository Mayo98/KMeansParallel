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

__host__ void readPoints(float *h_xval, float *h_yval, int* h_clusters) {

    std::string tmp = "";
    string line;
    //std::cout << std::filesystem::current_path().string()  << std::endl;
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


                h_xval[j] = std::stod(tmp);

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
__global__ void clusterAssignment(const float *d_xval, const float *d_yval, int *d_clusterval, const float *d_centroidX, const float *d_centroidY, int K, bool *d_done){
    const int idx = blockIdx.x*blockDim.x + threadIdx.x;
    *d_done = true;
    if (idx >= N) return;
    float min_dist = INFINITY;
    int closest_centroid = 0;

    float dist;
    for(int i = 0; i < K; i++)
    {
        float sum = 0.0;
        sum += pow(d_centroidX[i] - d_xval[idx], 2.0);
        sum += pow(d_centroidY[i] - d_yval[idx], 2.0);
        dist = sqrt(sum);
        if(dist < min_dist){
            min_dist = dist;
            closest_centroid = i;
        }
    }
    //assegno id-cluster al thread corrente
    if( d_clusterval[idx] != closest_centroid) {
        d_clusterval[idx] = closest_centroid;
        *d_done = false;
    }
}

__global__ void clusterPointsSum(float* d_xval, float* d_yval, int* d_clusterVal, float* d_clusterSumX, float* d_clusterSumY, int* d_clusterSize){
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
    float centroids_shifts_sum = 0;
    float max_tollerance = 0.0001;
    //alloco memoria Host

    float *h_xval = (float *) malloc(N * sizeof(float));
    float *h_yval = (float *) malloc(N * sizeof(float));
    int *h_clusterVal = (int *) malloc(N * sizeof(int));
    //float*h_centroids =
    float *h_centroidX = (float *) malloc(K * sizeof(float));
    float *h_centroidY = (float *) malloc(K * sizeof(float));
    float *h_prevCentroidX = (float *) malloc(K * sizeof(float));
    float *h_prevCentroidY = (float *) malloc(K * sizeof(float));

    int *h_clusterSize = (int *) malloc(K * sizeof(int));
    float *h_clusterSumX = (float *) malloc(K * sizeof(float));
    float *h_clusterSumY = (float *) malloc(K * sizeof(float));
    bool *h_done = (bool *)malloc(sizeof(bool));

    //alloco memoria Device
    float *d_xval;
    float *d_yval;
    //float*h_centroids =
    int *d_clusterVal;
    float *d_centroidX;
    float *d_centroidY;
    int *d_clusterSize;
    float *d_clusterSumX;
    float *d_clusterSumY;
    bool *d_done;
    cudaMalloc(&d_xval, total_points * sizeof(float));
    cudaMalloc(&d_yval, total_points * sizeof(float));
    cudaMalloc(&d_centroidX, K * sizeof(float));
    cudaMalloc(&d_centroidY, K * sizeof(float));
    cudaMalloc(&d_clusterVal, N * sizeof(int));
    cudaMalloc(&d_clusterSize, K * sizeof(int));
    cudaMalloc(&d_clusterSumX, K * sizeof(float));
    cudaMalloc(&d_clusterSumY, K * sizeof(float));
    cudaMalloc(&d_done, sizeof(bool));

    readPoints(h_xval, h_yval, h_clusterVal);       //leggo i punti dataset

    dimensions = 2;   //dimensione punto

    //std::cout<<"dimensione all_points: "<<total_points<<" dimensione clusters: "<<all_points.getDimClusters()<<" punto 1: "<<all_points.getYval(1)<< std::endl;
    // Inizializzo Clusters
    std::vector<int> used_pointIds;  //mem punti gia usati per init cluster

    bool exit = false;
    float x, y;

    used_pointIds = indexGenerator(K, N);
    for (int i = 0; i < K; i++) {

        h_centroidX[i] = h_xval[used_pointIds[i]];
        h_centroidY[i] = h_yval[used_pointIds[i]];
        //x = all_points.getXval(used_pointIds[i]);
        //y = all_points.getYval(used_pointIds[i]);
        h_clusterVal[used_pointIds[i]] = i;
        h_clusterSize[i] = 0;
    }
    std::cout << "Clusters Inizializzati = " << std::endl
              << std::endl;

    std::cout << "Eseguo Clustering K-Means..." << std::endl;

    int iter = 1;
    *h_done = true;
    //copia host -> device

    cudaMemcpy(d_xval, h_xval, N * sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(d_yval, h_yval, N * sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(d_centroidX, h_centroidX, K * sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(d_centroidY, h_centroidY, K * sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(d_clusterSize, h_clusterSize, K * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_clusterSumX, h_centroidX, K * sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(d_clusterSumY, h_centroidY, K * sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(d_done, h_done,sizeof(bool), cudaMemcpyHostToDevice);

    while (true) {

        std::cout << "Iter - " << iter << "/" << iters << std::endl;
        bool done = true;
        clusterAssignment<<<(N + TPB - 1) / TPB, TPB>>>(d_xval, d_yval, d_clusterVal, d_centroidX, d_centroidY, K, d_done);

        //setto a 0 la sommma dei punti appartenenti ai cluster, per iniziare aggiornamento centroide
        cudaMemset(d_clusterSumX, 0.0, K * sizeof(float));
        cudaMemset(d_clusterSumY, 0.0, K * sizeof(float));

        clusterPointsSum<<<(N + TPB - 1) / TPB, TPB>>>(d_xval, d_yval, d_clusterVal, d_clusterSumX, d_clusterSumY,
                                                       d_clusterSize);

        //copia variabili aggiornate da device -> host

        cudaMemcpy(h_clusterSumX, d_clusterSumX, K * sizeof(float), cudaMemcpyDeviceToHost);
        cudaMemcpy(h_clusterSumY, d_clusterSumY, K * sizeof(float), cudaMemcpyDeviceToHost);
        cudaMemcpy(h_clusterSize, d_clusterSize, K * sizeof(int), cudaMemcpyDeviceToHost);
        cudaMemcpy(h_done, d_done, sizeof(bool), cudaMemcpyDeviceToHost);

        cudaMemset(d_clusterSize, 0, K * sizeof(int));

        //Ricalcolo nuovi centroidi per ogni cluster
        for (int i = 0; i < K; i++) {
            h_centroidX[i] = h_clusterSumX[i] / h_clusterSize[i];
            h_centroidY[i] = h_clusterSumY[i] / h_clusterSize[i];

        }

        cudaMemcpy(d_centroidX, h_centroidX, K * sizeof(float), cudaMemcpyHostToDevice);
        cudaMemcpy(d_centroidY, h_centroidY, K * sizeof(float), cudaMemcpyHostToDevice);




        float deltaX, deltaY = 0;
        for(int i= 0; i < K; i++)
        {
                deltaX = h_centroidX[i] - h_prevCentroidX[i];
                deltaY = h_centroidY[i] - h_prevCentroidY[i];
                double percentageShiftX = std::abs((deltaX / h_centroidX[i]) * 100.0);
                double percentageShiftY = std::abs((deltaY / h_centroidY[i]) * 100.0);
                std::cout<< "perX: "<< percentageShiftX << " percY: "<< percentageShiftY<< std::endl;
                if(percentageShiftX > max_tollerance || percentageShiftY > max_tollerance)
                {
                    done = false;

                }
            }

        for (int i = 0; i < K; i++) {
            h_prevCentroidX[i] = h_centroidX[i];
            h_prevCentroidY[i] = h_centroidY[i];

        }
        if (done || iter >= iters) {
            std::cout << "Clustering completed in iteration : " << iter << std::endl
                      << std::endl;
            break;
        }
        iter++;

    }
    cudaMemcpy(h_clusterVal, d_clusterVal, N * sizeof(int), cudaMemcpyDeviceToHost);


    ///scrivo i risultati clustering su file
    //scrittura centroidi
    std::ofstream outfile;
    //std::cout << std::filesystem::current_path().string() << std::endl;
    //std::cout << output_dir + "/" + "clusters.txt" << std::endl;
    outfile.open(output_dir + "/" + "clusters.txt");
    if (outfile.is_open()) {
        for (int i = 0; i < K; i++) {
            //std::cout <<  i << " cluster contiene: "<< clusters[i].getSize() <<std::endl;
            std::cout << "Cluster " << i << " centroid : ";

            std::cout << h_centroidX[i] << " "<< h_centroidY[i] << std::endl;    // Output console

            outfile << h_centroidX[i] << " " << h_centroidY[i]; // Output file
            outfile << std::endl;
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

    cudaFree(d_xval);
    cudaFree(d_yval);
    cudaFree(d_clusterVal);
    cudaFree(d_centroidX);
    cudaFree(d_centroidY);
    cudaFree(d_clusterSize);

    free(h_xval);
    free(h_yval);
    free(h_clusterVal);
    free(h_centroidX);
    free(h_centroidY);
    free(h_clusterSize);

}



float averageParallelExecutions(int K, int iters, std::string output_dir, std::string input_dir)
{
    float mediaS, mediaP;
    float sum;

    for(int i  = 0; i < 2; i++ ) {
        auto start = std::chrono::high_resolution_clock::now();
        ParallelKMeans kmeans(K, iters, output_dir, input_dir);
        //kmeans.run_parallel2(all_points);
        kmeans.run();
        //float end = omp_get_wtime( );
        auto end = std::chrono::high_resolution_clock::now();

        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
        std::cout << "Tempo di esecuzione parallela: " << duration.count() << " millisecondi" << std::endl;
        sum += static_cast<float>(duration.count());
        std::cout << "<<------------------------------>>" << std::endl;

    }
    mediaP = static_cast<float>(sum)/2;
    return mediaP;
}

float averageSeqExecutions(int K, int iters, std::string output_dir, std::string input_dir)
{
    float mediaS;
    float sum;

    for(int i  = 0; i < 2; i++ ) {
        auto start = std::chrono::high_resolution_clock::now();
        SequentialKMeans kmeans(K, iters, output_dir, input_dir);
        kmeans.run();
        //float end = omp_get_wtime( );
        auto end = std::chrono::high_resolution_clock::now();

        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
        std::cout << "Tempo di esecuzione Sequenziale: " << duration.count() << " millisecondi" << std::endl;
        sum += static_cast<float>(duration.count());
        std::cout << "<<------------------------------>>" << std::endl;

    }
    mediaS = static_cast<float>(sum)/2;
    return mediaS;
}

int main() {
    std::string output_dir = "../cmake-build-debug/cluster_details";   //dir output
    int K = 3;                               //numero cluster
    std::string input_dir= "input1.txt";

    // Avvio il clustering
    int iters = 100;


    auto mediaS = averageSeqExecutions(K, iters, output_dir, input_dir);
    auto mediaP = averageParallelExecutions(K, iters, output_dir, input_dir);
    std::cout << "Media esecuzione Sequenziale : " << mediaS << std::endl;
    std::cout << "Media esecuzione Parallela : " << mediaP << std::endl;

    float speedup = static_cast<float>(mediaS) / static_cast<float>(mediaP);
    std::cout << "Speedup: " << speedup << std::endl;

    return 0;
}
