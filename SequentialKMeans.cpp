//
// Created by Giacomo Magistrato on 28/02/24.
//

#include "SequentialKMeans.h"

#include "Cluster.h"
#include "Point.h"
#include <omp.h>
#include <cmath>
#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <chrono>
#include <filesystem>

#include <algorithm>
static std::vector<int> used_pointIds;
static bool first = true;


std::vector<int>index_Generator(int K, int total_points){
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

int SequentialKMeans::getNearestClusterId(int idx) {

    float min_dist = INFINITY;
    int NearestClusterId;
    float dist;
    for (int i = 0; i < K; i++) {
        float sum = 0.0;
        sum += pow(clusters[i].getCentroidByPos(0) - all_points.getXval(idx), 2.0);
        sum += pow(clusters[i].getCentroidByPos(1) - all_points.getYval(idx), 2.0);
        dist = sqrt(sum);
        if (dist < min_dist) {
            min_dist = dist;
            NearestClusterId = i;
        }

    }
    return NearestClusterId;
}

SequentialKMeans::SequentialKMeans(int K, int iterations, std::string output_dir, std::string input_dir, std::vector<int>used_pointIds):all_points(input_dir) {
    this->K = K;
    this->iters = iterations;
    this->output_dir = output_dir;
    this->used_pointIds = used_pointIds;
    //this->all_points = Point(input_dir);
}
///////////////////////////////// SEQUENZIALE \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ TODO: termina prevCentroids

void SequentialKMeans::run() {
    total_points = all_points.getDimensions();  //num totale di punti

    dimensions = 2;   //dimensione punto
    float max_tollerance = 0.0001;
    std::vector<float>prevCentroidsX;
    std::vector<float>prevCentroidsY;
    //std::cout<<"dimensione all_points: "<<total_points<<" dimensione clusters: "<<all_points.getDimClusters()<<" punto 1: "<<all_points.getYval(1)<< std::endl;
    // Inizializzo Clusters
    //std::vector<int> used_pointIds;  //mem punti gia usati per init cluster

    float x, y;

    //used_pointIds = index_Generator(K, total_points);
    for (int i = 0; i < K; i++) {

        x = all_points.getXval(used_pointIds[i]);
        y = all_points.getYval(used_pointIds[i]);
        all_points.setCluster(used_pointIds[i], i);
        //Cluster cluster(i,used_pointIds[i], x, y);  //creo un cluster avente centroide il punto attuale
        Cluster cluster(i, x, y);
        //cluster.addPoint();
        std::cout<<"init: x:"<< cluster.getCentroidByPos(0)<<" y: "<<cluster.getCentroidByPos(1)<<std::endl;
        prevCentroidsX.push_back(x);
        prevCentroidsY.push_back((y));
        clusters.push_back(cluster);
    }

    std::cout << "Clusters Inizializzati = " << clusters.size() << std::endl
              << std::endl;

    std::cout << "Eseguo Clustering K-Means..." << std::endl;

    int iter = 1;
    std::vector<float>clusterSumX;
    std::vector<float>clusterSumY;
    for(int i = 0; i<K; i++){
        clusterSumX.push_back(0);
        clusterSumY.push_back(0);
    }
    while (true) {
        auto start = std::chrono::high_resolution_clock::now();

        std::cout << "Iter - " << iter << "/" << iters << std::endl;
        bool done = true;

        // Aggiungo punti al cluster più vicino
        int j = 0;
//#pragma omp parallel for reduction(&&: done) num_threads(16)
/*
        for (int i = 0; i < total_points; i++) {
            int currentClusterId = all_points.getCluster(i);   //ottengo clusterID punto corrente
            int nearestClusterId = getNearestClusterId(i);  //cluster più vicino

            //se il cluster del punto != dal cluster più vicino, setto sul punto il cluster giusto
            if (currentClusterId != nearestClusterId) {
                j++;
                all_points.setCluster(i, nearestClusterId);
                clusters[nearestClusterId].addPoint(i);
                if(currentClusterId != -1) {  //rimuovo l'indice punto dal cluster precedente
                    clusters[currentClusterId].removePoint(i);
                }
                //done = false;

            }
        }
        */
        for (int i = total_points-1; i > 0; i--) {
            int currentClusterId = all_points.getCluster(i);   //ottengo clusterID punto corrente
            int nearestClusterId = getNearestClusterId(i);  //cluster più vicino

            //se il cluster del punto != dal cluster più vicino, setto sul punto il cluster giusto
            if (currentClusterId != nearestClusterId) {
                j++;
                all_points.setCluster(i, nearestClusterId);
                //clusters[nearestClusterId].addPoint();
                if(currentClusterId != -1) {  //rimuovo l'indice punto dal cluster precedente
                    //clusters[currentClusterId].removePoint();
                }
                done = false;

            }
        }

/*

        // Ricalcolo nuovi Centroidi
        int somma = 0;
        for (int i = 0; i < K; i++) {

            int ClusterSize = clusters[i].getClusterSize();

            float sumX = 0.0, sumY = 0.0;

            if (ClusterSize > 0) {
                for (int p = 0; p < ClusterSize; p++) {
                    sumX += all_points.getXval(clusters[i].getIdByPos(p));
                    sumY += all_points.getYval(clusters[i].getIdByPos(p));
                }
                clusters[i].setCentroid(sumX / ClusterSize, sumY / ClusterSize);
                std::cout<<clusters[i].getCentroidByPos(0)<<"  "<<clusters[i].getCentroidByPos(1)<<std::endl;
            }
            //std::cout<<"Cluster "<<i <<": "  <<clusters[i].getSize() << std::endl;
            somma += clusters[i].getSize();
        }

*/
        // Ricalcolo nuovi Centroidi

        //float *clusterSumX = (float *) malloc(K * sizeof(float));
        //float *clusterSumY = (float *) malloc(K * sizeof(float));
        for(int i = 0; i< K; i++)
        {
            clusterSumX[i] = 0;
            clusterSumY[i] = 0;
        }

        int clusterId;
        for (int i = total_points-1; i > 0; i--)
        {
            x = all_points.getXval(i);
            y = all_points.getYval(i);
            clusterId = all_points.getCluster(i);
            clusterSumX[clusterId] += x;
            clusterSumY[clusterId] += y;
            clusters[clusterId].addPoint();
        }

        for(int i = 0; i <K; i++)
        {
            //std::cout<< "cluster: "<<clusters[i].getId()<<" centroidX: "<<clusters[i].getCentroidByPos(0)<<" centroidY: "<<clusters[i].getCentroidByPos(1)<<" "<<clusters[i].getClusterSize()<<std::endl;
            clusters[i].setCentroid(clusterSumX[i]/clusters[i].getClusterSize(), clusterSumY[i]/clusters[i].getClusterSize());
            //std::cout<<"nuovi: x:"<< clusters[i].getCentroidByPos(0)<<" y: "<<clusters[i].getCentroidByPos(1)<<std::endl;
        }
        //std::cout<<"ci sono "<< std::endl;

        /*
        int somma = 0;
        for (int i = 0; i < K; i++) {

            int ClusterSize = clusters[i].getClusterSize();

            float sumX = 0.0, sumY = 0.0;

            if (ClusterSize > 0) {
                for (int p = 0; p < ClusterSize; p++) {
                    sumX += all_points.getXval(clusters[i].getIdByPos(p));
                    sumY += all_points.getYval(clusters[i].getIdByPos(p));
                }
                clusters[i].setCentroid(sumX / ClusterSize, sumY / ClusterSize);
                std::cout<<clusters[i].getCentroidByPos(0)<<"  "<<clusters[i].getCentroidByPos(1)<<std::endl;
            }
            //std::cout<<"Cluster "<<i <<": "  <<clusters[i].getSize() << std::endl;
            somma += clusters[i].getSize();
        }
        */
        //std::cout<<"elementi Clusters "<< somma << std::endl;
        float deltaX, deltaY = 0;
        for(int i= 0; i < K; i++)
        {
            deltaX = clusters[i].getCentroidByPos(0) - prevCentroidsX[i];
            deltaY = clusters[i].getCentroidByPos(1) - prevCentroidsY[i];
            double percentageShiftX = std::abs((deltaX / clusters[i].getCentroidByPos(0)) * 100.0);
            double percentageShiftY = std::abs((deltaY / clusters[i].getCentroidByPos(1))*100);
            //std::cout<< "perX: "<< percentageShiftX << " percY: "<< percentageShiftY<< std::endl;
            if(percentageShiftX > max_tollerance || percentageShiftY > max_tollerance)
            {
                done = false;

            }
        }
        for (int i = 0; i < K; i++) {
            prevCentroidsX[i] = clusters[i].getCentroidByPos(0);
            prevCentroidsY[i] = clusters[i].getCentroidByPos(1);
            clusters[i].setSize(0);
        }

        //interrompo se clustering completo o num iterazioni raggiunto
        if (done || iter >= iters) {
            std::cout << "Clustering completed in iteration : " << iter << std::endl
                      << std::endl;
            break;
        }
        iter++;
        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
        std::cout << "Tempo di esecuzione run():  " << duration.count()  << "  millisecondi" << std::endl;
    }


    // Scrivo i risultati del clustering su file (centroidi finali)
    std::ofstream outfile;
    std::cout << std::filesystem::current_path().string() << std::endl;
    std::cout<<output_dir+ "/" +"clustersS.txt"<<std::endl;
    outfile.open(output_dir + "/" +"clustersS.txt");
    if (outfile.is_open()) {
        for (int i = 0; i < K; i++) {
            //std::cout <<  i << " cluster contiene: "<< clusters[i].getSize() <<std::endl;
            std::cout << "Cluster " << clusters[i].getId() << " centroid : ";

            for (int j = 0; j < dimensions; j++) {
                std::cout << clusters[i].getCentroidByPos(j) << " ";    // Output console
                outfile << clusters[i].getCentroidByPos(j) << " "; // Output file
            }
            std::cout << std::endl;
            outfile << std::endl;
        }
        outfile.close();
    } else {
        std::cout << "Error: Unable to write to clusters.txt";
    }


    outfile.open(output_dir + "/" + "clusteringS.txt");
    for (int i = 0; i < total_points; i++) {
        // indice cluster è ID-1
        if (outfile.is_open()) {

            outfile << all_points.getXval(i) << " " << all_points.getYval(i) << " "
                    << all_points.getCluster(i)  ;// Output to file
        }

        outfile << std::endl;

    }
    outfile.close();
}