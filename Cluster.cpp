//
// Created by Giacomo Magistrato on 28/02/24.
//


#include "Cluster.h"
#include "Point.h"
#include <iostream>
#include <fstream>
Cluster::Cluster(int clusterId, float centX, float centY)
{
    this->clusterId = clusterId;
    this->centroid.push_back(centX);
    this->centroid.push_back(centY);
    this->clusterSize++;
}
/*
void Cluster::addPoint(int idxP)
{
    pointsIds.push_back(idxP);
    //points[idxP] = std::make_pair(x, y);
    //pointsX.push_back(x);
    //pointsY.push_back(y);
}
*/
void Cluster::addPoint()
{
    this->clusterSize++;
}
/*
void Cluster::removePoint(int idx)
{

    int size = pointsIds.size();
    std::cout<<"punti init : "<< size<<std::endl;
    for (int i = 0; i < size; i++)
    {
        if (pointsIds[i] == idx)
        {
            pointsIds.erase(pointsIds.begin() + i);
            std::cout<<"punti final : "<< pointsIds.size()<<std::endl;
            return;
        }
    }
    return;
}
*/
void Cluster::removePoint()
{
    this->clusterSize--;
}
//void Cluster::removeAllPoints() { pointsX.clear(); pointsY.clear(); }

int Cluster::getId() { return clusterId; }

float Cluster::getIdByPos(int pos) { return pointsIds[pos]; }

/////Metodi con mappa
/*
   void Cluster::removePointbyIdx(int idx) {
    points.erase(idx);
    //pointsX.erase(pointsX.begin() + idx);
    //pointsY.erase(pointsY.begin() + idx);
}
  float Cluster::getSumX() {
    float sum = 0.0;
    for(auto&coppia: points){
        sum+= coppia.second.first;
    }
    return sum;
}
float Cluster::getSumY() {
    float sum = 0.0;
    for(auto&coppia: points){
        sum+= coppia.second.second;
    }
    return sum;
}

void Cluster::saveResult(std::ofstream& outfile, std::string path, std::string path2) {

    std::cout << "Cluster " << this->getId() << " centroid : ";
    outfile.open(path, std::ios_base::app);
    if (outfile.is_open()) {
        for (int i = 0; i < centroid.size(); i++) {
            std::cout << getCentroidByPos(i) << " ";    // Output console
            outfile << getCentroidByPos(i) << " "; // Output file
        }
        std::cout << std::endl;
        outfile << std::endl;

    }else {
        std::cout << "Error: Unable to write to clusters.txt";
    }
    outfile.open(path2, std::ios_base::app);
    if(outfile.is_open()) {
        for (auto &elem: points) {
            if (outfile.is_open()) {
                outfile << elem.second.first << " " << elem.second.second << " " << this->getId()
                        << std::endl;
            }
        }
    }else{std::cout << "Error: Unable to write to clusters.txt";}
    outfile.close();


}
 */ /////////

int Cluster::getSize() { return pointsIds.size(); }
int Cluster::getClusterSize() {return clusterSize;}
float Cluster::getCentroidByPos(int pos) { return centroid[pos]; }
void Cluster::setSize(int size ) {this->clusterSize = size;}
void Cluster::setCentroid(float x, float y) { this->centroid[0] = x; this->centroid[1] = y; }
