//
// Created by Giacomo Magistrato on 28/02/24.
//


#include "Cluster.h"
#include "Point.h"
#include <iostream>
#include <fstream>
Cluster::Cluster(int clusterId, int idxP, double centX, double centY)
{
    this->clusterId = clusterId;
    this->centroid.push_back(centX);
    this->centroid.push_back(centY);
    this->addPoint(idxP);
}

void Cluster::addPoint(int idxP)
{
    pointsIds.push_back(idxP);
    //points[idxP] = std::make_pair(x, y);
    //pointsX.push_back(x);
    //pointsY.push_back(y);
}

bool Cluster::removePoint(int idx)
{

    int size = pointsIds.size();

    for (int i = 0; i < size; i++)
    {
        if (pointsIds[i] == idx)
        {
            pointsIds.erase(pointsIds.begin() + i);
            return true;
        }
    }
    return false;
}

//void Cluster::removeAllPoints() { pointsX.clear(); pointsY.clear(); }

int Cluster::getId() { return clusterId; }

double Cluster::getIdByPos(int pos) { return pointsIds[pos]; }

/////Metodi con mappa
/*
   void Cluster::removePointbyIdx(int idx) {
    points.erase(idx);
    //pointsX.erase(pointsX.begin() + idx);
    //pointsY.erase(pointsY.begin() + idx);
}
  double Cluster::getSumX() {
    double sum = 0.0;
    for(auto&coppia: points){
        sum+= coppia.second.first;
    }
    return sum;
}
double Cluster::getSumY() {
    double sum = 0.0;
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

double Cluster::getCentroidByPos(int pos) { return centroid[pos]; }

void Cluster::setCentroid(double x, double y) { this->centroid[0] = x; this->centroid[1] = y; }
