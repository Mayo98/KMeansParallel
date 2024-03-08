//
// Created by Giacomo Magistrato on 28/02/24.
//

#ifndef KMEANSPARALLEL_CLUSTER_H
#define KMEANSPARALLEL_CLUSTER_H
#include "Point.h"
#include <vector>
#include <string>
#include <unordered_map>

class Cluster {
private:
    int clusterId;
    std::vector<float> centroid;
    std::vector<int>pointsIds;
    int clusterSize = 0;

public:
    //Cluster(int clusterId, int idxP, float centX, float centY);
    Cluster(int clusterId, float centX, float centY);
    //void addPoint(int idx);
    void addPoint();
    //void removePoint(int idx);
    void removePoint();
    void setSize(int size);
    int getClusterSize();
    int getId();

    float getIdByPos(int pos);
    int getSize();
    float getCentroidByPos(int pos);

    void setCentroid(float x, float y) ;
};



#endif //KMEANSPARALLEL_CLUSTER_H
