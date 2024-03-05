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
    std::vector<double> centroid;
    std::vector<int>pointsIds;

public:
    Cluster(int clusterId, int idxP, double centX, double centY);

    void addPoint(int idx);
    bool removePoint(int idx);


    int getId();

    double getIdByPos(int pos);
    int getSize();
    double getCentroidByPos(int pos);

    void setCentroid(double x, double y) ;
};



#endif //KMEANSPARALLEL_CLUSTER_H
