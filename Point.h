//
// Created by Giacomo Magistrato on 28/02/24.
//

#ifndef KMEANSPARALLEL_POINT_H
#define KMEANSPARALLEL_POINT_H

//#ifndef KMEANS_PARALLEL_POINT_H
//#define KMEANS_PARALLEL_POINT_H
#ifndef POINT_H
#define POINT_H

#include <vector>
#include <string>

class Point
{
private:
    int dimensions ;
    //std::vector<float> values;
    void lineToVec(std::string& line);

    std::vector<float>xval;
    std::vector<float>yval;
    std::vector<int>clusters;

public:
    Point(std::string dir);
    int getDimClusters();
    int getDimensions();
    int getCluster(int idx);
    void setCluster(int idx, int val);
    int getNumPoints();
    float getXval(int pos);
    void initCluster(int dimensions);
    float getYval(int pos);
};

#endif  // POINT_H



#endif //KMEANSPARALLEL_POINT_H
