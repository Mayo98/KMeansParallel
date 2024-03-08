//
// Created by Giacomo Magistrato on 28/02/24.
//

#ifndef KMEANSPARALLEL_SEQUENTIALKMEANS_H
#define KMEANSPARALLEL_SEQUENTIALKMEANS_H



#include "Cluster.h"
#include "Point.h"
#include <omp.h>
#include <cmath>
#include <iostream>
#include <vector>
#include <string>

class SequentialKMeans {
private:
    int K, iters, dimensions{}, total_points{};
    std::vector<Cluster> clusters;
    std::string input_dir;
    std::string output_dir;
    Point all_points;


    int getNearestClusterId(int pointIdx);
    int getNearestClusterId_Parallel(Point point);

public:
    SequentialKMeans(int K, int iterations, std::string output_dir, std::string input_dir, std::vector<int>used_pointIds);
    std::vector<int>used_pointIds;
    void run();
    void run_parallel(std::vector<Point> &all_points);
    void run_parallel2(std::vector<Point> &all_points);
};



#endif //KMEANSPARALLEL_SEQUENTIALKMEANS_H
