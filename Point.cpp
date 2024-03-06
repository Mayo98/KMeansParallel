//
// Created by Giacomo Magistrato on 28/02/24.
//

#include "Point.h"

#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <filesystem>
using namespace std;

void Point::lineToVec(std::string& dir) {
    std::vector<float> values;
    std::string tmp = "";
    string line;
    std::cout << std::filesystem::current_path().string() <<"Ciao" << std::endl;
    ifstream infile("../cmake-build-debug/input1.txt");
    if (!infile.is_open()) {
        cout << "Error: Failed to open file." << endl;
        return;
    }
    while (getline(infile, line)) {
        for (int i = 0; i < static_cast<int>(line.length()); i++) {
            if ((48 <= static_cast<int>(line[i]) && static_cast<int>(line[i]) <= 57) || line[i] == '.' ||
                line[i] == '+' || line[i] == '-' || line[i] == 'e') {
                tmp += line[i];

            } else if (!tmp.empty()) {

                values.push_back(std::stod(tmp));
                xval.push_back(std::stod(tmp));
                tmp = "";
            }
        }
        if (!tmp.empty()) {

            values.push_back(std::stod(tmp));
            yval.push_back(std::stod(tmp));
            tmp = "";
        }
        this->clusters.push_back(-1); //inizializzo anche il vettore clusters a -1
    }
    std::cout<<"dim clusters: "<<clusters.size()<<std::endl;

    cout << "\nDataset fetched!" << endl
         << endl;
    std::cout << "Punti totali letti : " << xval.size() << std::endl;

    /* // Esco se ho più cluster che punti
     if ((int) xval.size() < K) {
         cout << "Errore: Numero di clusters maggiore del numero di punti." << endl;
         return 1;
     }*/

}
/*
 for (int i = 0; i < static_cast<int>(line.length()); i++)
{
if ((48 <= static_cast<int>(line[i]) && static_cast<int>(line[i]) <= 57) || line[i] == '.' || line[i] == '+' || line[i] == '-' || line[i] == 'e')
{
    tmp += line[i];
}
else if (!tmp.empty())
{
    float value = std::stod(tmp);
    xval.push_back(value);
    tmp = ""; // Reimposta tmp dopo aver salvato il valore in xval

    // Aggiungi qui la logica per salvare il valore in yval, se necessario
}
}

if (!tmp.empty())
{
    float value = std::stod(tmp);
    yval.push_back(value);
    tmp = ""; // Reimposta tmp dopo aver salvato il valore in yval
}

*/



Point::Point(std::string dir)
{
    this->lineToVec(dir);
    dimensions = xval.size();

    //this->initCluster(dimensions);
}

int Point::getNumPoints()
{
    return xval.size();
}
void Point::initCluster(int dimensions){
    clusters[dimensions] = {0};
    for(int i = 0; i < clusters.size();i++)
    {
        std::cout<<"clus: "<< clusters[i]<<std::endl;
    }
}
int Point::getDimensions()
{
    return dimensions;
}

int Point::getCluster(int idx)
{
    return clusters[idx];
}

int Point::getDimClusters(){return clusters.size();}
void Point::setCluster(int idx, int val)
{
    clusters[idx] = val;        //setto il nuovo id del cluster per il punto
}

float Point::getXval(int pos)
{
    return xval[pos];
}
float Point::getYval(int pos)
{
    return yval[pos];
}

