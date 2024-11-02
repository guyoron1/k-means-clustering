#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define DEFAULT_ITER 200
#define EPSILON 0.001
#define ERROR_MESSAGE "An Error Has Occurred"

typedef struct{
    double* centroid;
    double* prevCentroid;
    double** dataPointList;
    int dataPointsSize;
}CLUSTER;

void clearClusterDataPointList(CLUSTER *cluster){
    cluster -> dataPointsSize = 0;
}

void updateClusterCentroid(CLUSTER *cluster, int d){
    int i = 0;
    int assignedPoints = cluster-> dataPointsSize;

    if (assignedPoints == 0){
        printf(ERROR_MESSAGE);
        return;
    }
  
    for (i = 0; i < d; i++){
        double sum=0.0;
        int j;
        for (j = 0; j < assignedPoints; j++){
            sum += cluster -> dataPointList[j][i];
        }
        cluster->prevCentroid[i] = cluster->centroid[i];
        cluster -> centroid[i] = sum/(assignedPoints);
    }
}

void addDataPointToCluster(CLUSTER *cluster, double* dataPoint){
    cluster -> dataPointList[cluster->dataPointsSize] = dataPoint;
    cluster -> dataPointsSize++;
}

double calcDistance(double* vector1, double* vector2, int d){
    double sum = 0.0;
    int i=0;
    for (i = 0; i < d; i++){
        sum += pow((vector1[i] - vector2[i]), 2);
    }
    return sqrt(sum);
}

int hasClustersConverged(int K, int d, CLUSTER* clusters){
    int i;
    for (i = 0; i < K; i++){
        if (calcDistance(clusters[i].centroid, clusters[i].prevCentroid, d) >= EPSILON){
            return 0;
        }
    }
    return 1;
}

void assignClusterToDataPoint(CLUSTER* clusters, double** dataPoints, int N, int K, int d){
    int pointIndex;
    for (pointIndex = 0;  pointIndex < N; pointIndex++){
        int clusterIndex;
        int minClusterIdx = 0;
        double minDistance = calcDistance(dataPoints[pointIndex], clusters[0].centroid, d);
        for (clusterIndex = 1; clusterIndex < K; clusterIndex++){
            double distance = calcDistance(dataPoints[pointIndex], clusters[clusterIndex].centroid, d);
            if (distance < minDistance){
                minDistance = distance;
                minClusterIdx = clusterIndex;
            }
        
        }
        addDataPointToCluster(&clusters[minClusterIdx], dataPoints[pointIndex]);
    }
}

void freeDataPoints(double** dataPoints, int N){
    int i;
    for (i = 0; i < N; i++){
        free(dataPoints[i]);
    }
    free(dataPoints);
}

void freeClusters(CLUSTER* clusters, int K) {
    int i;
    for (i = 0; i < K; i++) {
        free(clusters[i].centroid);
        free(clusters[i].prevCentroid);
        free(clusters[i].dataPointList);
    }
    free(clusters);
}

CLUSTER* initializeClusters(int K, int N, int d, double** dataPoints){
    CLUSTER* clusters = (CLUSTER*)malloc(K * sizeof(CLUSTER));
    int i = 0;
    if (!clusters) {
        printf(ERROR_MESSAGE);
        return NULL;
    }

    for(i=0; i<K; i++){
        int j = 0;
        clusters[i].centroid = (double*)malloc(d * sizeof(double));
        clusters[i].prevCentroid =(double*) malloc(d * sizeof(double));
        clusters[i].dataPointList = (double**)malloc(N *sizeof(double*));
        clusters[i].dataPointsSize = 0;
        if (!clusters[i].centroid || !clusters[i].prevCentroid || !clusters[i].dataPointList) {
            printf(ERROR_MESSAGE);
            freeClusters(clusters, i+1);
            return NULL;
        }

        for (j = 0; j < d; j++) {
            clusters[i].centroid[j] = dataPoints[i][j];
            clusters[i].prevCentroid[j] = dataPoints[i][j];
        }
    }
    return clusters;
}

void printClusters(int K, int d, CLUSTER* clusters){
    int clusterIdx;
    for (clusterIdx = 0; clusterIdx < K; clusterIdx++){
        int i;
        for (i = 0; i < d; i++){
            printf("%.4f", clusters[clusterIdx].centroid[i]);
            if (i<d-1){
                printf(",");
            } 
        }
        if (clusterIdx < K-1){
            printf("\n");
        }
    }
}

int solveKMeanAlgo(int K, int N, int d, int iter, double** dataPoints){
    int counter = 0;
    CLUSTER* clusters = initializeClusters(K, N, d, dataPoints);
    if (!clusters){
        return 1;
    }

    for (counter = 0; counter < iter; counter++){
        int i=0;
        assignClusterToDataPoint(clusters, dataPoints, N, K, d);
    
        for (i = 0; i < K; i++){
            updateClusterCentroid(&clusters[i], d);
            clearClusterDataPointList(&clusters[i]);
        }

        if (hasClustersConverged(K, d, clusters)){
            break;
        }
    }
    printClusters(K, d, clusters);
    freeClusters(clusters, K);   
    return 0;
}

double** allocateDataPointsMemory(int d, int N){
    double** dataPoints;
    int i = 0;

    dataPoints = (double**)calloc(N, sizeof(double*));
    if (!dataPoints){
        printf(ERROR_MESSAGE);
        return NULL;
    }
    
    for (i = 0; i < N; i++){
        dataPoints[i] = (double*)calloc(d, sizeof(double));
        if (!dataPoints[i]){
            freeDataPoints(dataPoints, i);
            printf(ERROR_MESSAGE);
            return NULL;
        }
    }
    return dataPoints;  
}

double** retrieveDataPoints(int N, int d){
    int i = 0;
    int j = 0;
    double** dataPoints = allocateDataPointsMemory(d, N);

    if (!dataPoints){
        return NULL;
    }

    for (i = 0; i < N; i++){
        for (j = 0; j < d; j++){
            if (scanf("%lf,", &(dataPoints[i][j]))!=1){
                scanf("%lf\n", &(dataPoints[i][j]));
            }    
        }
    }

    return dataPoints;
}

int isInteger(const char *str) {
    if (!*str)
        return 0;

    while (*str) {
        char c = *str;
        if (c < '0' || c > '9')
        {
           return 0;
        }
        str++;
    }
    return 1;
}

int* retrieveAndValidateInputs(int argc, char** argv) {
    int N = 0;
    int K = 0;
    int d = 0;
    int iter = DEFAULT_ITER;
    int* inputs = NULL;

    if (argc < 4) {
        printf(ERROR_MESSAGE);
        return NULL;
    }
    
    if (!isInteger(argv[2]) || atoi(argv[2]) <= 1) {
        printf("Invalid number of points!");
        return NULL;
    }
    N = atoi(argv[2]);

    if (!isInteger(argv[1]) || atoi(argv[1]) <= 1 || atoi(argv[1]) >= N) {
        printf("Invalid number of clusters!");
        return NULL;
    }
    K = atoi(argv[1]);

    if (!isInteger(argv[3]) || atoi(argv[3]) < 1) {
        printf("Invalid dimension of point!");
        return NULL;
    }
    d = atoi(argv[3]);

    if (argc > 4) {
        iter = atoi(argv[4]);
        if (!isInteger(argv[4]) || iter <= 1 || iter >= 1000) {
            printf("Invalid maximum iteration!");
            return NULL;
        }
    }

    inputs = (int*)malloc(4 * sizeof(int));
    if (!inputs) {
        printf(ERROR_MESSAGE);
        return NULL;
    }

    inputs[0] = K;
    inputs[1] = N;
    inputs[2] = d;
    inputs[3] = iter;

    return inputs;
}

int main(int argc, char** argv){
    double** dataPoints = NULL;
    int res = 0;
    int* inputs = retrieveAndValidateInputs(argc, argv);

    if (!inputs){
        return 1;
    }
    dataPoints = retrieveDataPoints(inputs[1], inputs[2]);
    if (!dataPoints){
        free(inputs);
        return 1;
    }

    res = solveKMeanAlgo(inputs[0], inputs[1], inputs[2], inputs[3], dataPoints);

    freeDataPoints(dataPoints, inputs[1]);
    free(inputs);
    return res;
}