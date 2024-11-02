import sys
import numpy as np

DEAFULT_ITER = 200
DECIMAL_DIGITS = 4
ERROR_MESSAGE = "An Error Has Occurred"
epsilon = 0.001

class Cluster:
    def __init__(self, centroid:np.array):
        self.centroid = centroid
        self.prev_centroid = None
        self.data_points_list = []

    def assign_data_point(self, data_point):
        self.data_points_list.append(data_point)
    
    def clear_data_points_list(self):
        self.data_points_list.clear()

    def update_centroid(self):
        self.prev_centroid = self.centroid
        self.centroid = np.mean(self.data_points_list, axis=0)

    def print_centroid(self):
        d = len(self.centroid)
        for i in range(0, d):
            value = format(self.centroid[i], f'.{DECIMAL_DIGITS}f')
            if i < d-1:
                print(value, end=",")
            else:
                print(value)            

class KMeansAlgo:
    def __init__(self, K, N, d, iter, input_data):
        self.K = K
        self.N = N
        self.d = d
        self.iter = iter
        self.input_data = input_data

    def _has_clusters_converged(self, clusters):
        for cluster in clusters:
            if self._get_distance(cluster.prev_centroid, cluster.centroid) >= epsilon:
                return False
        return True


    def _get_distance(self, vector1, vector2):
        return np.linalg.norm(vector1-vector2)

    def _assign_cluster_to_datapoint(self, clusters, datapoints):
        for data_point in datapoints:
            distances = np.array([self._get_distance(data_point, cluster.centroid) for cluster in clusters])
            closest_cluster_idx = np.argmin(distances)
            clusters[closest_cluster_idx].assign_data_point(data_point)

    def _get_datapoints_from_file(self):
        with open(self.input_data, 'r') as file:
            datapoints = []
            counter = 0
            while True or counter<self.N:
                line = file.readline()
                if not line:
                    break
                point = np.fromstring(line, dtype=float, sep=',')
                datapoints.append(point)
                counter+= 1
            return datapoints

    def solve_kmeans_algo(self):
        data_points = self._get_datapoints_from_file()
        clusters = [Cluster(data_point) for data_point in data_points[:self.K]]
        for counter in range(0, self.iter):
            self._assign_cluster_to_datapoint(clusters, data_points)
            for cluster in clusters:
                cluster.update_centroid()
                cluster.clear_data_points_list()

            if self._has_clusters_converged(clusters):
                break
        return clusters


def print_clusters(clusters, d):
    for cluster in clusters:
        cluster.print_centroid()
        
def retrieve_argv_and_validate():
    if len(sys.argv) < 5:
        print(ERROR_MESSAGE)
        return None
    
    if not sys.argv[2].isnumeric() or int(sys.argv[2])<=1:
        print("Invalid number of points!")
        return None
    N = int(sys.argv[2])

    K = sys.argv[1]
    if not K.isnumeric() or int(K) <=1 or int(K) >= N:
        print("Invalid number of clusters!")
        return None
    K = int(K)
    
    d = sys.argv[3]
    if not d.isnumeric() or int(d)<1:
        print("Invalid dimension of point!")
    d = int(d)

    iter = DEAFULT_ITER
    fileName = ""

    if len(sys.argv) > 5:
        iter = sys.argv[4]
        if not iter.isnumeric() or int(iter) <= 1 or int(iter) >= 1000:
                print("Invalid maximum iteration!")
                return None
        iter = int(iter)
        fileName = sys.argv[5]    
    else:
        fileName = sys.argv[4]

    return K,N,d,iter,fileName

        
if __name__ == "__main__":
    parameters = retrieve_argv_and_validate()
    if parameters != None:
       kmeansAlgo = KMeansAlgo(*parameters)
       print_clusters(kmeansAlgo.solve_kmeans_algo(), parameters[2])
       


