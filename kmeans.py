import math
import sys

def euclidian_distance(v1, v2):
    sum = 0
    for i in range(len(v1)):
       sub = v1[i] - v2[i]
       sum += sub ** 2
    return math.sqrt(sum)

if __name__ == "__main__":

    K=-1
    iter = 400
    filename = ""
    epsilon = 0.001

    try:
        K = float(sys.argv[1])
        if not K.is_integer():
            raise ValueError
        K = int(K)
    except:
        print("Incorrect number of clusters!")
        sys.exit(1) 
    if(K <= 1): 
        print("Incorrect number of clusters!")
        sys.exit(1)

    content = sys.stdin.read()
    vectors = content.split("\n")
    vectors = vectors[:len(vectors)-1]

    N = len(vectors)
    if (K>=N):
        print("Incorrect number of clusters!")
        sys.exit(1)

    if len(sys.argv) > 2:
        try:
            iter = float(sys.argv[2])
            if not iter.is_integer():
                raise ValueError
            iter = int(iter)
        except:
            print("Incorrect maximum iteration!")
            sys.exit(1)
        if(iter>=1000 or iter<=1):
            print("Incorrect maximum iteration!")
            sys.exit(1)

    for i in range(len(vectors)):
        vectors[i] = [float(x) for x in vectors[i].split(",")]
    
    centroids = []
    for i in range(K):
        centroids.append(vectors[i].copy())
    
    for i in range(iter):
        assigned_to_cent = [ [] for _ in range(len(centroids)) ]
        convergence = True
        for v in vectors:
            mn = sys.float_info.max
            mn_j = 0
            for j in range(K):
                c = centroids[j]
                d = euclidian_distance(v, c)
                if (d<mn):
                    mn = d
                    mn_j = j
            assigned_to_cent[mn_j].append(v)
    
        for j in range(len(centroids)):
            if not assigned_to_cent[j]:
                continue
            mean = [0.0 for x in range(len(centroids[j]))]
            for assigned in assigned_to_cent[j]:
                for p in range(len(assigned)):
                    mean[p] += assigned[p]

            for p in range(len(mean)):
                mean[p] /= len(assigned_to_cent[j])

            if (euclidian_distance(centroids[j], mean))>=epsilon:
                convergence = False
            centroids[j] = mean
        if convergence:
            break

    for c in centroids:
        print(','.join(f"{x:.4f}" for x in c))
