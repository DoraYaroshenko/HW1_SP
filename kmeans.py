import math
import sys

def euclidian_distance(v1, v2):
    #להוסיף בדיקה לאורכי הוקטורים
    sum = 0
    for i in range(len(v1)):
       sub = v1[i] - v2[i]
       sum += sub ** 2
    return math.sqrt(sum)

if __name__ == "__main__":

    K=7
    iter = 400
    epsilon = 0.001
    #try:
        #K = (int)(sys.argv(1))
    #except:
        #print("Incorrect number or clusters!")
        #sys.exit(1)
    #if(K <= 1):
        #print("Incorrect number of clusters!")
        #sys.exit(1)
    
    #try:
        #iter = sys.argv(2)
    #finally:
        #try:
            #iter = (int)(iter)
        #except:
            #print("Incorrect maximum iteration!")
            #sys.exit(1)
    
    #if(iter>=1000 or iter<=1):
        #print("Incorrect maximum iteration!")
        #sys.exit(1)

    file = open('.vscode/tests/input_2.txt', 'r')
    content = file.read()
    vectors = content.split("\n")
    vectors = vectors[:len(vectors)-1]
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
            #mean = 0
            mean = [0.0 for x in range(len(centroids[j]))]
            #sum = 0
            #count = 0
            for assigned in assigned_to_cent[j]:
                for p in range(len(assigned)):
                    mean[p] += assigned[p]
                #d = euclidian_distance(centroids[j], assigned)
                #sum+=d
                #count+=1
            #mean = sum / count

            for p in range(len(mean)):
                mean[p] /= len(assigned_to_cent[j])

            if (euclidian_distance(centroids[j], mean))>=epsilon:
                convergence = False
            centroids[j] = mean
        if convergence:
            break

    for c in centroids:
        print(','.join(f"{x:.4f}" for x in c))


    #i = 0
    #for v in centroids:
        #i+=1
        #print(v)
        #if(i>10):
            #break