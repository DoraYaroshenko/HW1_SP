#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define eps 0.001

typedef struct
{
    float *coordinates;
    int dimension;  // uint
} vector;

typedef struct
{
    vector *centroid;
    vector *members;
    int num_of_members; // uint
} cluster;

typedef struct {
    vector *all_vectors;
    int num_vectors;
} all_vecs;

int checkConvergence(vector *v1, vector *v2);
void assignVectorTocluster(vector *v, cluster *clus);
void updateCentroid(cluster *clus);
double distance(vector *v1, vector *v2);
vector *sumVectors(vector *vectors);
vector *mulByScalar(vector *v, double scalar);
void emptycluster(cluster *clus);
cluster *initiateClusters(vector *vec_array, int num_of_clusters);
cluster *iterateAlgorithm(cluster *cluster_array, vector *all_vectors, int K, int N, int iters);
all_vecs getInput();
void errorHandling();
void printOutput(cluster *clus, int K, int D);
void freeMemory(cluster *clus, vector *all_vectors, int K, int N);

int checkConvergence(vector *v1, vector *v2)
{
    return distance(v1, v2) < eps;
}

void assignVectorTocluster(vector *v, cluster *clus)
{
    clus->members = (vector *)realloc(clus->members, ((clus->num_of_members) + 1) * sizeof(vector));
    clus->num_of_members++;
    if (clus->members == NULL)
    {
        errorHandling();
    }
    clus->members[clus->num_of_members - 1] = *v;
}

void updateCentroid(cluster *clus)
{
    double scalar = 1 / (clus->num_of_members);
    vector *sum_vector = sumVectors(clus->members);
    clus->centroid = mulByScalar(sum_vector, scalar);
    // TODO: here you have a memory leak! 
    // sum_vector is allocated in sumVectors 
    // then passed into mulByScalar 
    // but mulByScalar allocates another vector for the result and 
    // returns a pointer to the new allocated result vector
    // then no one touches sum_vector again and the underlying memory will never be freed
}

double distance(vector *v1, vector *v2)
{
    double v1_coordinates[] = *(v1->coordinates);
    double v2_coordinates[] = *(v2->coordinates);
    double sum = 0;
    for (int i = 0; i < v1->dimension; i++)
    {
        sum += pow(v1_coordinates[i] - v2_coordinates[i], 2);
    }
    return sqrt(sum);
}

vector *sumVectors(vector *vectors)
{
    vector *sum_vec = malloc(sizeof(vector));
    sum_vec->dimension = vectors[0].dimension;
    sum_vec->coordinates = calloc(sum_vec->dimension, sizeof(float));
    if (sum_vec == NULL)
    {
        errorHandling();
    }
    int num_of_vectors = sizeof(vectors) / sizeof(vector); // TODO: you cant sizeof(vectors)
    for (int i = 0; i < sum_vec->dimension; i++)
    {
        for (int j = 0; j < num_of_vectors; j++)
        {
            sum_vec->coordinates[i] += vectors[j].coordinates[i];
        }
    }
    return sum_vec;
}

vector *mulByScalar(vector *v, double scalar)
{
    vector *mul_vec = malloc(sizeof(vector));
    mul_vec->dimension = v->dimension;
    mul_vec->coordinates = calloc(v->dimension, sizeof(float));
    if (mul_vec == NULL)
    {
        errorHandling();
    }
    for (int i = 0; i < mul_vec->dimension; i++)
    {
        mul_vec->coordinates[i] = v->coordinates[i] * scalar;
    }
    return mul_vec;
}

void emptycluster(cluster *clus)
{
    free(clus->members);
}

cluster *initiateClusters(vector *all_vectors, int K)
{
    cluster *cluster_array = (cluster *)malloc(K * sizeof(cluster));
    if (cluster_array == NULL)
    {
        errorHandling();
    }
    for (int i = 0; i < K; i++)
    {
        cluster_array[i].centroid = &all_vectors[i];
    }
    return cluster_array;
}
cluster *iterateAlgorithm(cluster *cluster_array, vector *all_vectors, int K, int N, int iter)
{
    for (int i = 0; i < iter; i++)
    {
        int convergence_flag = 0;
        for (int j = 0; j < N; j++)
        {
            double min_dist = distance(&all_vectors[j], cluster_array[0].centroid);
            cluster *potencial_cluster = &cluster_array[0];
            for (int k = 1; k < K; k++)
            {
                double new_dist = distance(&all_vectors[j], cluster_array[k].centroid);
                if (new_dist < min_dist)
                {
                    min_dist = new_dist;
                    potencial_cluster = &cluster_array[k];
                }
            }
            assignVectorTocluster(&all_vectors[j], potencial_cluster);
        }
        for (int k = 0; k < K; k++)
        {
            vector *old_centroid = cluster_array[k].centroid;
            updateCentroid(&cluster_array[k]);
            convergence_flag += checkConvergence(old_centroid, cluster_array[k].centroid);
            emptycluster(&cluster_array[k]);
        }
        if (convergence_flag == K)
            break;
    }
    return cluster_array;
}

all_vecs getInput()
{
    float n;
    char c;
    all_vecs all_vectors;
    all_vectors.all_vectors = (vector *)malloc(sizeof(vector));
    if (all_vectors.all_vectors == NULL)
    {
        errorHandling();
    }
    vector curr_vector;
    curr_vector.coordinates = (float *)malloc(sizeof(float));
    if (curr_vector.coordinates == NULL)
    {
        errorHandling();
    }
    int i = 0, j = 0;
    int vector_counter = 0;
    while (scanf("%lf%c", &n, &c) == 2)
    {
        int dimension_counter = 0;
        int coordinates_counter = 0;
        if (c == '\n')
        {
            curr_vector.coordinates[j] = n;
            all_vectors.all_vectors[i] = curr_vector;
            all_vectors.all_vectors = (vector *)realloc(all_vectors.all_vectors, sizeof(vector) * (vector_counter + 1));
            if (all_vectors.all_vectors == NULL)
            {
                errorHandling();
            }
            vector new_vector;
            new_vector.coordinates = (float *)malloc(sizeof(float));
            new_vector.dimension = dimension_counter;
            if (new_vector.coordinates == NULL)
            {
                errorHandling();
            }
            curr_vector = new_vector;
            vector_counter++;
            i++;
            continue;
        }
        curr_vector.coordinates[j] = n;
        curr_vector.coordinates = (float *)realloc(curr_vector.coordinates, sizeof(float) * (coordinates_counter + 1));
        if (curr_vector.coordinates == NULL)
        {
            errorHandling();
        }
        coordinates_counter++;
        dimension_counter++;
        j++;
    }
    // You dont update all_vectors.num_vectors! 
    return all_vectors;
}

void errorHandling()
{
    printf("An Error Has Occured");
    exit(1);
}

void printOutput(cluster *clus, int K, int d)
{
    for (int i = 0; i < K; i++)
    {
        for (int j = 0; j < d; j++)
        {
            printf("Float: %.4f, Char: %c", clus[i].centroid->coordinates[j], ',');
        }
        printf('\n');
    }
}

void freeMemory(cluster *clus, vector *all_vectors, int K, int N)
{
    for (int i = 0; i < N; i++)
    {
        free(all_vectors[i].coordinates);
    }
    free(all_vectors);
    for (int i = 0; i < K; i++)
    {
        free(clus->members); // remember to also free the cluster's centeroid 
    }
    free(clus);
}

int main(int argc, char **argv)
{
    all_vecs all_vectors = getInput();
    int dimension = all_vectors.all_vectors[0].dimension;
    int K = atoi(argv[1]);
    int N = all_vectors.num_vectors;
    if (!(K > 1 && K < N))
    {
        printf("Incorrect number of clusters!");
        exit(1);
    }
    int iter = 400;
    if (argc > 3)
    {
        iter = atoi(argv[2]);
    }
    if (!(iter > 1 && iter < 1000))
    {
        printf("Incorrect maximum iteration!");
        exit(1);
    }
    cluster *cluster_array = initiateClusters(&all_vectors, K);
    cluster_array = iterateAlgorithm(cluster_array, &all_vectors, K, N, iter);
    printOutput(cluster_array, K, dimension);
    freeMemory(cluster_array, &all_vectors, K, N);
}
