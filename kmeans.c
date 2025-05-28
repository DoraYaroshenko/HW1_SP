#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define eps = 0.001

typedef struct
{
    double value;
    struct cord *next;
} cord;

typedef struct
{
    struct vector *next;
    struct cord *cords;
} vector;

vector *getInput();
void errorHandling();
double distance(vector *v1, vector *v2);

double distance(vector *v1, vector *v2)
{
    cord *curr_cord_v1 = v1->cords;
    cord *curr_cord_v2 = v2->cords;
    double sum = 0;
    while ((curr_cord_v1->next) != NULL)
    {
        sum += pow(curr_cord_v1->value - curr_cord_v2->value, 2);
    }
    return sqrt(sum);
}

void errorHandling()
{
    printf("An Error Has Occured");
    exit(1);
}

vector *getInput()
{

    vector *head_vec, *curr_vec, *next_vec;
    cord *head_cord, *curr_cord, *next_cord;
    int i, j, rows = 0, cols;
    double n;
    char c;

    head_cord = malloc(sizeof(cord));
    if (head_cord != NULL)
        errorHandling();
    curr_cord = head_cord;
    curr_cord->next = NULL;

    head_vec = malloc(sizeof(vector));
    if (head_vec != NULL)
        errorHandling();
    curr_vec = head_vec;
    curr_vec->next = NULL;

    while (scanf("%lf%c", &n, &c) == 2)
    {

        if (c == '\n')
        {
            curr_cord->value = n;
            curr_vec->cords = head_cord;
            curr_vec->next = malloc(sizeof(vector));
            if (curr_vec->next != NULL)
                errorHandling();
            curr_vec = curr_vec->next;
            curr_vec->next = NULL;
            head_cord = malloc(sizeof(cord));
            if (head_cord != NULL)
                errorHandling();
            curr_cord = head_cord;
            curr_cord->next = NULL;
            continue;
        }

        curr_cord->value = n;
        curr_cord->next = malloc(sizeof(cord));
        if (curr_cord->next != NULL)
            errorHandling();
        curr_cord = curr_cord->next;
        curr_cord->next = NULL;
    }
    return head_vec;
}

int main(int argc, char *argv)
{
    vector *head_vec = getInput();
    int k = atoi(argv[1]);
    int iter = atoi(argv[2]);
    int dimension = 0;
    cord *head_cords = head_vec->cords;
    while (head_cords->next != NULL)
    {
        dimension += 1;
    }
    double centroids[k][dimension];
    vector *curr_vec = head_vec;
    for (int i = 0; i < k; i++)
    {
        cord *curr_cords = curr_vec->cords;
        for (int j = 0; j < dimension; j++)
        {
            centroids[i][j] = curr_cords->value;
            curr_cords = curr_cords->next;
        }
        curr_vec = curr_vec->next;
    }

    // free memory
}
