#pragma clang diagnostic push
#pragma ide diagnostic ignored "cert-err34-c"
#include <stdio.h>
#include <string.h>

#define NUM_CITIES 197769
#define BUFF_SIZE 255

struct Point {
    double x, y;
};

void readCities(const char fileName[], struct Point cities[]);
void readTourLinkern(const char fileName[], int tour[]);

int main() {

    struct Point cities[NUM_CITIES];
    int tour[NUM_CITIES];
    readCities("cities.csv", cities);
    readTourLinkern("out6.proc", tour);

    return 0;
}

void readTourLinkern(const char fileName[], int tour[]) {
    FILE *fp;
    char buff[BUFF_SIZE];
    int cityA, cityB, dist;

    fp = fopen(fileName, "r");
    int i = 0;
    while (fgets(buff, BUFF_SIZE, fp)) {
        sscanf(buff, "%d %d %d", &cityA, &cityB, &dist);
        tour[i++] = cityA;
    }
    fclose(fp);
}

void readCities(const char fileName[], struct Point cities[]) {
    FILE *fp;
    char buff[BUFF_SIZE];
    int id;
    double x, y;

    fp = fopen(fileName, "r");
    // First line is header
    fgets(buff, BUFF_SIZE, fp);
    while (fgets(buff, BUFF_SIZE, fp)) {
        sscanf(buff, "%d,%lf,%lf", &id, &x, &y);
        cities[id].x = x;
        cities[id].y = y;
    }
    fclose(fp);
}


#pragma clang diagnostic pop