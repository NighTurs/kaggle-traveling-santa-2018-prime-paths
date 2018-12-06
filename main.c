#pragma clang diagnostic push
#pragma ide diagnostic ignored "cert-err34-c"

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdbool.h>

#define NUM_CITIES 197769
#define BUFF_SIZE 255
#define MAX_CAND 15


struct Point {
    double x, y;
};

struct Cand {
    int city;
    double dist;
};

struct Point cities[NUM_CITIES];
// At [i][0].city is size of candidates for city i
struct Cand candidates[NUM_CITIES][MAX_CAND + 1];
#define CAND_SIZE(x) (candidates[x][0].city)

void readCities(const char fileName[]);

void readTourLinkern(const char fileName[], int tour[]);

void buildCandidates(const char fileName[], int const tour[]);

double dist(double x1, double y1, double x2, double y2);

double distCity(int a, int b);

int main() {
    int tour[NUM_CITIES];
    readCities("cities.csv");
    readTourLinkern("out6.proc", tour);
    buildCandidates("my2.cand", tour);

    return 0;
}

double distCity(int a, int b) {
    return dist(cities[a].x, cities[a].y, cities[b].x, cities[b].y);
}

double dist(double x1, double y1, double x2, double y2) {
    return sqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2));
}

bool candContains(int a, int b) {
    for (int j = 1; j <= CAND_SIZE(a); j++) {
        if (candidates[a][j].city == b) {
            return true;
        }
    }
    return false;
}

void addCandidate(int a, int b) {
    candidates[a][0].city++;
    candidates[a][CAND_SIZE(a)].city = b;
    candidates[a][CAND_SIZE(a)].dist = distCity(a, b);
}

void buildCandidates(const char fileName[], int const tour[]) {
    FILE *fp;
    char buff[BUFF_SIZE];
    int cityA, cityB, num, pos, posOffset, alpha;

    fp = fopen(fileName, "r");
    // First line is NUM_OF_CITIES
    fgets(buff, BUFF_SIZE, fp);
    while (fgets(buff, BUFF_SIZE, fp)) {
        if (buff[0] == '-') {
            // Candidates are ending with -1 and then EOF
            break;
        }
        sscanf(buff, "%d %d %d %n", &cityA, &cityB, &num, &pos);
        cityA--;
        cityB--;
        candidates[cityA][0].city = num;
        for (int i = 1; i <= num; i++) {
            sscanf(buff + pos, "%d %d %n", &cityB, &alpha, &posOffset);
            cityB--;
            pos += posOffset;
            candidates[cityA][i].city = cityB;
            candidates[cityA][i].dist = distCity(cityA, cityB);
        }
    }
    for (int i = 0; i < NUM_CITIES; i++) {
        int from = tour[i];
        int to;
        if (i == NUM_CITIES - 1) {
            to = 0;
        } else {
            to = tour[i + 1];
        }
        if (!candContains(from, to)) {
            addCandidate(from, to);
        }
        if (!candContains(to, from)) {
            addCandidate(to, from);
        }
    }
    for (int i = 0; i < NUM_CITIES; i++) {
        for (int h = 1; h <= CAND_SIZE(i); h++) {
            int to = candidates[i][h].city;
            if (!candContains(to, i)) {
                addCandidate(to, i);
            }
        }
    }
    fclose(fp);
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

void readCities(const char fileName[]) {
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