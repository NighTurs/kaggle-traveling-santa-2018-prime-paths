#pragma clang diagnostic push
#pragma ide diagnostic ignored "cert-err34-c"

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdbool.h>

#define NUM_CITIES 197769
#define BUFF_SIZE 255
#define MAX_CAND 15
#define PENALTY 1.1


struct Point {
    double x, y;
};

struct Cand {
    int city;
    double dist;
};

struct Node {
    int cityId;
    struct Node *from, *to;
    int step, backStep;
    double cumFw[10], cumBack[10];
};

int tour[NUM_CITIES];
int primes[NUM_CITIES];
struct Point cities[NUM_CITIES];
struct Node nodes[NUM_CITIES];
// At [i][0].city is size of candidates for city i
struct Cand candidates[NUM_CITIES][MAX_CAND + 1];
#define CAND_SIZE(x) (candidates[x][0].city)

double getTourCost(struct Node nodes[]);

void buildNodes(int const tour[], struct Node nodes[]);

void fillPrimes();

void readCities(const char fileName[]);

void readTourLinkern(const char fileName[], int tour[]);

void buildCandidates(const char fileName[], int const tour[]);

double dist(double x1, double y1, double x2, double y2);

double distCity(int a, int b);

int main() {
    readCities("cities.csv");
    readTourLinkern("out6.proc", tour);
    buildCandidates("my2.cand", tour);
    fillPrimes();
    buildNodes(tour, nodes);
    printf("%.5lf", getTourCost(nodes));

    return 0;
}

double getTourCost(struct Node nodes[]) {
    return nodes[0].cumFw[NUM_CITIES % 10];
}

void buildNodes(int const tour[], struct Node nodes[]) {
    int curCity = tour[1];
    int tourPos = 2;
    int prev = 0;
    bool fullLoop = false;
    do {
        int to;
        if (curCity == 0) {
            to = tour[1];
            fullLoop = true;
        } else if (tourPos != NUM_CITIES) {
            to = tour[tourPos];
        } else {
            to = 0;
        }
        struct Node *node = &nodes[curCity];
        node->cityId = curCity;
        node->from = &nodes[prev];
        node->to = &nodes[to];
        node->step = tourPos;
        struct Node *prevNode = &nodes[prev];
        double d = distCity(curCity, prev);
        for (int m = 0; m < 10; m++) {
            double cum = prevNode->cumFw[m == 0 ? 9 : m - 1];
            node->cumFw[m] += cum + ((m == 0 && !primes[prev]) ? d * PENALTY : d);
        }
        prev = curCity;
        curCity = to;
        tourPos++;
    } while (!fullLoop);

    nodes[0].step = 1;

    curCity = tour[NUM_CITIES - 1];
    prev = 0;
    int step = 2;
    do {
        struct Node *node = &nodes[curCity];
        node->backStep = step;
        struct Node *prevNode = &nodes[prev];
        double d = distCity(prev, curCity);
        for (int m = 0; m < 10; m++) {
            double cum = prevNode->cumBack[m == 0 ? 9 : m - 1];
            node->cumBack[m] += cum + ((m == 0 && !primes[prev]) ? d * PENALTY : d);
        }
        step++;
        prev = curCity;
        curCity = node->from->cityId;
    } while (curCity != tour[NUM_CITIES - 1]);
}

void fillPrimes() {
    memset(primes, 1, sizeof(primes));

    primes[0] = 0;
    primes[1] = 0;
    for (int i = 2; i < sqrt(NUM_CITIES) + 1; i++) {
        if (!primes[i]) {
            continue;
        }
        for (int h = i + i; h < NUM_CITIES; h += i) {
            primes[h] = 0;
        }
    }
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