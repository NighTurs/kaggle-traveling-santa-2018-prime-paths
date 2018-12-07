#pragma clang diagnostic push
#pragma ide diagnostic ignored "cert-err34-c"

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdbool.h>
#include <stdlib.h>
#include <time.h>
#include <signal.h>

#define NUM_CITIES 197769
#define BUFF_SIZE 255
#define MAX_CAND 15
#define PENALTY 1.1
#define MAX_K 20
#define ILLEGAL_OPT -1e6
#define E 1e-9


typedef struct {
    double x, y;
} Point;

typedef struct {
    int city;
    double dist;
} Cand;

typedef struct NodeStruct {
    int cityId;
    struct NodeStruct *from, *to;
    int step, backStep;
    double cumFw[10], cumBack[10];
} Node;

typedef struct {
    int v;
    Node *n;
} PStruct;

typedef struct {
    Node (*nodes)[];
    Node *startNode;
    Node *t[MAX_K * 2 + 1];
    Node *tBest[MAX_K * 2 + 1];
    int incl[MAX_K * 2 + 1];
    int inclBest[MAX_K * 2 + 1];
    double dist[MAX_K * 2 + 1];
    PStruct p[MAX_K * 2 + 1];
    int q[MAX_K * 2 + 1];
    int cycle[MAX_K * 2 + 1];
    int size[MAX_K * 2 + 1];
    double maxGain;
    int maxK;
    int bestK;
    int bestRev;
    int minT;
    bool doReverse;
    bool isFindMax;
} KOptData;

int tour[NUM_CITIES];
int primes[NUM_CITIES];
Point cities[NUM_CITIES];
Node nodes[NUM_CITIES];
// At [i][0].city is size of candidates for city i
Cand candidates[NUM_CITIES][MAX_CAND + 1];
#define CAND_SIZE(x) (candidates[x][0].city)

void kOptStart(KOptData *data, int const order[], int orderSize);

void kOptMovRec(KOptData *data, int k);

double gainKOptMove(KOptData *data, int k);

void findPermutation(KOptData *data, int k);

double calcPenalty(Node *city, int step);

double cummDiff(KOptData *data, Node *a, Node *b, bool forward, int inStep);

bool isAdded(KOptData *data, Node *t1, Node *t2, int k);

bool isDeleted(KOptData *data, Node *t1, Node *t2, int k);

double getTourCost(Node nodes[]);

void calcNewNodes(Node *t[], const int incl[], int k, int revStart, Node oldNodes[], Node newNodes[], int tour[]);

void buildNodes(int const tour[], Node nodes[]);

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
    printf("%.5lf\n", getTourCost(nodes));

    KOptData data;
    data.nodes = &nodes;
    data.startNode = &nodes[0];
    data.maxK = 2;
    data.isFindMax = true;
    data.maxGain = 0;
    data.doReverse = true;
    int order[NUM_CITIES];
    for (int i = 0; i < NUM_CITIES; i++) {
        order[i] = i;
    }

    clock_t start, end;
    start = clock();
    kOptStart(&data, order, NUM_CITIES);
    end = clock();
    printf("%.5lf\n", ((double) (end - start)) / CLOCKS_PER_SEC);

    printf("%.5lf\n", data.maxGain);

    calcNewNodes(data.tBest, data.inclBest, data.bestK, data.bestRev, nodes, nodes, tour);

    printf("%.5lf\n", getTourCost(nodes));

    return 0;
}

void kOptStart(KOptData *data, int const order[], int orderSize) {
    for (int i = 0; i < orderSize; i++) {
        Node *t1 = &(*data->nodes)[order[i]];
        for (int x1 = 0; x1 < 2; x1++) {
            Node *t2 = x1 == 0 ? t1->from : t1->to;
            data->t[1] = t1;
            data->t[2] = t2;
            data->minT = t1->cityId;
            kOptMovRec(data, 2);
            if (!data->isFindMax && data->maxGain > E) {
                return;
            }
        }
    }
}

void kOptMovRec(KOptData *data, int k) {
    Node *t1 = data->t[1];
    Node *t2 = data->t[2 * k - 2];
    double d;
    for (int x3 = 1; x3 <= CAND_SIZE(t2->cityId); x3++) {
        Cand *t3Cand = &(candidates[t2->cityId][x3]);
        Node *t3 = &(*data->nodes)[t3Cand->city];
        if (t3->cityId < data->minT || t3 == t2->from || t3 == t2->to || isAdded(data, t2, t3, k - 2)) {
            continue;
        }
        data->t[2 * k - 1] = t3;
        data->incl[2 * k - 1] = 2 * k - 2;
        data->incl[2 * k - 2] = 2 * k - 1;
        data->dist[2 * k - 1] = t3Cand->dist;
        data->dist[2 * k - 2] = t3Cand->dist;
        for (int x4 = 0; x4 < 2; x4++) {
            Node *t4 = x4 == 0 ? t3->from : t3->to;
            if (t4->cityId < data->minT || isDeleted(data, t3, t4, k - 2)) {
                continue;
            }
            data->t[2 * k] = t4;
            data->incl[1] = 2 * k;
            data->incl[2 * k] = 1;
            d = distCity(t1->cityId, t4->cityId);
            data->dist[1] = d;
            data->dist[2 * k] = d;
            double gain = gainKOptMove(data, k);

            if (gain > data->maxGain) {
                data->maxGain = gain;
                memcpy(data->tBest + 1, data->t + 1, 2 * k * sizeof(Node *));
                memcpy(data->inclBest + 1, data->incl + 1, 2 * k * sizeof(int));
                data->bestK = k;
                data->bestRev = data->incl[0];
            }
            if (!data->isFindMax && data->maxGain > E) {
                return;
            }
        }
    }
}

double gainKOptMove(KOptData *data, int k) {
    findPermutation(data, k);

    int count = 1;
    int i = data->t[data->p[k * 2].v] != data->startNode ? 1 : k * 2;
    int startI = i;
    int endI = i == k * 2 ? 1 : k * 2;
    Node *prev = data->startNode;
    Node *cur = data->t[data->p[i].v];
    double cost = 0;
    int step = 0;
    bool forward = true;

    do {
        if (prev != data->startNode) {
            cost += cummDiff(data, prev, cur, forward, step);
        }
        step += abs(prev->step - cur->step) + 1;
        int nextT = data->incl[data->p[i].v];
        Node *next = data->t[nextT];
        cost += data->dist[nextT] * calcPenalty(cur, step);
        prev = next;
        i = data->q[nextT];
        forward = (i & 1) == 0;
        // odd -> -=1 even -> +=1
        i = i ^ 1;
        if (i > 2 * k || i < 1) {
            break;
        }
        count++;
        cur = data->t[data->p[i].v];
    } while (true);
    if (count != k) {
        return ILLEGAL_OPT;
    }

    double forwardGain =
            getTourCost(*data->nodes) - cost -
            cummDiff(data, data->t[data->p[endI].v], data->t[data->p[startI].v], forward, step);

    if (!data->doReverse) {
        data->incl[0] = 0;
        return forwardGain;
    }

    if (data->t[data->p[k * 2].v] == data->startNode) {
        i = 1;
    } else {
        i = k * 2;
    }
    startI = i;
    endI = i == k * 2 ? 1 : k * 2;
    cur = data->t[data->p[i].v];
    cost = 0;
    step = 0;
    prev = data->startNode;
    do {
        if (step != 0) {
            cost += cummDiff(data, prev, cur, forward, step);
        }
        if (step == 0) {
            if (startI == 1) {
                step += cur->step;
            } else {
                step += cur->backStep;
            }
        } else {
            step += abs(prev->backStep - cur->backStep) + 1;
        }
        int nextT = data->incl[data->p[i].v];
        Node *next = data->t[nextT];
        cost += data->dist[nextT] * calcPenalty(cur, step);
        prev = next;
        i = data->q[nextT];
        forward = (i & 1) == 0;
        // odd -> -=1 even -> +=1
        i = i ^ 1;
        if (i > 2 * k || i < 1) {
            break;
        }
        count++;
        cur = data->t[data->p[i].v];
    } while (true);

    double backwardGain =
            getTourCost(*data->nodes) - cost -
            cummDiff(data, data->t[data->p[endI].v], data->t[data->p[startI].v], forward, step);

    if (forwardGain > backwardGain) {
        data->incl[0] = 0;
        return forwardGain;
    }
    data->incl[0] = 1;
    return backwardGain;
}

bool isAdded(KOptData *data, Node *t1, Node *t2, int k) {
    int i = 2 * k;
    while ((i -= 2) > 0) {
        if ((t1 == data->t[i] && t2 == data->t[i + 1]) || (t1 == data->t[i + 1] && t2 == data->t[i])) {
            return true;
        }
    }
    return false;
}

bool isDeleted(KOptData *data, Node *t1, Node *t2, int k) {
    int i = 2 * k + 2;
    while ((i -= 2) > 0) {
        if ((t1 == data->t[i - 1] && t2 == data->t[i]) || (t1 == data->t[i] && t2 == data->t[i - 1])) {
            return true;
        }
    }
    return false;
}

int pCmp(const void *a, const void *b) {
    return (((PStruct *) a)->n->step - ((PStruct *) b)->n->step);
}

void findPermutation(KOptData *data, int k) {
    int i, j;
    // index from 1 to be consistent with paper
    for (i = j = 1; j <= k; i += 2, j++) {
        if (data->t[i]->from == data->t[i + 1]) {
            data->p[j].v = i + 1;
            data->p[j].n = data->t[i + 1];
        } else {
            data->p[j].v = i;
            data->p[j].n = data->t[i];
        }
    }
    // sorts by steps in increasing order
    qsort(data->p + 1, k, sizeof(PStruct), pCmp);

    for (j = 2 * k; j >= 2; j -= 2) {
        // from always on odd positions
        data->p[j - 1].v = i = data->p[j / 2].v;
        data->p[j].v = (i & 1) == 1 ? (i + 1) : (i - 1);
    }
    // q is used to lookup t_i position in p
    for (i = 1; i <= 2 * k; i++) {
        data->q[data->p[i].v] = i;
    }
}

double calcPenalty(Node *city, int step) {
    return (!primes[city->cityId] && step % 10 == 0) ? PENALTY : 1;
}

double cummDiff(KOptData *data, Node *a, Node *b, bool forward, int inStep) {
    if (a == b) {
        return 0;
    }
    int bStep = b->cityId == 0 ? NUM_CITIES + 1 : b->step;
    int aStep = a->cityId == 0 ? NUM_CITIES + 1 : a->step;
    if ((bStep - aStep > 0) == forward) {
        int stepDiff = abs(bStep - aStep);
        if (forward) {
            return b->cumFw[(inStep + stepDiff) % 10] - a->cumFw[inStep % 10];
        } else {
            if (a->cityId == 0) {
                return b->cumBack[(b->backStep - 1) % 10];
            } else {
                return b->cumBack[(inStep + stepDiff) % 10] - a->cumBack[inStep % 10];
            }
        }
    } else {
        Node *startNode = data->startNode;
        if (forward) {
            return b->cumFw[(bStep - 1) % 10] + startNode->cumFw[(NUM_CITIES) % 10] - a->cumFw[
                    (aStep - 1) % 10];
        } else {
            if (b->cityId == 0) {
                return startNode->cumBack[(startNode->backStep - 1) % 10] - a->cumBack[
                        (a->backStep - 1) % 10];
            } else {
                return b->cumBack[(b->backStep - 1) % 10] + startNode->cumBack[
                        (startNode->backStep - 1) % 10] - a->cumBack[(a->backStep - 1) % 10];
            }
        }
    }
}

void calcNewNodes(Node *t[], const int incl[], int k, int revStart, Node oldNodes[], Node newNodes[], int tour[]) {
    Node *cur = &oldNodes[0];
    int reverse = revStart;
    int c[k * 2 + 1];
    memset(c, 0, sizeof(int) * (k * 2 + 1));
    int pos = 0;
    if (reverse) {
        for (int i = 1; i <= k * 2; i += 2) {
            if (t[i]->cityId == 0) {
                if (cur->cityId != 0) {
                    exit(SIGABRT);
                }
                tour[pos++] = cur->cityId;
                if (cur->to == t[i + 1]) {
                    cur = cur->from;
                } else {
                    cur = cur->to;
                    reverse = false;
                }
            }
            if (t[i + 1]->cityId == 0) {
                if (cur->cityId != 0) {
                    exit(SIGABRT);
                }
                tour[pos++] = cur->cityId;
                if (cur->to == t[i]) {
                    cur = cur->from;
                } else {
                    cur = cur->to;
                    reverse = false;
                }
            }
        }
    }

    do {
        new_nodes_label:
        tour[pos++] = cur->cityId;
        for (int i = 2; i <= k * 2; i += 2) {
            if (c[i]) {
                continue;
            }
            Node *t1 = t[i];
            Node *t2 = t[i - 1];
            int nt = incl[i];
            Node *t3 = t[nt];
            Node *t4 = t[((nt - 1) ^ 1) + 1];
            if (cur != t1 && cur != t3) {
                continue;
            }
            c[i] = 1;
            if (cur == t1) {
                cur = t3;
                reverse = cur->to == t4;
            } else {
                cur = t1;
                reverse = cur->to == t2;
            }
            goto new_nodes_label;
        }
        if (reverse) {
            cur = cur->from;
        } else {
            cur = cur->to;
        }
    } while (cur->cityId != 0);
    buildNodes(tour, newNodes);
}

double getTourCost(Node nodes[]) {
    return nodes[0].cumFw[NUM_CITIES % 10];
}

void buildNodes(int const tour[], Node nodes[]) {
    int curCity = tour[1];
    int tourPos = 2;
    int prev = 0;
    memset(&nodes[0].cumFw, 0, sizeof(nodes[0].cumFw));
    memset(&nodes[0].cumBack, 0, sizeof(nodes[0].cumBack));

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
        Node *node = &nodes[curCity];
        node->cityId = curCity;
        node->from = &nodes[prev];
        node->to = &nodes[to];
        node->step = tourPos;
        Node *prevNode = &nodes[prev];
        double d = distCity(curCity, prev);
        for (int m = 0; m < 10; m++) {
            double cum = prevNode->cumFw[m == 0 ? 9 : m - 1];
            node->cumFw[m] = cum + ((m == 0 && !primes[prev]) ? d * PENALTY : d);
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
        Node *node = &nodes[curCity];
        node->backStep = step;
        Node *prevNode = &nodes[prev];
        double d = distCity(prev, curCity);
        for (int m = 0; m < 10; m++) {
            double cum = prevNode->cumBack[m == 0 ? 9 : m - 1];
            node->cumBack[m] = cum + ((m == 0 && !primes[prev]) ? d * PENALTY : d);
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