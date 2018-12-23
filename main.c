#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wsign-conversion"
#pragma ide diagnostic ignored "cert-err34-c"

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdbool.h>
#include <stdlib.h>
#include <sys/time.h>
#include <signal.h>
#include <limits.h>
#include <pthread.h>

#define NUM_CITIES 197769
#define TIME_CHECK_FREQ 10000
#define BUFF_SIZE 255
#define WRITE_BUFF_SIZE 4096
#define MAX_CAND 15
#define PENALTY 1.1
#define MAX_K 20
#define ILLEGAL_OPT -1e6
#define E 1e-9
#define SHORT_CYCLE_SAFE 100
#define CYCLE_HIST_SIZE 30

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
    int city;
    double maxGain;
} CityGain;

typedef struct {
    Node *added1, *added2;
    Node *deleted1, *deleted2;
} NodeChange;

typedef struct {
    CityGain cityGain[NUM_CITIES];
    int *order;
    int orderSize;
    Node *nodes;
    Node *startNode;
    Node *t[MAX_K * 2 + 1];
    Node *tBest[MAX_K * 2 + 1];
    int incl[MAX_K * 2 + 1];
    int inclBest[MAX_K * 2 + 1];
    double dist[MAX_K * 2 + 1];
    PStruct p[MAX_K * 2 + 1];
    int q[MAX_K * 2 + 1];
    Node *deleted[MAX_K * 2 + 1];
    int nDeleted;
    Node *added[MAX_K * 2 + 1];
    int nAdded;
    double maxGain;
    int maxK;
    int cycleMax;
    int bestK;
    int bestRev;
    int bestCycle;
    int timeLimit;
    long secStart;
    long secEnd;
    long long timeChecks;
    int minT;
    int nProcessed;
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

void improveTour(int nThreads, const char subFile[], int timeLimit, int cycleLen);

void initCityGain(CityGain cityGain[]);

void *kOptStart(void *arg);

void kOptMovRec(KOptData *data, int k);

void updateCityGains(KOptData *data, double gain, int k);

void writeBestOpt(KOptData *data, double gain, int k, int cycle);

double gainKOptMove(KOptData *data, int k);

bool timeLimitExceeded(KOptData *data);

int patchCycles(KOptData *data, int k);

void patchCyclesRec(KOptData *data, int k, int m, int M, int curCycle, PStruct p[], int cycle[], int size[]);

int findCycle(KOptData *data, Node *N, int k, PStruct p[], int cycle[]);

int shortestCycle(KOptData *data, int M, int k, PStruct p[], const int cycle[], int size[]);

int countCycles(KOptData *data, PStruct p[], const int q[], int cycle[], int k);

void findPermutation(KOptData *data, PStruct p[], int q[], int k);

double calcPenalty(Node *city, int step);

double cummDiff(KOptData *data, Node *a, Node *b, bool forward, int inStep);

void markAdded(KOptData *data, Node *t1, Node *t2);

void unmarkAdded(KOptData *data, Node *t1, Node *t2);

bool isAdded(KOptData *data, Node *t1, Node *t2);

void markDeleted(KOptData *data, Node *t1, Node *t2);

void unmarkDeleted(KOptData *data, Node *t1, Node *t2);

bool isDeleted(KOptData *data, Node *t1, Node *t2);

double getTourCost(Node nodes[]);

void calcNewNodes(Node *t[], const int incl[], int k, int revStart, Node oldNodes[], Node newNodes[], int tour[]);

void buildNodes(int const tour[], Node nodes[]);

void fillPrimes();

void readCities(const char fileName[]);

void readTourLinkern(const char fileName[], int tour[]);

void readTourSubmission(const char fileName[], int tour[]);

void writeSubmission(const char fileName[], Node nodes[]);

void buildCandidates(const char fileName[], int const tour[]);

double dist(double x1, double y1, double x2, double y2);

double distCity(int a, int b);

void shuffle(int *array, size_t n);

//args: original_tour, submission_in, submission_out, num_threads, timeLimit, cycleLength
int main(int argc, char **argv) {
    srand(1234);
    readCities("cities.csv");
    readTourLinkern(argv[1], tour);
    buildCandidates("my2.cand", tour);
    readTourSubmission(argv[2], tour);
    fillPrimes();
    buildNodes(tour, nodes);
    printf("%.5lf\n", getTourCost(nodes));
    improveTour(atoi(argv[4]), argv[3], atoi(argv[5]), atoi(argv[6]));
    return 0;
}

int CityGainCmp(const void *a, const void *b) {
    return (((CityGain *) b)->maxGain - ((CityGain *) a)->maxGain) > 0 ? 1 : -1;
}

void improveTour(int nThreads, const char subFile[], int timeLimit, int cycleLen) {
    struct timeval start, end;
    struct timeval start2, end2;
    pthread_t thread_id[nThreads];
    KOptData *datas[nThreads];
    CityGain *cityGain = (CityGain *) malloc(sizeof(CityGain) * NUM_CITIES);
    int order[NUM_CITIES - 1];
    int cycleHist[CYCLE_HIST_SIZE];
    memset(cycleHist, 0, sizeof(cycleHist));
    int histIdx = 0;
    cycleHist[histIdx] = cycleLen;

    for (int i = 0; i < nThreads; i++) {
        datas[i] = (KOptData *) malloc(sizeof(KOptData));
        datas[i]->nodes = nodes;
        datas[i]->startNode = &nodes[0];
        datas[i]->isFindMax = true;
        datas[i]->nDeleted = 0;
        datas[i]->nAdded = 0;
    }

    KOptData *bestData = (KOptData *) malloc(sizeof(KOptData));
    int maxGainTh = 0;
    int iter = 0;
    do {
        double curCost = getTourCost(nodes);
        for (int i = 0; i < NUM_CITIES - 1; i++) {
            order[i] = i + 1;
        }
        shuffle(order, NUM_CITIES - 1);
        int orderSize = NUM_CITIES - 1;
        int curCycleLen = 0;
        for (int i = 0; i < CYCLE_HIST_SIZE; i++) {
            if (curCycleLen < cycleHist[i]) {
                curCycleLen = cycleHist[i];
            }
        }
        int chunk = (NUM_CITIES - 1 + nThreads) / nThreads;
        curCycleLen += SHORT_CYCLE_SAFE;
        int time = 0;
        int processed = 0;
        double maxGain = 0;
        int itK = 4;
        gettimeofday(&start, NULL);
        while (time < timeLimit) {
            for (int i = 0; i < nThreads; i++) {
                KOptData *cur = datas[i];
                cur->order = order + (i * chunk);
                cur->orderSize = (i + 1) * chunk > NUM_CITIES - 1 ? NUM_CITIES - 1 - (i * chunk) : chunk;
                cur->maxK = itK;
                cur->cycleMax = curCycleLen;
                cur->maxGain = 0;
                cur->timeLimit = (timeLimit - time) / 2 > 0 ? (timeLimit - time) / 2 : 1;
                cur->doReverse = iter < 10;
                if (cur->nAdded != 0 || cur->nDeleted != 0) {
                    printf("Added Deleted leak");
                    exit(1);
                }
                initCityGain(cur->cityGain);
                pthread_create(&thread_id[i], NULL, kOptStart, (void *) cur);
            }

            processed = 0;
            double maxCurGain = 0;
            gettimeofday(&start2, NULL);
            for (int i = 0; i < nThreads; i++) {
                pthread_join(thread_id[i], NULL);
                processed += datas[i]->nProcessed;
                if (maxCurGain < datas[i]->maxGain) {
                    maxCurGain = datas[i]->maxGain;
                    maxGainTh = i;
                }
            }
            if (maxCurGain > maxGain) {
                maxGain = maxCurGain;
                *bestData = *datas[maxGainTh];
            }

            gettimeofday(&end2, NULL);
            time += end2.tv_sec - start2.tv_sec;

            memcpy(cityGain, datas[0]->cityGain, sizeof(CityGain) * NUM_CITIES);
            for (int i = 0; i < nThreads; i++) {
                for (int h = 0; h < NUM_CITIES; h++) {
                    if (cityGain[h].maxGain < datas[i]->cityGain[h].maxGain) {
                        cityGain[h].maxGain = datas[i]->cityGain[h].maxGain;
                    }
                }
            }

            int inOrderSize = orderSize;
            orderSize = 0;
            qsort(cityGain, NUM_CITIES, sizeof(CityGain), CityGainCmp);
            for (int i = 0; i < NUM_CITIES; i++) {
                if (cityGain[i].maxGain > E) {
                    order[(i % nThreads) * chunk + i / nThreads] = cityGain[i].city;
                    orderSize = i + 1;
                } else {
                    break;
                }
            }

            printf("----Time=%ld Gain=%.3lf SolK=%d ItK=%d Cycle=%d MaxCycle=%d InCities=%d Processed=%d OutCities=%d\n",
                   end2.tv_sec - start2.tv_sec, maxCurGain,
                   datas[maxGainTh]->bestK,
                   itK, datas[maxGainTh]->bestCycle, curCycleLen, inOrderSize, processed, orderSize);
            fflush(stdout);
            itK++;
        }
        gettimeofday(&end, NULL);

        if (maxGain < E) {
            break;
        }

        if (bestData->bestCycle > 0) {
            histIdx++;
            if (histIdx >= CYCLE_HIST_SIZE) {
                histIdx = 0;
            }
            cycleHist[histIdx] = bestData->bestCycle;
        }

        calcNewNodes(bestData->tBest, bestData->inclBest, bestData->bestK, bestData->bestRev, nodes, nodes, tour);
        double newCost = getTourCost(nodes);
        printf("Iter=%d Time=%ld Gain=%.3lf GainDiff=%.3lf Cost=%.3lf K=%d CycleLen=%d MaxCycle=%d\n",
               iter,
               end.tv_sec - start.tv_sec,
               bestData->maxGain,
               curCost - newCost - bestData->maxGain,
               newCost,
               bestData->bestK,
               bestData->bestCycle,
               curCycleLen
        );
        fflush(stdout);
        writeSubmission(subFile, nodes);
        iter++;
    } while (true);

    for (int i = 0; i < nThreads; i++) {
        free(datas[i]);
    }
}

void initCityGain(CityGain cityGain[]) {
    for (int i = 0; i < NUM_CITIES; i++) {
        cityGain[i].maxGain = 0;
        cityGain[i].city = i;
    }
}

void *kOptStart(void *arg) {
    KOptData *data = (KOptData *) arg;
    Node **t = data->t;
    struct timeval start;
    gettimeofday(&start, NULL);
    data->secStart = start.tv_sec;
    data->timeChecks = 1;
    for (int i = 0; i < data->orderSize; i++) {
        Node *t1 = &data->nodes[data->order[i]];
        for (int x1 = 0; x1 < 2; x1++) {
            Node *t2 = x1 == 0 ? t1->to : t1->from;
            if (t2->cityId == 0) {
                continue;
            }
            t[1] = t1;
            t[2] = t2;
            data->minT = t1->cityId;
            markDeleted(data, t1, t2);
            kOptMovRec(data, 2);
            unmarkDeleted(data, t1, t2);
            if (!data->isFindMax && data->maxGain > E) {
                data->nProcessed = i + 1;
                return NULL;
            }
            if (timeLimitExceeded(data)) {
                data->nProcessed = i + 1;
                return NULL;
            }
        }
    }
    data->nProcessed = data->orderSize;
    return NULL;
}

void kOptMovRec(KOptData *data, int k) {
    Node **t = data->t;
    int *incl = data->incl;
    double *dist = data->dist;
    Node *nodes = data->nodes;
    int minT = data->minT;
    bool isFindMax = data->isFindMax;
    int maxK = data->maxK;
    Node *t1 = t[1];
    Node *t2 = t[2 * k - 2];
    double d;
    for (int x3 = 1; x3 <= CAND_SIZE(t2->cityId); x3++) {
        Cand *t3Cand = &(candidates[t2->cityId][x3]);
        Node *t3 = &nodes[t3Cand->city];
        if (t3->cityId < minT || t3 == t2->from || t3 == t2->to || isAdded(data, t2, t3)) {
            continue;
        }
        t[2 * k - 1] = t3;
        incl[2 * k - 1] = 2 * k - 2;
        incl[2 * k - 2] = 2 * k - 1;
        dist[2 * k - 1] = t3Cand->dist;
        dist[2 * k - 2] = t3Cand->dist;
        markAdded(data, t2, t3);
        for (int x4 = 0; x4 < 2; x4++) {
            Node *t4 = x4 == 0 ? t3->to : t3->from;
            if (t4->cityId < minT || isDeleted(data, t3, t4)) {
                continue;
            }
            t[2 * k] = t4;
            incl[1] = 2 * k;
            incl[2 * k] = 1;
            d = distCity(t1->cityId, t4->cityId);
            dist[1] = d;
            dist[2 * k] = d;
            markDeleted(data, t3, t4);
            if (!isDeleted(data, t1, t4)) {
                double gain = gainKOptMove(data, k);

                if (gain > E) {
                    updateCityGains(data, gain, k);
                }
                if (gain > data->maxGain) {
                    writeBestOpt(data, gain, k, 0);
                }
                if (!isFindMax && data->maxGain > E) {
                    unmarkAdded(data, t2, t3);
                    unmarkDeleted(data, t3, t4);
                    return;
                }
                if (gain == ILLEGAL_OPT) {
                    if (k + 2 <= maxK && t4 != t1) {
                        markAdded(data, t4, t1);
                        patchCycles(data, k);
                        unmarkAdded(data, t4, t1);
                        if (!data->isFindMax && data->maxGain > E) {
                            unmarkAdded(data, t2, t3);
                            unmarkDeleted(data, t3, t4);
                            return;
                        }
                    }
                }
            }
            if (k < maxK) {
                kOptMovRec(data, k + 1);
                if (!isFindMax && data->maxGain > E) {
                    unmarkAdded(data, t2, t3);
                    unmarkDeleted(data, t3, t4);
                    return;
                }
            }
            if (timeLimitExceeded(data)) {
                unmarkAdded(data, t2, t3);
                unmarkDeleted(data, t3, t4);
                return;
            }
            unmarkDeleted(data, t3, t4);
        }
        unmarkAdded(data, t2, t3);
    }
}

int patchCycles(KOptData *data, int k) {
    PStruct p[k * 2 + 1];
    int q[k * 2 + 1];
    int cycle[k * 2 + 1];
    int size[k * 2 + 1];
    findPermutation(data, p, q, k);

    int M, i;
    Node *s2;

    M = countCycles(data, p, q, cycle, k);

    if (data->maxK < k + M) {
        return M;
    }
    int curCycle = shortestCycle(data, M, k, p, cycle, size);
    if (data->cycleMax < size[curCycle]) {
        return M;
    }

    for (i = 0; i < k; i++) {
        if (cycle[p[2 * i].v] != curCycle) {
            continue;
        }
        Node *sStart = data->t[p[2 * i].v];
        Node *sStop = data->t[p[2 * i + 1].v];
        for (Node *s1 = sStart; s1 != sStop; s1 = s2) {
            s2 = s1->to;
            if (s1->cityId == 0 || s2->cityId == 0) {
                continue;
            }
            data->t[2 * k + 1] = s1;
            data->t[2 * k + 2] = s2;

            markDeleted(data, s1, s2);
            patchCyclesRec(data, k, 2, M, curCycle, p, cycle, size);
            unmarkDeleted(data, s1, s2);
            if (!data->isFindMax && data->maxGain > E) {
                return M;
            }
            if (timeLimitExceeded(data)) {
                return M;
            }
        }
    }
    return M;
}

void patchCyclesRec(KOptData *data, int k, int m, int M, int curCycle, PStruct p[], int cycle[], int size[]) {
    Node **t = data->t;
    int *incl = data->incl;
    double *dist = data->dist;
    double d;
    int newCycle, i;
    int cycleAdj[1 + 2 * k];
    Node *s1 = t[2 * k + 1];
    Node *s2 = t[i = 2 * (k + m) - 2];
    incl[incl[i] = i + 1] = i;
    for (int x3 = 1; x3 <= CAND_SIZE(s2->cityId); x3++) {
        Cand *s3Cand = &(candidates[s2->cityId][x3]);
        Node *s3 = &data->nodes[s3Cand->city];
        if (s3->cityId == 0 || s3 == s2->from || s3 == s2->to || isAdded(data, s2, s3)
            || (newCycle = findCycle(data, s3, k, p, cycle)) == curCycle) {
            continue;
        }
        t[2 * (k + m) - 1] = s3;
        dist[2 * (k + m) - 2] = s3Cand->dist;
        dist[2 * (k + m) - 1] = s3Cand->dist;
        markAdded(data, s2, s3);
        for (int x4 = 0; x4 < 2; x4++) {
            Node *s4 = x4 == 0 ? s3->to : s3->from;
            if (s4->cityId == 0 || isDeleted(data, s3, s4)) {
                continue;
            }
            t[2 * (k + m)] = s4;
            d = distCity(s4->cityId, s1->cityId);
            markDeleted(data, s3, s4);
            if (M > 2) {
                for (i = 1; i <= 2 * k; i++) {
                    cycleAdj[i] = cycle[i] == newCycle ? curCycle : cycle[i];
                }
                patchCyclesRec(data, k, m + 1, M - 1, curCycle, p, cycleAdj, size);
                if (!data->isFindMax && data->maxGain > E) {
                    unmarkAdded(data, s2, s3);
                    unmarkDeleted(data, s3, s4);
                    return;
                }
                if (s4 == s1) {
                    unmarkDeleted(data, s3, s4);
                    continue;
                }
                incl[incl[2 * k + 1] = 2 * (k + m)] = 2 * k + 1;
                dist[2 * k + 1] = d;
                dist[2 * (k + m)] = d;
                markAdded(data, s4, s1);
                patchCycles(data, k + m);
                unmarkAdded(data, s4, s1);
                if (!data->isFindMax && data->maxGain > E) {
                    unmarkAdded(data, s2, s3);
                    unmarkDeleted(data, s3, s4);
                    return;
                }
            } else if (s4 != s1) {
                incl[incl[2 * k + 1] = 2 * (k + m)] = 2 * k + 1;
                dist[2 * k + 1] = d;
                dist[2 * (k + m)] = d;
                double gain = gainKOptMove(data, k + m);
                if (gain > E) {
                    updateCityGains(data, gain, k + m);
                }
                if (gain > data->maxGain) {
                    writeBestOpt(data, gain, k + m, size[curCycle]);
                }
                if (!data->isFindMax && data->maxGain > E) {
                    unmarkAdded(data, s2, s3);
                    unmarkDeleted(data, s3, s4);
                    return;
                }
            }
            unmarkDeleted(data, s3, s4);
        }
        unmarkAdded(data, s2, s3);
    }

    if (M != 2 || data->maxK < k + m + 1) {
        return;
    }

    for (int x3 = 1; x3 <= CAND_SIZE(s2->cityId); x3++) {
        Cand *s3Cand = &(candidates[s2->cityId][x3]);
        Node *s3 = &data->nodes[s3Cand->city];
        if (s3->cityId == 0 || s3 == s2->from || s3 == s2->to || isAdded(data, s2, s3)) {
            continue;
        }
        t[2 * (k + m) - 1] = s3;
        dist[2 * (k + m) - 2] = s3Cand->dist;
        dist[2 * (k + m) - 1] = s3Cand->dist;
        for (int x4 = 0; x4 < 2; x4++) {
            Node *s4 = x4 == 0 ? s3->to : s3->from;
            if (s4->cityId == 0 || isDeleted(data, s3, s4)) {
                continue;
            }
            t[2 * (k + m)] = s4;

            for (int x5 = 1; x5 <= CAND_SIZE(s4->cityId); x5++) {
                Cand *s5Cand = &(candidates[s4->cityId][x5]);
                Node *s5 = &data->nodes[s5Cand->city];
                if (s5->cityId == 0 || s5 == s4->from || s5 == s4->to || isAdded(data, s4, s5) ||
                    (newCycle = findCycle(data, s5, k, p, cycle)) == curCycle) {
                    continue;
                }

                t[2 * (k + m) + 1] = s5;
                incl[incl[2 * (k + m)] = 2 * (k + m) + 1] = 2 * (k + m);
                dist[2 * (k + m)] = s5Cand->dist;
                dist[2 * (k + m) + 1] = s5Cand->dist;
                for (int x6 = 0; x6 < 2; x6++) {
                    Node *s6 = x6 == 0 ? s5->to : s5->from;
                    if (s6->cityId == 0 || isDeleted(data, s5, s6) || isAdded(data, s6, s1)) {
                        continue;
                    }
                    t[2 * (k + m) + 2] = s6;
                    incl[incl[2 * k + 1] = 2 * (k + m) + 2] = 2 * k + 1;
                    d = distCity(s6->cityId, s1->cityId);
                    dist[2 * k + 1] = d;
                    dist[2 * (k + m) + 2] = d;
                    double gain = gainKOptMove(data, k + m + 1);
                    if (gain > E) {
                        updateCityGains(data, gain, k + m + 1);
                    }
                    if (gain > data->maxGain) {
                        writeBestOpt(data, gain, k + m + 1, size[curCycle]);
                    }
                    if (!data->isFindMax && data->maxGain > E) {
                        return;
                    }
                }
            }
        }
    }
}

void updateCityGains(KOptData *data, double gain, int k) {
    CityGain *cityGain = data->cityGain;
    Node **t = data->t;
    for (int i = 1; i <= k * 2; i++) {
        if (cityGain[t[i]->cityId].maxGain < gain) {
            cityGain[t[i]->cityId].maxGain = gain;
        }
    }
}

void writeBestOpt(KOptData *data, double gain, int k, int cycle) {
    data->maxGain = gain;
    memcpy(data->tBest + 1, data->t + 1, 2 * k * sizeof(Node *));
    memcpy(data->inclBest + 1, data->incl + 1, 2 * k * sizeof(int));
    data->bestK = k;
    data->bestRev = data->incl[0];
    data->bestCycle = cycle;
}

bool isBetween(Node *a, Node *b, Node *c) {
    return b->step <= c->step && b->step >= a->step;
}

int findCycle(KOptData *data, Node *N, int k, PStruct p[], int cycle[]) {
    /* Binary search */
    int low = 1, high = k;
    while (low < high) {
        int mid = (low + high) / 2;
        if (isBetween(data->t[p[2 * low].v], N, data->t[p[2 * mid + 1].v])) {
            high = mid;
        } else {
            low = mid + 1;
        }
    }
    return cycle[p[2 * low].v];
}

int countCycles(KOptData *data, PStruct p[], const int q[], int cycle[], int k) {
    int i, j, count = 0;

    for (i = 1; i <= 2 * k; i++) {
        cycle[i] = 0;
    }
    for (i = 1; i <= 2 * k; i++) {
        if (cycle[p[i].v] == 0) {
            count++;
            j = i;
            do {
                cycle[p[j].v] = count;
                j = q[data->incl[p[j].v]];
                cycle[p[j].v] = count;
                if ((j ^= 1) > 2 * k) {
                    j = 1;
                }
            } while (j != i);
        }
    }
    return count;
}

int shortestCycle(KOptData *data, int M, int k, PStruct p[], const int cycle[], int size[]) {
    int i, minCycle = -1, minSize = INT_MAX;
    for (i = 1; i <= M; i++) {
        size[i] = 0;
    }
    p[0].v = p[2 * k].v;
    for (i = 0; i < 2 * k; i += 2) {
        size[cycle[p[i].v]] += abs(data->t[p[i].v]->step - data->t[p[i + 1].v]->step) + 1;
    }
    for (i = 1; i <= M; i++) {
        if (size[i] < minSize) {
            minSize = size[i];
            minCycle = i;
        }
    }
    return minCycle;
}

bool timeLimitExceeded(KOptData *data) {
    if (data->timeLimit == 0) {
        return false;
    }
    data->timeChecks++;
    if (data->timeChecks % TIME_CHECK_FREQ == 0) {
        struct timeval t;
        gettimeofday(&t, NULL);
        data->secEnd = t.tv_sec;
    }
    if (data->secEnd - data->secStart >= data->timeLimit) {
        return true;
    } else {
        return false;
    }
}

double gainKOptMove(KOptData *data, int k) {
    Node **t = data->t;
    int *incl = data->incl;
    double *dist = data->dist;
    PStruct *p = data->p;
    int *q = data->q;
    Node *startNode = data->startNode;

    findPermutation(data, p, q, k);

    int count = 1;
    int i = 1;
    int startI = i;
    int endI = k * 2;
    Node *prev = startNode;
    Node *cur = t[p[i].v];
    double cost = 0;
    int step = 0;
    bool forward = true;

    do {
        if (prev != startNode) {
            cost += cummDiff(data, prev, cur, forward, step);
        }
        step += abs(prev->step - cur->step) + 1;
        int nextT = incl[p[i].v];
        Node *next = t[nextT];
        cost += dist[nextT] * calcPenalty(cur, step);
        prev = next;
        i = q[nextT];
        forward = (i & 1) == 0;
        // odd -> -=1 even -> +=1
        i = i ^ 1;
        if (i > 2 * k || i < 1) {
            break;
        }
        count++;
        cur = t[p[i].v];
    } while (true);
    if (count != k) {
        return ILLEGAL_OPT;
    }

    double forwardGain =
            getTourCost(data->nodes) - cost -
            cummDiff(data, t[p[endI].v], t[p[startI].v], forward, step);

    if (!data->doReverse) {
        incl[0] = 0;
        return forwardGain;
    }


    i = k * 2;
    startI = i;
    endI = 1;
    cur = t[p[i].v];
    cost = 0;
    step = 0;
    prev = startNode;
    do {
        if (step != 0) {
            cost += cummDiff(data, prev, cur, forward, step);
        }
        if (step == 0) {
            step += cur->backStep;
        } else {
            step += abs(prev->backStep - cur->backStep) + 1;
        }
        int nextT = incl[p[i].v];
        Node *next = t[nextT];
        cost += dist[nextT] * calcPenalty(cur, step);
        prev = next;
        i = q[nextT];
        forward = (i & 1) == 0;
        // odd -> -=1 even -> +=1
        i = i ^ 1;
        if (i > 2 * k || i < 1) {
            break;
        }
        count++;
        cur = t[p[i].v];
    } while (true);

    double backwardGain =
            getTourCost(data->nodes) - cost -
            cummDiff(data, t[p[endI].v], t[p[startI].v], forward, step);

    if (forwardGain > backwardGain) {
        incl[0] = 0;
        return forwardGain;
    }
    incl[0] = 1;
    return backwardGain;
}

void markAdded(KOptData *data, Node *t1, Node *t2) {
    data->nAdded++;
    data->added[data->nAdded * 2 - 1] = t1;
    data->added[data->nAdded * 2] = t2;
}

void unmarkAdded(KOptData *data, Node *t1, Node *t2) {
    data->nAdded--;
}


bool isAdded(KOptData *data, Node *t1, Node *t2) {
    Node **added = data->added;
    for (int i = 1; i <= data->nAdded; i++) {
        if ((added[i * 2 - 1] == t1 && added[i * 2] == t2) || (added[i * 2 - 1] == t2 && added[i * 2] == t1)) {
            return true;
        }
    }
    return false;
}

void markDeleted(KOptData *data, Node *t1, Node *t2) {
    data->nDeleted++;
    data->deleted[data->nDeleted * 2 - 1] = t1;
    data->deleted[data->nDeleted * 2] = t2;
}

void unmarkDeleted(KOptData *data, Node *t1, Node *t2) {
    data->nDeleted--;
}

bool isDeleted(KOptData *data, Node *t1, Node *t2) {
    Node **deleted = data->deleted;
    for (int i = 1; i <= data->nDeleted; i++) {
        if ((deleted[i * 2 - 1] == t1 && deleted[i * 2] == t2) || (deleted[i * 2 - 1] == t2 && deleted[i * 2] == t1)) {
            return true;
        }
    }
    return false;
}

int pCmp(const void *a, const void *b) {
    return (((PStruct *) a)->n->step - ((PStruct *) b)->n->step);
}

void findPermutation(KOptData *data, PStruct p[], int q[], int k) {
    Node **t = data->t;
    int i, j;
    // index from 1 to be consistent with paper
    for (i = j = 1; j <= k; i += 2, j++) {
        if (t[i]->from == t[i + 1]) {
            p[j].v = i + 1;
            p[j].n = t[i + 1];
        } else {
            p[j].v = i;
            p[j].n = t[i];
        }
    }
    // sorts by steps in increasing order
    qsort(p + 1, k, sizeof(PStruct), pCmp);

    for (j = 2 * k; j >= 2; j -= 2) {
        // from always on odd positions
        p[j - 1].v = i = p[j / 2].v;
        p[j].v = (i & 1) == 1 ? (i + 1) : (i - 1);
    }
    // q is used to lookup t_i position in p
    for (i = 1; i <= 2 * k; i++) {
        q[p[i].v] = i;
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

void readTourSubmission(const char fileName[], int tour[]) {
    FILE *fp;
    char buff[BUFF_SIZE];
    int city;

    fp = fopen(fileName, "r");
    int i = 0;
    fgets(buff, BUFF_SIZE, fp);
    while (fgets(buff, BUFF_SIZE, fp)) {
        sscanf(buff, "%d", &city);
        if (i < NUM_CITIES) {
            tour[i++] = city;
        }
    }
    fclose(fp);
}

void writeSubmission(const char fileName[], Node nodes[]) {
    FILE *fp = fopen(fileName, "w");
    char file_buffer[WRITE_BUFF_SIZE + 64];
    int buffer_count = 0;

    Node *node = &nodes[0];
    buffer_count += sprintf(&file_buffer[buffer_count], "Path\n");
    do {
        buffer_count += sprintf(&file_buffer[buffer_count], "%d\n", node->cityId);
        node = node->to;
        if (buffer_count >= WRITE_BUFF_SIZE) {
            fwrite(file_buffer, WRITE_BUFF_SIZE, 1, fp);
            buffer_count -= WRITE_BUFF_SIZE;
            memcpy(file_buffer, &file_buffer[WRITE_BUFF_SIZE], buffer_count);
        }
    } while (node->cityId != 0);
    buffer_count += sprintf(&file_buffer[buffer_count], "0");

    if (buffer_count > 0) {
        fwrite(file_buffer, 1, buffer_count, fp);
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

void shuffle(int *array, size_t n) {
    if (n > 1) {
        size_t i;
        for (i = 0; i < n - 1; i++) {
            size_t j = i + rand() / (RAND_MAX / (n - i) + 1);
            int t = array[j];
            array[j] = array[i];
            array[i] = t;
        }
    }
}

#pragma clang diagnostic pop