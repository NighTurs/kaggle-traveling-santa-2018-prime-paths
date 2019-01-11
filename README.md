## 3rd palce solution to [Traveling Santa 2018 - Prime Paths](https://www.kaggle.com/c/traveling-santa-2018-prime-paths)

Code is a real mess, that is not how you should use branches... Still I don't think I will ever clean it, so will just give hints how to use this mess. Will do that by giving examples how to run it for every stage I listed on my [kaggle post](https://www.kaggle.com/c/traveling-santa-2018-prime-paths/discussion/77324).
You may have troubles running it on Windows because of gettimeofday (comment it out maybe..).

#### First stage:

I will just provide one of my .par files for LKH. Something like this with maybe different seed and multiple reruns and changes to PATCHING MOVE_TYPE params gave me my raw score.

#### Second stage:

1. Switch to master branch.
2. Make sure to have cities.csv and .cand file in the same directory.
3. Comment 'readTourSubmission(argv[2], tour)'. Oh shame on me... 
4. gcc -Wall -O3 -lm -lpthread main.c
5. Run with './a.out out9.proc submission.csv submission.csv 2 100 197769'. Arguments: raw TSP, submission.csv to continue from, output submission.csv, number of threads, timelimit seconds, CycleSize.
6. To restart you may want to undo step 3, and maybe change itK in code.

#### Third stage:

1. Switch to pure_gain branch (don't question the name).
2. Same as first stage
3. Compile as in a first stage
4. Run with './a.out out9.proc submission.csv submission.csv 2 10 1000 10'. Arguments same as in a first stage but last one is maximum K for opt to look into.
5. Here I would often find myself changing 'datas[i]->pureGainLimit = -1;' to someting like 0, -3. Shame on me for not having it as argument.

#### Fourth and fith stage:

1. Switch to gpx branch.
2. Apart from files for other stages you will need .tsp file here
3. Sudden switch to C++ to incorporate [GPX2 code](https://github.com/rtinos/gpx2), run make
4. Run with './a.out submission.csv submission.csv submission.csv 2 0.1 1000'. Instead of raw TSP for first argument here I have some other submission.csv I wanted to get additional candidates from.
5. Here I will change often pureLimit[] variable with is basically 'The positive gain criterion' for different K in opts. Also number of kicks in 'for (int kick = 0; kick < 11; kick++)'. And maxStep minStep at the same place to kick at smalled or larget chunks of tour.
 






 



