#include <stdio.h>
#include <pthread.h>
#include <math.h>
#include <stdlib.h>
#include <sys/time.h>

//! Uncomment if running on mac
#include "../pthread_barrier.h"

// Parameters
int boardSize;
double precision;
int nThreads;

// Matrices
double **board, **results;
// Arrays for coordination and communication.
int *rowStart, *rowEnd, *threadDone;
// Signal for threads to run or not
int run;

pthread_barrier_t barrier;

int printBoard() {
    for (int i = 0; i < boardSize; ++i) {
        for (int j = 0; j < boardSize; ++j) {
            printf("%.4f ", board[i][j]);
        }
        printf("\n\n");
    }
    return 0;
}

void *threadFunction(void *arg) {
    const int id = *(int *) arg;
//    Start and end indices for this thread.
    const int start = rowStart[id];
    const int end = rowEnd[id];
    const int target = (end - start) * (boardSize - 2);

    int nCorrect;
    double newValue, currentValue, diff;
    while (run == 1) {
        nCorrect = 0;
        for (int i = start; i < end; ++i) {
            for (int j = 1; j < boardSize - 1; ++j) {
//                Save old value
                currentValue = board[i][j];

//                Make and save average
                newValue = (board[i][j + 1] +
                            board[i + 1][j] +
                            board[i][j - 1] +
                            board[i - 1][j]) / 4;
                results[i][j] = newValue;

//                Mark if value was correct
                diff = fabs(newValue - currentValue);
                if (diff < precision) { ++nCorrect; }
            }
        }

//        Mark all values were correct.
        threadDone[id] = (nCorrect == target);

//        Wait for other threads to finish
//        printf("%d: Done. Waiting for first barrier...\n", id);
        pthread_barrier_wait(&barrier);

//        Wait for main thread to check results
//        printf("%d: Waiting for main to check...\n", id);
        pthread_barrier_wait(&barrier);
    };
    return NULL;
}

int checkBoard() {
    double result;
    double difference;
    for (int i = 1; i < boardSize - 1; ++i) {
        for (int j = 1; j < boardSize - 1; ++j) {
            result = (
                             board[i][j + 1] +
                             board[i + 1][j] +
                             board[i][j - 1] +
                             board[i - 1][j]
                     ) / 4;
            difference = fabs(result - board[i][j]);
            if (difference >= precision) {
                printf("Invalid value found at (%d, %d).\n", j, i);
                printf("Precision: %f\n", precision);
                printf("Average of neighbours: %f\n", result);
                printf("Value of centre cell: %f\n", board[i][j]);
                printf("| %f - %f | = %.10f\n", board[i][j], result, difference);
                printf("       %.4f\n\n", board[i - 1][j]);
                printf("%.4f %.4f %.4f\n\n", board[i][j - 1], board[i][j], board[i][j + 1]);
                printf("       %.4f\n\n", board[i + 1][j]);
                return -1;
            }
        }
    }
//    printf("Answer correct.\n");
    return 0;
}

int cleanup() {
    for (int i = 0; i < boardSize; ++i) {
        free(board[i]);
        free(results[i]);
    }
    free(board);
    free(results);
    free(rowStart);
    free(rowEnd);
    free(threadDone);
    pthread_barrier_destroy(&barrier);
    return 0;
}

int init() {
//    Allocate 2D array of doubles for board
    board = malloc(boardSize * sizeof(*board));
    results = malloc(boardSize * sizeof(*results));
    for (int i = 0; i < boardSize; ++i) {
        board[i] = malloc(boardSize * sizeof(*board[i]));
        results[i] = malloc(boardSize * sizeof(*results[i]));
    }

//    Build initial values for the board (and results)
    for (int i = 0; i < boardSize; ++i) {
        board[0][i] = 1.;
        results[0][i] = 1.;
    }

    for (int i = 1; i < boardSize; ++i) {
        board[i][0] = 1.;
        results[i][0] = 1.;
        for (int j = 1; j < boardSize; ++j) {
            board[i][j] = 0.;
            results[i][j] = 0.;
        }
    }


//    Initialise barrier
    int errorCode = pthread_barrier_init(&barrier, NULL, nThreads + 1);
    if (errorCode != 0) {
        printf("Error making barrier");
        cleanup();
        exit(errorCode);
    }


    rowStart = malloc(nThreads * sizeof(int));
    rowEnd = malloc(nThreads * sizeof(int));

//    Allocate flags and mark threads as not done.
    threadDone = malloc(nThreads * sizeof(int));
    for (int i = 0; i < nThreads; ++i) {
        threadDone[i] = 0;
    }

    run = 1;
    return 0;
}

// Swaps the pointers of 2 boards.
void swap(double ***a, double ***b) {
    double **tmp = *a;
    *a = *b;
    *b = tmp;
}

int main(int argc, char *argv[]) {
//    size precision threads
    if (argc != 4) {
        printf("Incorrect number of arguments: 3 required.");
        exit(1);
    }

    boardSize = (int) strtol(argv[1], NULL, 10);
    precision = pow(10, (int) strtol(argv[2], NULL, 10));
    nThreads = (int) strtol(argv[3], NULL, 10);

    struct timeval tval_start, tval_seqEnd, tval_paraEnd, tval_seqDuraton, tval_paraDuration;

//    printf("Relaxing %dx%d grid to precision %f with %d threads...\n", boardSize, boardSize, precision, nThreads);
    gettimeofday(&tval_start, NULL);
    init();

//    Allocate thread numbers
    const int defaultHeight = (boardSize - 2) / nThreads;
    int remainder = (boardSize - 2) % nThreads;

    int start, end;
    int currentHeight = 0;


    for (int i = 0; i < nThreads; ++i) {
        start = currentHeight + 1;
        end = start + defaultHeight;

        if (remainder > 0) {
            --remainder;
            ++end;
        }

        currentHeight += (end - start);

        rowEnd[i] = end;
        rowStart[i] = start;
    }


    pthread_t threads[nThreads];
    int thread_ids[nThreads];
    int errorCode;
    for (int i = 0; i < nThreads; ++i) {
        thread_ids[i] = i;
        errorCode = pthread_create(&threads[i], NULL, threadFunction, &thread_ids[i]);
        if (errorCode != 0) {
            printf("Error %d when creating worker thread #%d.", errorCode, i);
            cleanup();
            exit(errorCode);
        }
    }

    gettimeofday(&tval_seqEnd, NULL);

    while (run == 1) {
//        printf("MAIN: Waiting for first barrier...\n");
        pthread_barrier_wait(&barrier);
//        printf("MAIN: Checking results...\n");

        swap(&board, &results);

        run = 0;
        for (int i = 0; i < nThreads; ++i) {
            if (threadDone[i] == 0) {
                run = 1;
                break;
            }
        }
        pthread_barrier_wait(&barrier);
    }

//    printf("MAIN: Done.\n");

//    Wait for all threads to finish
    for (int i = 0; i < nThreads; ++i) {
        pthread_join(threads[i], NULL);
    }

    gettimeofday(&tval_paraEnd, NULL);
    timersub(&tval_seqEnd, &tval_start, &tval_seqDuraton);
    timersub(&tval_paraEnd, &tval_seqEnd, &tval_paraDuration);


    printBoard();

//    printf("\nDuration: %lds\n", duration);
//    printf("%d,%d,%ld\n", boardSize, nThreads, duration);
    printf("Seq: %ld.%06ld, Para: %ld.%06ld, size: %d, nThreads: %d\n", (long int)tval_seqDuraton.tv_sec, (long int)tval_seqDuraton.tv_usec, (long int)tval_paraDuration.tv_sec, (long int)tval_paraDuration.tv_usec, boardSize, nThreads);
    int result = checkBoard();
    cleanup();
    if (result == -1) {
        printf("Error - Board not correct.");

        exit(-1);
    }

    return result;
}
