#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <sys/time.h>


// Parameters
int boardSize;
double precision;

// Matrices
double **board, **results;


int printBoard() {
    for (int i = 0; i < boardSize; ++i) {
        for (int j = 0; j < boardSize; ++j) {
            printf("%.4f ", board[i][j]);
        }
        printf("\n\n");
    }
    return 0;
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
                printf("| %f - %f | = %.10f\n", board[i][j], result,
                       difference);
                printf("       %.4f\n\n", board[i - 1][j]);
                printf("%.4f %.4f %.4f\n\n", board[i][j - 1], board[i][j],
                       board[i][j + 1]);
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
    return 1;
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
    if (argc != 3) {
        printf("Incorrect number of arguments: 3 required, %d provided.", argc);
        exit(1);
    }

    boardSize = (int) strtol(argv[1], NULL, 10);
    precision = pow(10, (int) strtol(argv[2], NULL, 10));

    struct timeval tval_start, tval_seqEnd, tval_seqDuraton;

//    printf("Relaxing %dx%d grid to precision %f with 1 thread...\n", boardSize, boardSize, precision, nThreads);
    gettimeofday(&tval_start, NULL);
    init();
    double newValue, currentValue, diff;
    int run = 1;
    while (run == 1) {
        run = 0;
        swap(&board, &results);
        for (int i = 1; i < boardSize - 1; ++i) {
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
                if (diff >= precision) { run = 1; }
            }
        }
    };

//    printf("MAIN: Done.\n");

//    Wait for all threads to finish

    gettimeofday(&tval_seqEnd, NULL);
    timersub(&tval_seqEnd, &tval_start, &tval_seqDuraton);

//    printBoard();

    printf("%d,%d,%ld.%06ld\n", boardSize, 1, (long int) tval_seqDuraton.tv_sec,
           (long int) tval_seqDuraton.tv_usec);

    int result = checkBoard();

    if (result == -1) {
        printf("Error - Board not correct.");
        exit(-1);
    }

    cleanup();

    return result;
}
