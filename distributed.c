#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <mpi.h>

// Tags for messages between processors
const int TAG_SEND_END = 950;
const int TAG_SEND_START = 951;

// Rank of Root processor
const int ROOT = 0;

//    How many iterations to perform before checking?
const int CHECK_FREQ = 1;

int printBoard(double **board, int size) {
    for (int i = 0; i < size; ++i) {
        for (int j = 0; j < size; ++j) {
            printf("%.4f ", board[i][j]);
        }
        printf("\n\n");
    }
    return 0;
}

int checkBoard(double **board, int boardSize, double precision) {
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
    printf("Answer correct.\n");
    return 0;
}

void printBuff(double *buff, int size) {
    for (int i = 0; i < size; ++i) {
        printf("%.2f ", buff[i]);
    }
    printf("\n");
}

void printIntArr(int *buff, int size) {
    for (int i = 0; i < size; ++i) {
        printf("%d ", buff[i]);
    }
    printf("\n");
}

int allValuesEqual(const int *arr, int value, int size) {
//    Return 0 if all values are equal. -1 otherwise
    for (int i = 0; i < size; ++i) {
        if (arr[i] != value) {
            return -1;
        }
    }
    return 0;
}

int main(int argc, char **argv) {
    struct timeval tval_start, tval_seqEnd, tval_seqDuraton;
    double **board;
    double **results;

//    Stores the start/end indexes for each worker
    int *rowStart;
    int *rowEnd;


    int errNum, nProcs, myId;


//    Get number of instances, and my id
    errNum = MPI_Init(&argc, &argv);

    if (errNum != 0) {
        printf("Error while calling MPI_Init()!");
        exit(errNum);
    }


    errNum = MPI_Comm_rank(MPI_COMM_WORLD, &myId);

    if (errNum != 0) {
        printf("Error while calling MPI_Comm_rank()!");
        exit(errNum);
    }

    errNum = MPI_Comm_size(MPI_COMM_WORLD, &nProcs);

    if (errNum != 0) {
        printf("Error while calling MPI_Comm_size()!");
        exit(errNum);
    }

    //    size precision threads
    if (argc != 3) {
        printf("Incorrect number of arguments: 3 required, %d provided.", argc);
        exit(1);
    }

    int boardSize = (int) strtol(argv[1], NULL, 10);
    const int power = (int) strtol(argv[2], NULL, 10);
    double precision = pow(10, power);

    if (myId == ROOT) {
        gettimeofday(&tval_start, NULL);
    }

    board = malloc(boardSize * sizeof(*board));
    results = malloc(boardSize * sizeof(*results));
    for (int i = 0; i < boardSize; ++i) {
        board[i] = malloc(boardSize * sizeof(*board[i]));
        results[i] = malloc(boardSize * sizeof(*results[i]));
    }

    rowStart = malloc(nProcs * sizeof(int));
    rowEnd = malloc(nProcs * sizeof(int));

//    printf("Building board\n");

    for (int i = 0; i < boardSize; ++i) {
        board[0][i] = 1;
        results[0][i] = 1;
    }

    for (int i = 1; i < boardSize; ++i) {
        board[i][0] = 1;
        results[i][0] = 1;

        for (int j = 1; j < boardSize; ++j) {
            board[i][j] = 0;
            results[i][j] = 0;
        }
    }

    //    Allocate work for each thread.
    const int defaultHeight = (boardSize - 2) / nProcs;

//    Number of extra rows that will need allocating
    int remainder = (boardSize - 2) % nProcs;

    int start, end;
    int currentHeight = 0;

    for (int i = 0; i < nProcs; ++i) {
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
//    printf("%d: %d-%d\n", myId, rowStart[myId], rowEnd[myId]);

//    For gathering results
    const int sendBuffSize = (rowEnd[0] - rowStart[0]) *
                             boardSize;
    double *sendBuff = malloc(sendBuffSize * sizeof(double));

    const int recvBuffSize = sendBuffSize * nProcs;
    double *recvBuff = malloc(recvBuffSize * sizeof(double));

//    Find the indexes this proc should work between
    const int my_start = rowStart[myId];
    const int my_end = rowEnd[myId];

//    Number of cells that this proc will work on
    const int target = (my_end - my_start) * (boardSize - 2);

//    Keeps track of how many cells are correct
    int nCorrect = 0;

//    Number of iterations
    int count = 0;

    int done;
    int *doneResults = malloc(nProcs * sizeof(int));

    double currentValue, newValue, diff;
    while (1) {
//        Reset for this iteration
        done = 0;
        nCorrect = 0;

//        Averaging
        for (int i = my_start; i < my_end; ++i) {
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

//        Check if we're done yet
//        printf("%d: Correct=%d; target=%d\n", myId, nCorrect, target);
        done = (nCorrect == target);

//        Swap pointers of result and board so we use new values next iteration
        double **tmp = board;
        board = results;
        results = tmp;

        if ((count % CHECK_FREQ) == 0) {
//            printf("%d: Gathering for check...\n", myId);
            MPI_Gather(&done, 1, MPI_INT, doneResults, 1,
                       MPI_INT, ROOT, MPI_COMM_WORLD);
            int stop;
            if (myId == ROOT) {
//                printf("%d: Checking...\n", myId);
                int result = allValuesEqual(doneResults, 1, nProcs);
                if (result == 0) {
//                    printIntArr(doneResults, nProcs);
                    stop = 1;
                } else {
                    stop = 0;
                }
//                printf("%d: Broadcast result (%d)...\n", myId, stop);
            }
            MPI_Bcast(&stop, 1, MPI_INT, ROOT, MPI_COMM_WORLD);
            if (stop == 1) {
//                All processors must have buffers of same length, so we will
//                take the maximum workload length (which the first
//                processors will always have), and any procs with less
//                workload can fill the buffer with NAN to indicate it's empty.


                int col = 0;
                int row = my_start;
                for (int i = 0; i < sendBuffSize; ++i) {
                    if (row == my_end) {
                        sendBuff[i] = NAN;
                        continue;
                    }

                    sendBuff[i] = board[row][col];

                    ++col;
//                    Check if we need to wrap to next line
                    if (col == boardSize) {
                        col = 0;
                        ++row;
                    }
                }

                MPI_Gather(sendBuff, sendBuffSize, MPI_DOUBLE,
                           recvBuff, sendBuffSize, MPI_DOUBLE,
                           ROOT, MPI_COMM_WORLD);
                free(sendBuff);

                if (myId == ROOT) {

                    col = 0;
                    row = rowStart[0];

                    for (int i = 0; i < recvBuffSize; ++i) {
//                        Ignore empty values. Skip forward until we reach a
//                        real value.
                        if (isnan(recvBuff[i])) {
                            continue;
                        }

                        board[row][col] = recvBuff[i];
                        ++col;

//                        Check if we need to wrap to next line
                        if (col == boardSize) {
                            col = 0;
                            ++row;
                        }
                    }

//                    printBoard(board, boardSize);
                    checkBoard(board, boardSize, precision);
                }
                free(recvBuff);
                break;
            }

        }

//      Processing is not done, swap rows with neighbours
//
//      Threads that are in the middle of the board will run all of the
//      following code.
//      The top thread will only run the first if-block, and the bottom thread
//      will only run the second if-block

        MPI_Request REQUEST_RECV_START;

        if (myId != nProcs - 1) {
//            Send my bottom edge, don't wait
            MPI_Send(board[my_end - 1], boardSize, MPI_DOUBLE,
                     myId + 1, TAG_SEND_END, MPI_COMM_WORLD);

//            Receive my neighbor's top edge, don't wait.
//            We will wait after the following if-block. Prevents deadlock.
            MPI_Irecv(board[my_end], boardSize, MPI_DOUBLE, myId + 1,
                      TAG_SEND_START, MPI_COMM_WORLD, &REQUEST_RECV_START);
        }

        if (myId != 0) {
//            Send my top edge, don't wait
            MPI_Send(board[my_start], boardSize, MPI_DOUBLE, myId - 1,
                     TAG_SEND_START, MPI_COMM_WORLD);

//            Receive and wait for my neighbor's bottom edge
            MPI_Recv(board[my_start - 1], boardSize, MPI_DOUBLE,
                     myId - 1, TAG_SEND_END, MPI_COMM_WORLD,
                     NULL);
        }


        if (myId != nProcs - 1) {
//            Wait for neighbor's top edge
            MPI_Wait(&REQUEST_RECV_START, NULL);
        }

        ++count;
    }
    if (myId == ROOT) {
        gettimeofday(&tval_seqEnd, NULL);
        timersub(&tval_seqEnd, &tval_start, &tval_seqDuraton);
        printf("%d,%d,%ld.%06ld\n", boardSize, nProcs,
               (long int) tval_seqDuraton
                       .tv_sec,
               (long int) tval_seqDuraton.tv_usec);
    }

//    End
    errNum = MPI_Finalize();
    if (errNum != 0) {
        printf("Error while calling MPI_Finalize()!");
        exit(errNum);
    }
}
