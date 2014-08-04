#include <mpi.h>
#include <vector>
#define DIETAG     10000000
#define MASTER_RANK 0
#include <limits.h>

using namespace std;

const int npar = 3;
const int nmet = npar + 1;

void master();
void slave();

int main(int argc, char* argv[]) {
    int myrank;
    MPI_Init(&argc, &argv);   /* initialize MPI */
    MPI_Comm_rank(
      MPI_COMM_WORLD,         /* always use this */
      &myrank);               /* process rank, 0 thru N-1 */
    if (myrank == 0) {
        // Example of how to get upper bound for tags
        void* v;
        int flag;
        int vval;
        MPI_Comm_get_attr( MPI_COMM_WORLD, MPI_TAG_UB, &v, &flag );
        if (!flag) {
            fprintf( stderr, "Could not get TAG_UB\n" );fflush(stderr);
        } else {
            vval = *(int*)v;
            if (vval < 32767) {
                fprintf( stderr, "Got too-small value (%d) for TAG_UB\n", vval ); fflush(stderr);
            } else {
                fprintf( stderr, "MAX_INT: %d\n", INT_MAX); fflush(stderr);
                fprintf( stderr, "Got %d\n", vval ); fflush(stderr);
            }
        }

        master();
    } else {
        slave();
    }
    MPI_Finalize();           /* cleanup MPI */
    return 0;
}

void master() {
    int ntasks = 10;
    int scheduledTasks=0;
    int mpi_size, rank;
    vector<double> result(ntasks*nmet);

    double rec_buffer[nmet];
    MPI_Status status;
    MPI_Comm_size(
      MPI_COMM_WORLD,   /* always use this */
      &mpi_size);       /* #processes in application */
    
    double work[ntasks*npar];
    for (int i = 0; i<ntasks*npar; ++i) work[i] = i;

    // Seed the slaves.
    for (rank = 1; rank < mpi_size; ++rank) {
        MPI_Send(&work[scheduledTasks*npar], /* message buffer */
                 npar,                       /* number of elements */
                 MPI_DOUBLE,                 /* data item is a double */
                 rank,                       /* destination process rank */
                 scheduledTasks,             /* user chosen message tag */
                 MPI_COMM_WORLD);            /* always use this */
        scheduledTasks++;
    }

    
    // Receive a result from any slave and dispatch a new work
    // request work requests have been exhausted.
    
    while ( scheduledTasks < ntasks ) { 
        MPI_Recv(&rec_buffer,    /* message buffer */
                 nmet,           /* nmet data items */
                 MPI_DOUBLE,     /* of type double */
                 MPI_ANY_SOURCE, /* receive from any sender */
                 MPI_ANY_TAG,    /* any type of message */
                 MPI_COMM_WORLD, /* always use this */
                 &status);       /* received message info */
        
        for (int m = 0; m< nmet; m++) result[status.MPI_TAG*nmet + m] = rec_buffer[m];

        MPI_Send(&work[scheduledTasks*npar], 
                 npar, 
                 MPI_DOUBLE, 
                 status.MPI_SOURCE,
                 scheduledTasks, 
                 MPI_COMM_WORLD);
        scheduledTasks++;
    }
    // Receive results for outstanding work requests.
    
    cerr << "almost done1\n";
    for (rank = 1; rank < mpi_size; ++rank) {
        MPI_Recv(&rec_buffer, 
                 nmet, 
                 MPI_DOUBLE, 
                 MPI_ANY_SOURCE,
                 MPI_ANY_TAG, 
                 MPI_COMM_WORLD, 
                 &status);
        for (int m = 0; m< nmet; m++) result[status.MPI_TAG*nmet + m] = rec_buffer[m];
    }
     
    // Tell all the slaves to exit.
    for (rank = 1; rank < mpi_size; ++rank) {
        MPI_Send(0, 0, MPI_INT, rank, DIETAG, MPI_COMM_WORLD);
    }

    for (int r = 0; r<ntasks; r++) {
        for (int c = 0; c<nmet; c++) {
            cout << result[r*nmet + c] << '\t';
        }
        cout << endl;
    }
}

void slave() {
    double work[npar];
    double result[nmet];
    MPI_Status status;
    while (1) {
        MPI_Recv(&work, 
                 npar, 
                 MPI_DOUBLE, 
                 MASTER_RANK, 
                 MPI_ANY_TAG,
                 MPI_COMM_WORLD, 
                 &status);
        
        // make sure to reset the buffers as needed
        for (int i = 0; i<nmet; i++) result[i] = 0.0;
        
        // Check the tag of the received message.
        if (status.MPI_TAG == DIETAG) {
            return;
        }

        for (int m = 0; m< npar; m++) {
            result[m] = work[m]*work[m];
            result[npar] += work[m];
        }

        MPI_Send(&result, 
                 nmet, 
                 MPI_DOUBLE, 
                 MASTER_RANK, 
                 status.MPI_TAG, 
                 MPI_COMM_WORLD);
    }
}
