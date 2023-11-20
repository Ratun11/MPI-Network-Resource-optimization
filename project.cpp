#include <iostream>
#include <iomanip>
#include <sstream>
#include <random>
#include <exception>
#include <vector>
#include <mpi.h>
#include <algorithm>

using namespace std;

template <class T>
T enter_next_param(const string prompt) {
    T num;
    while (1) {
        string input = "";
        cout << prompt << endl;
        getline(cin, input);
        stringstream myStream(input);
        if (myStream >> num) break;
        cout << "Invalid Input" << endl << endl;
    }
    return num;
}

void get_parameters(int *list_size, int *num_bins, double *max_bin_val, double *min_bin_val, bool *num_list_out, int argc, char *argv[], int rank, int numprocs) {
    if (argc == 1) {
        if (rank == 0) {
            *list_size = enter_next_param<int>("Enter number of numbers in list:");
            *num_bins = enter_next_param<int>("Enter number of bins in histogram:");
            *max_bin_val = enter_next_param<double>("Enter Upper Range of histogram:");
            *min_bin_val = enter_next_param<double>("Enter Lower Range of histogram:");
            char print_out = enter_next_param<char>("Suppress numbers array printing (y or n):");
            if (print_out == 'y' || print_out == 'Y') *num_list_out = false;
        }
        MPI_Bcast(list_size, 1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(num_bins, 1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(max_bin_val, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Bcast(min_bin_val, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    }
    else if (argc == 5 || argc == 6) {
        *list_size = atoi(argv[1]);
        *num_bins = atoi(argv[2]);
        *max_bin_val = atof(argv[3]);
        *min_bin_val = atof(argv[4]);
        if (argc == 6) *num_list_out = false;
    }
    else {
        if (rank == 0) {
            cout << "Usage: mpirun -np [num MPI processes] histogram_MPI [list size num_bins max_bin_val min_bin_val]" << endl;
        }
        MPI_Finalize();
        exit(0);
    }
}

#define MEAN  50.0f
#define STDEV 20.0f

#define LOWER_BOUND 0.1f
#define	UPPER_BOUND 1000.0f

void create_list(vector<double>& numbers, int list_size) {
    default_random_engine generator(123546);
    uniform_real_distribution<double> distribution(LOWER_BOUND, UPPER_BOUND);

    for (int i = 0; i < list_size; i++) {
        numbers[i] = distribution(generator);
    }
}
void create_list2(vector<double>& numbers, int list_size){
  float min = 10.0, max = 50.0;
  std::random_device rd;
  std::mt19937 eng(rd());
  std::uniform_real_distribution<> distr(min,max);

  std::generate_n(std::back_inserter(numbers),list_size,[&]() {return distr(eng);});
  
}

void create_list3(vector<double>& numbers, int list_size){
    for(int i=0;i<list_size;i++){
      numbers.push_back(rand()%(100-0+1)+0);
 }
}

void create_list4(double* numbers, int list_size) {
    default_random_engine generator(123546);
    uniform_real_distribution<double> distribution(LOWER_BOUND, UPPER_BOUND);

    for (int i = 0; i < list_size; i++) {
        numbers[i] = distribution(generator);
    }
}

void print_matrix2(double* numbers,int list_size) {
    for (int i = 0;i<list_size;i++) {
        cout << numbers[i] << endl;
    }
}
void print_matrix(const vector<double>& numbers) {
    for (const double& num : numbers) {
        cout << num << endl;
    }
}

//void scatter(vector<double>& numbers,vector<double>& group,int num_size, int root, int rank,int numprocs)
/*void scatter(double *numbers,double *group,int num_size, int root, int rank,int numprocs)
{
   MPI_Status status;
   int tag = 234;
   
   // determine number of elements in subarray groups to be processed by
   // each MPI process assuming a perfectly even distribution of elements 
   const int number_elements_per_section = num_size/numprocs;

   // if root MPI process send portion of numbers array to each of the
   // the other MPI processes as well as make a copy of the portion
   // of the numbers array that is slated for the root MPI process
   if (rank==root) {
      int begin_element=0;

      for(int mpiproc=0;mpiproc<numprocs;mpiproc++) {
         
         int rem = num_size%numprocs;
	      float quo = static_cast<float>(num_size)/numprocs;
      	int number_elements_per_section = (mpiproc < rem) ? ceil(quo) : floor(quo);
   
         // in MPI root process case just copy the appropriate subsection
         // locally from the numbers array over to the group array
         if (mpiproc==root) {
            for (int i=0;i<number_elements_per_section;i++) 
               group[i]=numbers[i+begin_element];
         }
         // if not the root process send the subsection data to
         // the next MPI process
         else {
            MPI_Send(&numbers[begin_element],number_elements_per_section,
               MPI_DOUBLE,mpiproc,tag,MPI_COMM_WORLD);
         }
         // point to next unsent or uncopied data in numbers array
         begin_element += number_elements_per_section;
      }
   }
   // if a non root process just receive the data
   else {
      int rem = num_size%numprocs;
	   float quo = static_cast<float>(num_size)/numprocs;
	   int number_elements_per_section = (rank < rem) ? ceil(quo) : floor(quo);
	
      MPI_Recv(group,number_elements_per_section,MPI_DOUBLE,
               root,MPI_ANY_TAG,MPI_COMM_WORLD,&status);  
   }
}*/

void scatter(double *numbers,double *group,int num_size, int root, int rank,int numprocs)
{
   MPI_Status status;
   int tag = 234;

   const int number_elements_per_section = num_size/numprocs;

   if (rank==root) {
      int begin_element=0;

      for(int mpiproc=0;mpiproc<numprocs;mpiproc++) {

         if (mpiproc==root) {
            for (int i=0;i<number_elements_per_section;i++)
               group[i]=numbers[i+begin_element];
         }
         else {
            MPI_Send(&numbers[begin_element],number_elements_per_section,
               MPI_DOUBLE,mpiproc,tag,MPI_COMM_WORLD);
         }
         // point to next unsent or uncopied data in numbers array
         begin_element += number_elements_per_section;
      }
   }
   // if a non root process just receive the data
   else {
      MPI_Recv(group,number_elements_per_section,MPI_DOUBLE,
               root,MPI_ANY_TAG,MPI_COMM_WORLD,&status);
   }
}

int main(int argc, char *argv[]) {
    int numprocs, rank;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int list_size, num_bins;
    double max_val, min_val;
    bool print_flg = true;

    get_parameters(&list_size, &num_bins, &max_val, &min_val, &print_flg, argc, argv, rank, numprocs);
    const double bin_range = (max_val - min_val) / num_bins;

    vector<double> numbers(list_size);   
    vector<int> histogram;

     const char* band[13]={"mm","W_band","V_band","Ka_band","K_band","Ku_band","X_band","C_band","S_band","L_band","UHF","VHF","HF"};

    if (rank == 0) {
       numbers.resize(list_size);
        histogram.resize(num_bins);

        create_list(numbers, list_size);

        if (print_flg) {
            cout << "numbers matrix =" << endl;
            print_matrix(numbers);
            cout << endl;
        }
    }

    const int number_elements_per_section = list_size / numprocs;
    const int num_leftover_elements = list_size % numprocs;

    vector<int> sendcounts(numprocs);
    vector<int> displs(numprocs);

    for (int i = 0; i < numprocs; i++) {
        sendcounts[i] = number_elements_per_section + ((i < num_leftover_elements) ? 1 : 0);
        displs[i] = i * number_elements_per_section + min(i, num_leftover_elements);
    }

    const int numbers_local_array_sz = sendcounts[rank];
    vector<double> numbers_local(numbers_local_array_sz);
   vector<int> histogram_local(num_bins);

    MPI_Scatterv(numbers.data(), sendcounts.data(), displs.data(), MPI_DOUBLE, numbers_local.data(), numbers_local_array_sz, MPI_DOUBLE, 0, MPI_COMM_WORLD);
   // scatter(numbers,numbers_local,list_size,0,rank,numprocs);

    double t1 = MPI_Wtime();

    for (int i = 0; i < num_bins; i++) {
        histogram_local[i] = 0;
    }

    for (int i = 0; i < numbers_local_array_sz; i++) {
        double num = numbers_local[i];
        if (num >= min_val && num < max_val) {
            num -= min_val;
            int bin_index = static_cast<int>((num / bin_range));
            if (bin_index >= 0 && bin_index < num_bins) {
                histogram_local[bin_index]++;
            }
        }
    }
   double t2 = MPI_Wtime();

    MPI_Reduce(histogram_local.data(), histogram.data(), num_bins, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
   

    if (rank == 0) {
        cout << "Histogram Information" << endl;
        for (int i = 0; i < num_bins; i++) {
            double bin_center = min_val + bin_range * (i + 0.5);
            cout << fixed << setprecision(8) << band[i] << "  " << histogram[i] << endl;
       
 }
        cout << numprocs << " MPI Process OpenMPI Implementation" << endl << flush;
       cout << "Execution Time = "<< t2-t1 << "seconds"<< endl << flush;
    }


    if (rank == 0) {
        numbers.clear();
        histogram.clear();
    }

    MPI_Finalize();
}
