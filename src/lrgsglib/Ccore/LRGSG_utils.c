#include "LRGSG_utils.h"
//
/**
 * @brief Prints the current working directory to the specified output stream.
 *
 * This function retrieves the current working directory and prints it to the
 * specified output stream. If the output stream is invalid or an error occurs
 * while retrieving the directory, an appropriate error message is printed.
 *
 * @param output_stream (FILE*) The output stream to print to (e.g., stdout, stderr).
 */
void print_cwd(FILE *output_stream) {
    if (output_stream == NULL) {
        fprintf(stderr, "Invalid output stream.\n");
        return;
    }
    char cwd[1024];
    if (getcwd(cwd, sizeof(cwd)) == NULL) {
        fprintf(stderr, "Error getting current working directory: %s\n", strerror(errno));
    } else {
        fprintf(output_stream, "Current working directory: %s\n", cwd);
    }
}
/**
 * @brief Fills a section of an array with a specified value.
 *
 * This function sets all elements in a specified range of an array to a given value.
 *
 * @param arr   (double*) Pointer to the array to be filled.
 * @param start (size_t) Starting index of the section to fill (inclusive).
 * @param end   (size_t) Ending index of the section to fill (exclusive).
 * @param value (double) The value to set in the specified range.
 *
 * @note The function does not perform bounds checking. Ensure `start < end` 
 *       and `end` is within the bounds of the array to avoid undefined behavior.
 */
void fill_array_with_value(double *arr, size_t start, size_t end, double value) {
    for (size_t i = start; i < end; i++)
        arr[i] = value;
}
/**
 * @brief Flips the spin at the given index in the spin array.
 *
 * This function negates the value of the spin at the specified index in the array.
 *
 * @param nd  (size_t) Index of the spin to flip.
 * @param s   (spin_tp) Array of spins.
 */
void flip_spin(size_t nd, spin_tp s) {
    *(s + nd) = -*(s + nd);
}
/**
 * @brief Calculates the total magnetization of the spin system.
 * 
 * @param N   (size_t) Number of spins in the system.
 * @param s   (spin_tp) Pointer to the array of spins.
 * 
 * @return    (double) The total magnetization of the system.
 */
double calc_ext_magn(size_t N, spin_tp s) {
    double M = 0.;
    for (size_t i = 0; i < N; i++)
        M += *(s + i);
    return M;
}
/**
 * @brief Calculates the average magnetization of the spin system.
 * 
 * @param N   (size_t) Number of spins in the system.
 * @param s   (spin_tp) Pointer to the array of spins.
 * 
 * @return    (double) The average magnetization of the system.
 */
double calc_magn(size_t N, spin_tp s) {
    return calc_ext_magn(N, s) / N;
}
/**
 * @brief Calculates the squared magnetization of the spin system.
 *
 * This function computes the sum of the squares of the spins in the system.
 *
 * @param N (size_t) Number of spins in the system.
 * @param s (spin_tp) Pointer to the array of spins.
 *
 * @return (double) The squared magnetization of the system.
 */
double calc_ext_magn2(size_t N, spin_tp s) {
    double m2 = 0.;
    for (size_t i = 0; i < N; i++)
        m2 += *(s + i) * *(s + i);
    return m2;
}
/**
 * @brief Calculates the cluster magnetization of a subset of spins.
 *
 * This function computes the average magnetization of a cluster of spins
 * specified by their indices.
 *
 * @param cli_l (size_t) Length of the cluster (number of spins in the cluster).
 * @param cli (size_tp) Array of indices representing the cluster.
 * @param s (spin_tp) Pointer to the array of spins.
 *
 * @return (double) The average magnetization of the cluster.
 */
double calc_clust_magn(size_t cli_l, size_tp cli, spin_tp s) {
    double clm = 0.;
    for (size_t i = 0; i < cli_l; i++)
        clm += s[cli[i]];
    return clm / cli_l;
}
/**
 * @brief Calculates the weighted magnetization of a node's neighbors.
 * 
 * This function computes the sum of the weighted magnetization of a node's neighbors
 * in a spin system. The weights are determined by the edge weights between the node
 * and its neighbors.
 * 
 * @param nd (size_t) The index of the current node.
 * @param n_nn (size_t) The number of neighbors of the current node.
 * @param s (spin_tp) The array representing the spin states of all nodes.
 * @param neighs (size_tp *) A 2D array where each row contains the indices of the neighbors for each node.
 * @param edgl (double_p *) A 2D array where each row contains the edge weights for the neighbors of each node.
 * 
 * @return (double) The weighted magnetization of the neighbors of the node.
 */
double neigh_weight_magn(size_t nd, size_t n_nn, spin_tp s, size_tp *neighs,
                         double_p *edgl) {
    double sum = 0.;
    for (size_t i = 0; i < n_nn; i++)
        sum += *(*(edgl + nd) + i) * *(s + *(*(neighs + nd) + i));
    return sum;
}
/**
 * @brief Calculates the weighted magnetization of a node's neighbors.
 *
 * This function computes the sum of the weighted spins of a node's neighbors
 * based on the provided neighbor data structure.
 *
 * @param ne (NodeEdges) The structure containing the neighbors and their weights.
 * @param s (spin_tp) Pointer to the array of spins.
 * @param n_nn (size_t) The number of neighbors of the node.
 *
 * @return (double) The weighted magnetization of the node's neighbors.
 */
double neighWeight_magn(NodeEdges ne, spin_tp s, size_t n_nn) {
    double sum = 0.;
    for (size_t i = 0; i < n_nn; i++)
        sum += *(ne.weights + i) * *(s + *(ne.neighbors + i));
    return sum;
}
/**
 * @brief Calculates the total energy of the spin system.
 *
 * This function computes the total energy of the spin system based on the
 * interactions between spins and their neighbors.
 *
 * @param N (size_t) The number of spins in the system.
 * @param s (spin_tp) Pointer to the array of spins.
 * @param nlen (size_tp) Array containing the number of neighbors for each spin.
 * @param ne (NodesEdges) Array of neighbor data structures.
 *
 * @return (double) The total energy of the spin system.
 */
double calc_totEnergy(size_t N, spin_tp s, size_tp nlen, NodesEdges ne) {
    double sum = 0.;
    for (size_t i = 0; i < N; i++) {
        double tmp = *(s + i) * neighWeight_magn(ne[i], s, nlen[i]);
        sum += tmp;
    }
    return - sum / (2 * N);
}
/**
 * @brief Calculates the full energy of the spin system.
 *
 * This function computes the total energy of the spin system by considering
 * all interactions between spins and their neighbors, including edge weights.
 *
 * @param N (size_t) The number of spins in the system.
 * @param s (spin_tp) Pointer to the array of spins.
 * @param nlen (size_tp) Array containing the number of neighbors for each spin.
 * @param neighs (size_tp*) Pointer to the array of neighbor indices for each spin.
 * @param edgl (double_p*) Pointer to the array of edge weights for each spin.
 *
 * @return (double) The total energy of the spin system.
 */
double calc_energy_full(size_t N, spin_tp s, size_tp nlen, size_tp *neighs, double_p *edgl) {
    double sum = 0., tmp = 0.;
    for (size_t i = 0; i < N; i++) {
        tmp = *(s + i) * neigh_weight_magn(i, *(nlen + i), s, neighs, edgl);
        sum += tmp;
    }
    return -sum / N;
}
/**
 * @brief Checks if flipping any spin in the system would increase the energy (stability at T = 0).
 * 
 * @param N      (size_t) Number of spins in the system.
 * @param s      (spin_tp) Array of spins representing the system.
 * @param nlen   (size_t*) Array containing the number of neighbors for each spin.
 * @param ne     (NodesEdges) Array of neighbor data structures.
 * 
 * @return       (bool) True if the system is stable, false otherwise.
 */
bool glauber_isStableAtZeroTemp(size_t N, spin_tp s, size_tp nlen, NodesEdges ne) {
    for (size_t i = 0; i < N; i++) {
        double DE = 2 * (*(s + i)) * neighWeight_magn(ne[i], s, nlen[i]);
        if (DE < 0) {
            return false;
        }
    }
    return true;
}
/**
 * @brief Generate a logarithmically-spaced vector
 *
 * This function creates a vector of doubles where each value is logarithmically
 * spaced between the start and stop values.
 *
 * @param start The base-10 logarithm of the start value of the desired range
 * @param stop The base-10 logarithm of the stop value of the desired range
 * @param num The number of values to be generated in the vector
 * @return A dynamically allocated array of doubles representing the logarithmically-spaced vector
 *
 * @note The returned array must be freed by the caller using free() when no longer needed.
 */
extern double* logspace(double start, double stop, int num) {
    double* vec = (double*)malloc(num * sizeof(double));
    double step = (stop - start) / (num - 1);

    for (int i = 0; i < num; i++) {
        vec[i] = pow(10, start + i * step);
    }

    return vec;
}
/**
 * @brief Generate a logarithmically-spaced vector of integers.
 *
 * This function creates a vector of integers where each value is logarithmically
 * spaced between the start and stop values.
 *
 * @param stop The base-10 logarithm of the stop value of the desired range
 * @param num Pointer to the number of values to be generated in the vector. The value may be adjusted within the function.
 * @return A dynamically allocated array of integers representing the logarithmically-spaced vector
 *
 * @note The returned array must be freed by the caller using free() when no longer needed.
 */
extern int* logspace_int(double stop, int* num) {
    int* vec = (int*)malloc((*num) * sizeof(int));
    double step = stop / (*num - 1);

    for (int i = 0; i < *num; i++) {
        vec[i] = (int)round(pow(10, i * step));
    }
    // Check if the first two points are the same
    if (vec[0] == vec[1]) {
        int num_max = 1 + (int)(1 / (log10(2) / stop));
        num_max = num_max < *num ? num_max : *num;
        // Free the previously allocated memory
        free(vec);
        // Update the num value
        *num = num_max;
        vec = (int*)malloc((*num) * sizeof(int));
        step = stop / (*num - 1);

        for (int i = 0; i < *num; i++) {
            vec[i] = (int)round(pow(10, i * step));
        }
    }

    return vec;
}
/** perform the sum of a floating point array 
 * @param n (size_t) the number of vector components
 * @param v (double *) the floaring point array
 * @return (double) sum(v)
 */
extern double sum_vs(size_t n, double *v)
{
    double s = 0.;
    for (size_t i = 0; i < n; i++)
        s += *(v + i);
    return s;
}

/** perform the sum of squared components of an array
 * @param n (size_t) the number of vector components
 * @param v (double *) the array
 * @return (double) sum(v)
 */
extern double sum_vs_2(size_t n, double *v)
{
    double s = 0.;
    for (size_t i = 0; i < n; i++)
        s += v[i] * v[i];
    return s;
}

/**
 * @brief get the ReLU(x) (or softplus) of input unsigned 32-bit integer.
 * 
 * @param x (uint32_t) an integer
 * @return (int32_t) max(0, x).
 */
extern uint32_t softplus_u32(int32_t x)
{
    if(x > 0)
    {
        return x;
    }
    return 0;
}
/**
 * @brief Calculates the fraction of equal entries in two int8_t arrays.
 * 
 * The arrays are assumed to contain only +1 and -1 values. The function
 * counts the number of positions where the entries are equal and returns
 * the fraction of such positions with respect to the total size.
 * 
 * @param array1 Pointer to the first int8_t array.
 * @param array2 Pointer to the second int8_t array.
 * @param size The size of the arrays (number of elements in each array).
 * @return float The fraction of equal entries (value between 0 and 1).
 */
extern float calculateFractionEqual(int8_t *array1, int8_t *array2, size_t size) {
    if (size == 0) {
        return 0.0f; // Avoid division by zero
    }
    
    size_t equalCount = 0;
    for (size_t i = 0; i < size; i++) {
        equalCount += (array1[i] * array2[i]) == 1; // Increment if product is 1
    }
    
    return (float)equalCount / size;
}

/* STRING RELATED FUNCTIONS ************************************************* */
//
/**
 * @brief check if two strings are the same
 * 
 * This function compares two strings and returns true if they are identical,
 * false otherwise.
 * 
 * @param str1 (const char *) the first string to compare
 * @param str2 (const char *) the second string to compare
 * @return true if the two strings are the same false if the two strings differ 
 * for some characters
 */
extern bool strsme(const char *str1, const char *str2)
{
    if (strcmp(str1, str2))
        return false;
    else
        return true;
}
/**
 * @brief copy a number of characters from source string to destination string
 * 
 * This function copies a specified number of characters from the source string
 * to the destination string. It ensures that the destination string is null-terminated.
 * 
 * @param nc (size_t) number of characters to copy
 * @param outs (char *) the destination string
 * @param psrc (const char *) the source string
 */
extern void strncpy_safe(size_t nc, char *outs, const char* psrc)
{
    strncpy(outs, psrc, nc);
    outs[nc - 1] = '\0';
}
/**
 * @brief create a copy of a string after a character is met in string
 * 
 * This function creates a new string that is a copy of the portion of the source
 * string that follows the last occurrence of the specified character.
 * 
 * @param __strdst (char *) the destination string
 * @param __strsrc (const char *) the source string
 * @param __chrctr (const char) the caracter to lookup
 */
extern void strrmac(char *__strdst, const char *__strsrc, const char __chrctr)
{
    char *__rmstr = strrchr(__strsrc, __chrctr);
    strcpy(__strdst, __rmstr);
}
/**
 * @brief create a copy of a string until a character is met in string
 * 
 * This function creates a new string that is a copy of the portion of the source
 * string that precedes the last occurrence of the specified character.
 * 
 * @param __strdst (char *) the destination string
 * @param __strsrc (const char *) the source string
 * @param __chrctr (const char) the caracter to lookup
 */
extern void strrmuc(char *__strdst, const char *__strsrc, const char __chrctr)
{
    char *__rmstr = strrchr(__strsrc, __chrctr);
    strcpy(__strdst, __strsrc);
    __strdst[strlen(__strsrc)-strlen(__rmstr)] = '\0';
}
/**
 * @brief chack whether a string begin with another 
 * 
 * This function checks if the source string begins with the specified prefix string.
 * 
 * @param __strsrc (char *) the string where to look
 * @param __strbegin (const char *) the string to be contained
 * @return (bool) true if __strsrc begins with __strbegin
 * @return (bool) false otherwise
 */
extern bool strbws(char *__strsrc, const char *__strbegin)
{
    return (bool) (!(strncmp(__strsrc, __strbegin, strlen(__strbegin))));
}
/**
 * @brief get a string from file and check it is nor empty nor there are errors
 * with files
 * 
 * This function reads a line from the specified file and checks for errors. If
 * an error occurs, the program exits with a failure status.
 * 
 * @param fc (FILE **) file from which to read
 * @param row (char *) the char pointer to store the row content
 */
extern void __fgets_row(FILE **fc, char *row)
{
    if ((fgets(row, STRL1024, *fc) == NULL))
    {
        perror(MSG_ERR_FGETS);
        exit(EXIT_FAILURE);
    }
    row[strlen(row) - 1] = '\0';
}
/**
 * @brief Prepends one string to another.
 * 
 * This function takes a source string `t` and prepends it to the destination string `s`. 
 * It assumes that `s` has enough allocated space to hold the combined result.
 *
 * @param s (char *) The destination string to which `t` will be prepended. Must have enough space allocated.
 * @param t (const char *) The source string to prepend to `s`.
 */
void prepend(char *s, const char *t)
{
    size_t len = strlen(t);
    memmove(s + len, s, strlen(s) + 1);
    memcpy(s, t, len);
}
/**
 * @brief acquire size_t data from a string
 * 
 * This function scans the input string for a size_t number and returns it.
 * 
 * @param s (const char *) the input string
 * @return (size_t) the sscanned number
 */
extern size_t strtozu(const char *s)
{
    char c;
    int scanned;
    size_t i;
    scanned = sscanf(s, "%zu%c", &i, &c);
    if (scanned == 1)
        return i;
    else if (scanned > 1)
    {
        perror("strtozu");
        return i;
    }
    if (c != '\0' || errno != 0)
    {
        perror("strtozu");
        exit(EXIT_FAILURE);
    }
    return 0;
}

/**
 * @brief acquire uint32_t data from a string
 * 
 * This function scans the input string for a uint32_t number and returns it.
 * 
 * @param s (const char *) the input string
 * @return (uint32_t) the scanned number
 */
extern uint32_t strtou32(const char *s)
{
    char *endptr;
    uintmax_t tmp;
    tmp = strtoumax(s, &endptr, 10);
    if (*endptr != '\0')
    {
        fprintf(stderr, MSG_WRN_SCNU32 "%s\n", endptr);
        return (uint32_t) tmp;
    }
    if (tmp < UINT32_MAX)
    {
        return (uint32_t) tmp;
    }
    else
    {
        fprintf(stderr, MSG_ERR_SCNU32);
        exit(EXIT_FAILURE);
    }
}
/* FILES I/O FUNCTIONS ****************************************************** */
//
/**
 * @brief check that input file pointer points to an existing file
 * 
 * This function checks if a file exists at the specified path.
 * 
 * @param n (const char *) string containing path to file
 * @return true if file exists
 * @return false if file does not exist
 */
extern bool __feexist(const char *fn)
{
    FILE *f;
    if ((f = fopen(fn, "r")))
        fclose(f);
    else
        return false;
    return true;
}
/**
 * @brief check that input file pointer points to a non existing file
 * 
 * This function checks if a file does not exist at the specified path.
 * 
 * @param n (const char *) string containing path to file
 * @return true if file does not exist
 * @return false if file exists
 */
extern bool __fnexist(const char *fn)
{
    return (!(__feexist(fn)));
}
/**
 * @brief open a file according to an operative mode allowed by fopen
 * 
 * This function opens a file with the specified mode and checks for successful
 * opening. If the file cannot be opened, the program exits with a failure status.
 * 
 * @param f (FILE **) FILE pointer
 * @param fn (const char *) file name string
 * @param md (const char *) opening mode
 */
extern void __fopen(FILE **f, const char *fn, const char *md)
{
    if ((*f = fopen(fn, md)) == NULL)
    {
        perror(fn);
        exit(EXIT_FAILURE);
    }
}
/**
 * @brief check an fread worked out correctly
 * 
 * This function checks if the number of elements read by fread matches the
 * expected count. If not, an error message is printed and the program exits.
 * 
 * @param __frdval (size_t) the number of elements read off by fread
 * @param __frdcnt (size_t) the requested number of read
 */
extern void __fread_check(size_t __frdval, size_t __frdcnt)
{
     if (__frdval != __frdcnt)
     {
        perror(MSG_ERR_FREAD);
        exit(EXIT_FAILURE);
    }
}
/**
 * @brief open a pipe according to an operative mode allowed by fopen
 * 
 * This function opens a pipe with the specified mode and checks for successful
 * opening. If the pipe cannot be opened, the program exits with a failure status.
 * 
 * @param p (FILE **) pipe pointer
 * @param fn (const char *) file name string
 * @param md (const char *) opening mode
 */
extern void __popen(FILE **p, const char *fn, const char *md)
{
    if ((*p = popen(fn, md)) == NULL)
    {
        perror(MSG_ERR_POPEN);
        exit(EXIT_FAILURE);
    }
}

/**
 * @brief check if allocation made with malloc, calloc and realloc has worked,
 * othewise print on stderr the errno
 * 
 * This function checks if a pointer allocated with malloc, calloc, or realloc
 * is NULL, indicating a failed allocation. If the allocation failed, an error
 * message is printed and the program exits.
 * 
 * @param __ptr (void *) allocated pointer
 */
void __challoc(void *__ptr)
{
    if (__ptr == NULL)
    {
        perror("Alloc fail");
        exit(EXIT_FAILURE);
    }
}
/**
 * @brief (m)allocate array and check sucessfull allocation, return the
 * allocated pointer
 * 
 * This function allocates memory for an array of the specified size and checks
 * if the allocation was successful. It returns a pointer to the allocated memory.
 * 
 * @param n (size_t) the number of bytes to allocate
 * @return (void*) allocated pointer
 */
void *__chMalloc(size_t n)
{
    void *p = malloc(n);
    __challoc(p);
    return p;
}
/**
 * @brief (c)allocate array and check sucessfull allocation, return the
 * allocated pointer
 * 
 * This function allocates and zero-initializes memory for an array of the specified
 * size and checks if the allocation was successful. It returns a pointer to the allocated
 * memory.
 * 
 * @param n (size_t) the number of bytes to allocate
 * @return (void*) allocated pointer
 */
void *__chCalloc(size_t __nmemb, size_t __size)
{
    void *p = calloc(__nmemb, __size);
    __challoc(p);
    return p;
}
/**
 * @brief clone the current process, then, in child process, execve the program 
 * file to replace the child with the desired program.
 * 
 * This function creates a new process by forking the current process. In the child
 * process, it replaces the current program with a new program specified by the file
 * path. If the process creation or program execution fails, an error message is printed
 * and the program exits.
 * 
 * @param argv 
 * @return pid (pid_t) 
 */
extern pid_t call(char* argv[]) {
    pid_t pid = fork();
    if (pid == 0)
    {
        char* envp[] = { NULL };
        
        execve(argv[0], argv, envp);
        perror("Error execve");
        exit(EXIT_FAILURE);
    }
    else
    {
        
        return pid;
    }
}
/**
 * @brief wait for all the children of the program to end
 * 
 * This function waits for all child processes of the current process to terminate.
 * It prints the termination status of each child process.
 */
extern void __wait_childs(void)
{
    int crps;
    int status;
    while ((crps = wait(&status)) > 0)
    {
        printf("%d: PID %d exited with status 0x%.4X\n",
               (int)getpid(), crps, status);
    }
}


/**
 * @brief generate random string
 * 
 * This function generates a random string of the specified size using characters
 * from a predefined set. The generated string is null-terminated.
 * 
 * @param str 
 * @return size
 */
char *rand_string(char *str, size_t size)
{
    const char charset[] = "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKILMNOPQRSTUVWXYZ";
    if (size) {
        --size;
        for (size_t n = 0; n < size; n++) {
            uint64_t key = RNG_u64() % (sizeof charset - 1);
            str[n] = charset[key];
        }
        str[size] = '\0';
    }
    return str;
}


extern void __make_adj_from_tmp(size_t i, size_t j, double tmp, double_p **adj)
{
    *(*(*adj + j) + i) = tmp;
    *(*(*adj + i) + j) = *(*(*adj + j) + i);
}

extern void __fill_adj__(FILE **f, size_t N, double_p **adj)
{
    double tmp;
    for (size_t i = 0; i < N; i++)
    {
        for (size_t j = i; j < N; j++)
        {
            __fread_check(fread(&tmp, sizeof(tmp), 1, *f), 1);
            __make_adj_from_tmp(i, j, tmp, adj);
            printf("%.4g ", *(*(*adj + i) + j));
        }
        printf("\n");
    }
}

extern void __fill_edgl_read__(FILE **f, size_t N, double_p **edgl, 
    size_tp **neighs, size_tp *neigh_len)
{
    size_t node_i, cntr, last = 0;
    double w_ij;
    for (size_t i = 0; i < N; i++)
    {
        *(*edgl + i) = __chMalloc(1 * sizeof(**edgl));
        *(*neighs + i) = __chMalloc(1 * sizeof(**neighs));
    }
    cntr = 0;
    node_i = 0;
    #ifdef BIN_EDGE_LIST_FORMAT
    size_t i, j;
    while (fread(&i, sizeof(size_t), 1, *f) == 1 &&
        fread(&j, sizeof(size_t), 1, *f) == 1 &&
        fread(&w_ij, sizeof(double), 1, *f) == 1)
    {
    #else
    for(size_t i, j; fscanf(*f, "%zu %zu %lf", &i, &j, &w_ij) != EOF;)
    {
    #endif
        last = i;
        if (i != node_i)
        {
            *(*neigh_len + i-1) = cntr;
            node_i++;
            // printf("%zu, %zu, %zu\n", i-1, *(*neigh_len + i-1), cntr);
            cntr = 0;
        }
        *(*(*neighs + i) + cntr) = j;
        *(*(*edgl + i) + cntr) = w_ij;
        *(*edgl + i) = realloc(*(*edgl + i), (++cntr + 1) * sizeof(**edgl));
        *(*neighs + i) = realloc(*(*neighs + i), (cntr + 1) * sizeof(**neighs));
    }
    *(*neigh_len + last) = cntr;
}

extern void __fill_edgl_make__(FILE **f, size_t N, double_p **adj, double_p **edgl, size_tp **neighs, size_tp *neigh_len)
{
    size_t cntr;
    for (size_t i = 0; i < N; i++)
    {
        cntr = 0;
        *(*edgl + i) = __chMalloc(1 * sizeof(**edgl));
        *(*neighs + i) = __chMalloc(1 * sizeof(**neighs));
        for (size_t j = 0; j < N; j++)
        {
            if (fabs(*(*(*adj + j) + i)) > 0.)
            {
                *(*(*neighs + i) + cntr) = j;
                *(*(*edgl + i) + cntr) = *(*(*adj + j) + i);
                *(*edgl + i) = realloc(*(*edgl + i), (++cntr + 1) * sizeof(**edgl));
                *(*neighs + i) = realloc(*(*neighs + i), (cntr + 1) * sizeof(**neighs));
                fprintf(*f, "%zu %zu %lf\n", i, j, *(*(*adj + j) + i));
            }
        }
        *(*neigh_len + i) = cntr;
        // printf("%zu, %zu, %zu\n", i, *(*neigh_len + i), cntr);
    }
}



// Function to read edges from a binary file
extern Edges __read_bin_EdgeList__(const char *filename, size_t *edge_count) {
    FILE *file;
    Edge *edges;
    //
    __fopen(&file, filename, "rb");
    fseek(file, 0, SEEK_END);
    long file_size = ftell(file);
    fseek(file, 0, SEEK_SET);
    *edge_count = (size_t) file_size / sizeof(Edge);
    
    edges = __chMalloc(file_size);

    fread(edges, sizeof(Edge), *edge_count, file);
    fclose(file);

    return edges;
}

/**
 * @brief Processes edges from a binary file and organizes them into node-edge structures.
 *
 * This function reads edges from a binary file, calculates the number of neighbors
 * for each node, and organizes the edges into a node-edge structure for efficient access.
 *
 * @param filename (const char*) The name of the binary file containing the edge list.
 * @param N (size_t) The number of nodes in the graph.
 * @param edges (Edges*) Pointer to the structure where the edges will be stored.
 * @param node_edges (NodesEdges*) Pointer to the structure where node-edge relationships will be stored.
 * @param neigh_len (size_tp*) Pointer to an array where the number of neighbors for each node will be stored.
 */
extern void process_edges(const char *filename, size_t N, Edges *edges, NodesEdges *node_edges, size_tp *neigh_len) {
    size_t edge_count;
    *edges = __read_bin_EdgeList__(filename, &edge_count);
    *node_edges = __chMalloc(N * sizeof(**node_edges));
    *neigh_len = __chCalloc(N, sizeof(**neigh_len));
    for (Edge *e = *edges; e < *edges + edge_count; e++) {
        (*neigh_len)[e->u]++;
        (*neigh_len)[e->v]++;
    }
    // Allocate the actual arrays now that we know the sizes
    for (size_t i = 0; i < N; i++) {
        if ((*neigh_len)[i] > 0) {
            (*node_edges)[i].neighbors = __chMalloc((*neigh_len)[i] * sizeof(*(*node_edges)[i].neighbors));
            (*node_edges)[i].weights = __chMalloc((*neigh_len)[i] * sizeof(*(*node_edges)[i].weights));
        }
        (*neigh_len)[i] = 0; // Reset to use as an index
    }
    for (Edge *e = *edges; e < *edges + edge_count; e++) {
        uint64_t u = e->u;
        uint64_t v = e->v;
        double w = e->w;

        (*node_edges)[u].neighbors[(*neigh_len)[u]] = v;
        (*node_edges)[u].weights[(*neigh_len)[u]] = w;
        (*neigh_len)[u]++;

        (*node_edges)[v].neighbors[(*neigh_len)[v]] = u;
        (*node_edges)[v].weights[(*neigh_len)[v]] = w;
        (*neigh_len)[v]++;
    }
}

/**
 * @brief Builds a string identifier by appending an underscore and the input string.
 *
 * This function takes an input string and appends it to an underscore to create
 * a string identifier. If the input string is empty or NULL, the output string
 * will be empty.
 *
 * @param _str_id (const char*) The input string to append to the identifier.
 * @param str_id (char*) The output buffer where the identifier will be stored.
 * @param str_sz (size_t) The size of the output buffer.
 */
extern void build_str_id(const char *_str_id, char *str_id, size_t str_sz) {
    if (_str_id && _str_id[0])
        snprintf(str_id, str_sz, "_%s", _str_id);
    else if (str_sz)
        str_id[0] = '\0';
}

/* DICTIONARY IMPLEMENTATION ************************************************ */
// //
// #define HASHSIZE 101
// static struct nlist *hashtab[HASHSIZE]; /* pointer table */

// struct nlist { /* table entry: */
//     struct nlist *next; /* next entry in chain */
//     char *name; /* defined name */
//     char *defn; /* replacement text */
// };


// /* hash: form hash value for string s */
// unsigned hash(char *s)
// {
//     unsigned hashval;
//     for (hashval = 0; *s != '\0'; s++)
//       hashval = *s + 31 * hashval;
//     return hashval % HASHSIZE;
// }

// /* lookup: look for s in hashtab */
// struct nlist *lookup(char *s)
// {
//     struct nlist *np;
//     for (np = hashtab[hash(s)]; np != NULL; np = np->next)
//         if (strcmp(s, np->name) == 0)
//           return np; /* found */
//     return NULL; /* not found */
// }

// char *strdup(char *);
// /* install: put (name, defn) in hashtab */
// struct nlist *install(char *name, char *defn)
// {
//     struct nlist *np;
//     unsigned hashval;
//     if ((np = lookup(name)) == NULL) { /* not found */
//         np = (struct nlist *) malloc(sizeof(*np));
//         if (np == NULL || (np->name = strdup(name)) == NULL)
//           return NULL;
//         hashval = hash(name);
//         np->next = hashtab[hashval];
//         hashtab[hashval] = np;
//     } else /* already there */
//         free((void *) np->defn); /*free previous defn */
//     if ((np->defn = strdup(defn)) == NULL)
//        return NULL;
//     return np;
// }

// char *strdup(char *s) /* make a duplicate of s */
// {
//     char *p;
//     p = (char *) malloc(strlen(s)+1); /* +1 for ’\0’ */
//     if (p != NULL)
//        strcpy(p, s);
//     return p;
// }

/**
 * @brief Calculates the normalized Ising energy correctly (range -1 to +1) - OPTIMIZED.
 *
 * This function computes the Ising energy normalized to the range [-1, +1] where:
 * - Energy = -1 when all neighboring spins are aligned (ferromagnetic ground state)
 * - Energy = +1 when all neighboring spins are anti-aligned (antiferromagnetic state)
 * 
 * The energy is calculated as: E = -sum(s_i * w_ij * s_j) / total_edges
 * where the sum is over all unique edges (avoiding double counting).
 *
 * OPTIMIZATIONS:
 * - Pre-calculates total edges in first pass to avoid conditionals in main loop
 * - Uses const pointers and cached values for faster memory access
 * - Minimizes array indexing operations
 * - Compiler-friendly loop structure for better vectorization
 *
 * @param N (size_t) The number of spins in the system.
 * @param s (spin_tp) Pointer to the array of spins.
 * @param nlen (size_tp) Array containing the number of neighbors for each spin.
 * @param ne (NodesEdges) Array of neighbor data structures.
 *
 * @return (double) The normalized Ising energy in range [-1, +1].
 */
double calc_ising_energy_normalized(size_t N, spin_tp s, size_tp nlen, NodesEdges ne) {
    // First pass: count total edges (done once, outside main loop)
    size_t total_edges = 0;
    for (size_t i = 0; i < N; i++) {
        total_edges += nlen[i];
    }
    total_edges >>= 1; // Divide by 2 since each edge is counted twice
    
    // Handle edge case
    if (total_edges == 0) {
        return 0.0;
    }
    
    // Second pass: calculate energy sum (optimized inner loop)
    double sum = 0.;
    const int8_t * const spin_arr = s; // Cache pointer for faster access
    
    for (size_t i = 0; i < N; i++) {
        const int8_t spin_i = spin_arr[i]; // Cache current spin
        const size_t n_neighbors = nlen[i];
        const size_t * const neighbors = ne[i].neighbors;
        const double * const weights = ne[i].weights;
        
        // Unroll-friendly inner loop with pointer arithmetic
        for (size_t j = 0; j < n_neighbors; j++) {
            const size_t neighbor_idx = neighbors[j];
            // Only count each edge once
            if (i < neighbor_idx) {
                sum += weights[j] * spin_i * spin_arr[neighbor_idx];
            }
        }
    }
    
    // Normalize: -sum/total_edges gives energy per edge
    return -sum / total_edges;
}
/**
 * @brief Fast calculation of normalized Ising energy - matches calc_totEnergy performance.
 *
 * This function computes the Ising energy with the same performance as calc_totEnergy
 * but with correct normalization and no double counting. It uses the same loop structure
 * as calc_totEnergy but divides by the actual number of edges instead of 2*N.
 *
 * PERFORMANCE OPTIMIZATIONS:
 * - Single pass through all nodes (same as calc_totEnergy)
 * - No conditionals in inner loop
 * - Direct memory access patterns
 * - Minimal function calls
 *
 * @param N (size_t) The number of spins in the system.
 * @param s (spin_tp) Pointer to the array of spins.
 * @param nlen (size_tp) Array containing the number of neighbors for each spin.
 * @param ne (NodesEdges) Array of neighbor data structures.
 *
 * @return (double) The normalized Ising energy in range [-1, +1].
 */
double calc_ising_energy_fast(size_t N, spin_tp s, size_tp nlen, NodesEdges ne) {
    double sum = 0.;
    size_t total_edges = 0;
    
    // Single pass: same structure as calc_totEnergy but count edges
    for (size_t i = 0; i < N; i++) {
        double tmp = s[i] * neighWeight_magn(ne[i], s, nlen[i]);
        sum += tmp;
        total_edges += nlen[i];
    }
    
    // Correct normalization: divide by actual edge count, not by 2*N
    total_edges >>= 1; // Each edge counted twice
    return (total_edges > 0) ? -sum / (2 * total_edges) : 0.0;
}