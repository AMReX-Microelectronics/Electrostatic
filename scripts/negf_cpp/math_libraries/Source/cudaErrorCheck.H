////////////////////////////////////////////////////////////////////////////////
// CUDA error checking
////////////////////////////////////////////////////////////////////////////////
#define checkCudaErrors(val) \
    __checkCudaErrors__((val), #val, __FILE__, __LINE__)

template <typename T>
inline void __checkCudaErrors__(T code, const char *func, const char *file,
                                int line)
{
    if (code)
    {
        fprintf(stderr, "CUDA error at %s:%d code=%d \"%s\" \n", file, line,
                (unsigned int)code, func);
        cudaDeviceReset();
        exit(EXIT_FAILURE);
    }
}
