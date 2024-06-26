#if !defined(OFT_INCLUDE_GUARD)
#define OFT_INCLUDE_GUARD

#if (!defined OFT_STACK && defined OFT_PROFILE)
#undef OFT_PROFILE
#endif

#if (defined OFT_STACK && defined DEBUG_STACK_MOD_IND)
#define STACK_SUB_IND 1
#define DEBUG_STACK_PUSH CALL oft_stack_push(DEBUG_STACK_MOD_IND,DEBUG_STACK_SUB_IND)
#define DEBUG_STACK_POP CALL oft_stack_pop
#else
#define DEBUG_STACK_PUSH
#define DEBUG_STACK_POP
#endif

#if !defined( OFT_OMP_VTHRESH )
#define OFT_OMP_VTHRESH 2000
#endif

#if !defined(HDF5_TLEN)
#define HDF5_TLEN 4
#endif

#if !defined(OFT_PATH_SLEN)
#define OFT_PATH_SLEN 80
#endif

#if !defined(OFT_HIST_SLEN)
#define OFT_HIST_SLEN 40
#endif

#endif