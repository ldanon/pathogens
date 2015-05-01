//typedef enum {ERROR = -1, FALSE, TRUE} LOGICAL;

#define BitSet(arg,posn) ((arg) | (1L << (posn)))
#define BitClr(arg,posn) ((arg) & ~(1L << (posn)))
#define BitFlp(arg,posn) ((arg) ^ (1L << (posn)))
#define BitTst(arg,posn) (((arg) & (1L << (posn))) > 0)

#define BOOL(x) (!(!(x)))






