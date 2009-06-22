#ifndef _MBACLOCK_H_
#define _MBACLOCK_H_

#include <time.h>
class MBAclock {
  //long    i = 600000L;
  clock_t start, finish;
  double  duration;
  
public:
  MBAclock(){start = clock();}
  ~MBAclock(){};
  double getInterval() {
    finish = clock(); 
    duration = (double)(finish - start) / CLOCKS_PER_SEC;
    start = finish;
    return duration;
  }
};

#endif
