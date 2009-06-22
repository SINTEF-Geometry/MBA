#ifndef _MYCLOCK_H_
#define _MYCLOCK_H_

#include <time.h>
class MyClock {
  //long    i = 600000L;
  clock_t start, finish;
  double  duration;
  
public:
  MyClock(){start = clock();}
  ~MyClock(){};
  double getInterval() {
    finish = clock(); 
    duration = (double)(finish - start) / CLOCKS_PER_SEC;
    start = finish;
    return duration;
  }
};

#endif
