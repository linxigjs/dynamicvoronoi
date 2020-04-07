#include "bucketedqueue.h"

#include "limits.h"
#include <stdio.h>
#include <stdlib.h>

template <class T>
BucketPrioQueue<T>::BucketPrioQueue() {
  clear();
}

template <class T>
bool BucketPrioQueue<T>::empty() {
  return (count==0);
}


template <class T>
void BucketPrioQueue<T>::push(int prio, T t) {
  buckets[prio].push(t);
  //如果buckets是空，或者prio是最小（最高）优先级，那prio就是第一个要pop的元素
  if (nextPop == buckets.end() || prio < nextPop->first) {
    nextPop = buckets.find(prio);
  }
  count++;
}

template <class T>
T BucketPrioQueue<T>::pop() {
  while (nextPop!=buckets.end() && nextPop->second.empty()) ++nextPop;

  T p = nextPop->second.front();
  nextPop->second.pop();
  if (nextPop->second.empty()) {
    typename BucketType::iterator it = nextPop;
    nextPop++;
    buckets.erase(it);
  }
  count--;
  return p;
}

