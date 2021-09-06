#ifndef THREADED_LIST_HPP_
#define THREADED_LIST_HPP_

#include <list>
#include <vector>

namespace cartdg
{

template <typename T>
class Threaded_list
{
  public:
  typedef typename std::list<T>::iterator iterator;
  typedef typename std::vector<iterator> it_vec;

  protected:
  std::list<T> data;
  int nt;
  std::vector<int> sizes;
  std::vector<iterator> interior_iterators;

  public:
  Threaded_list(int n_thread) {}
  Threaded_list() {}
  int size();
  int n_thread();
  iterator begin();
  iterator end();
  it_vec thread_begins() {it_vec empty; return empty;}
  it_vec thread_ends() {it_vec empty; return empty;}
  void emplace(int i_thread, ...);
  void erase(iterator position);
  void balance();
};

}
#endif
